"""Prepare and run mesh/time qualification and declared input scenarios."""
from pathlib import Path
from collections import deque
import itertools,json,re,shutil,subprocess,sys,os
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import prepare_model as pm
ROOT=pm.ROOT
ANSYS=os.environ.get('ANSYS_EXECUTABLE',r'C:\Program Files\ANSYS Inc\v241\ansys\bin\winx64\ANSYS241.exe')

def seed_state(source):
    if (source/'final_temperatures.npy').exists():return np.load(source/'final_temperatures.npy')
    n=len(np.load(source/'mesh.npz')['nodes'])
    with (source/'nodal_temperatures.csv').open() as f:tail=list(deque(f,maxlen=n))
    rows=np.array([[float(v) for v in line.split(',')] for line in tail]);assert np.array_equal(rows[:,1],np.arange(1,n+1))
    assert abs(rows[0,0]%5580)<1e-6
    return rows[:,2]

def interpolate_field(source,target_mesh):
    src=np.load(source/'mesh.npz');field=seed_state(source);nodes=src['nodes'];conn=src['connectivity'];tags=src['part_names']
    axes=[np.unique(nodes[:,j]) for j in range(3)];box=nodes[conn];low=box.min(axis=1);high=box.max(axis=1)
    cells={tuple(int(np.searchsorted(axes[j],p[j])) for j in range(3)):i for i,p in enumerate(low)}
    targetnodes=target_mesh['nodes'];targetconn=target_mesh['connectivity'];targettags=target_mesh['part_names'];incident=[set() for _ in targetnodes]
    for row,tag in zip(targetconn,targettags):
        for n in row:incident[int(n)].add(str(tag))
    corner=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    out=np.zeros(len(targetnodes))
    for nid,p in enumerate(targetnodes):
        candidates=[]
        for axis,value in zip(axes,p):
            k=int(np.searchsorted(axis,value));near=[q for q in [k-1,k] if 0<=q<len(axis) and abs(axis[q]-value)<1e-12]
            if near:
                q=near[0];candidates.append([i for i in [q-1,q] if 0<=i<len(axis)-1])
            else:candidates.append([max(0,min(len(axis)-2,k-1))])
        samples=[]
        for key in itertools.product(*candidates):
            e=cells.get(key)
            if e is None or str(tags[e]) not in incident[nid]:continue
            q=np.clip((p-low[e])/(high[e]-low[e]),0,1)
            weight=np.prod(np.where(corner,q,1-q),axis=1);samples.append(float(weight@field[conn[e]]))
        assert samples,(nid,p,incident[nid])
        out[nid]=np.mean(samples)
    return out

def configure(run,initial,dt=10,cycles=3,solver=2,changes=None):
    changes=changes or {};duration=cycles*5580;settling=(cycles-2)*5580
    assert cycles>=2 and dt>0 and 5580/dt==int(5580/dt)
    np.save(run/'initial_temperatures.npy',initial)
    text=(run/'model.inp').read_text(encoding='utf-8')
    # Replace any existing seed, retaining all topology and surface definitions.
    text=re.sub(r'^IC,.*\n','',text,flags=re.M);text=text.replace('TUNIF,300','')
    text=text.replace('SPCTEMP,2,2.7','\n'.join(f'IC,{i},TEMP,{t:.12g}' for i,t in enumerate(initial,1))+'\nSPCTEMP,2,2.7')
    text=re.sub(r'^RADOPT,.*$',f'RADOPT,,1E-4,{solver},1000,1E-5,1',text,flags=re.M)
    if 'hemicube_resolution' in changes:
        text=re.sub(r'^HEMIOPT,.*$',f'HEMIOPT,{int(changes["hemicube_resolution"])},1E-6',text,flags=re.M)
    text=re.sub(r'^RESCONTROL,.*\n','',text,flags=re.M)
    # Set the last solution portion explicitly. Two final cycles are recorded.
    first=text.index('AUTOTS,ON');tail=text[first:];tail=tail[:tail.index('FINISH')]
    solution=['AUTOTS,ON',f'DELTIM,{dt},{min(dt,1)},{dt}','OUTRES,ALL,NONE','NEQIT,100','SAVE']
    if (run/'cached_view_factors.vf').exists():solution.append('VFOPT,READ,cached_view_factors,vf,,BINA')
    if settling:solution.extend([f'TIME,{settling}','SOLVE'])
    solution.extend(['OUTRES,NSOL,ALL',f'TIME,{duration}','SOLVE'])
    text=text[:first]+'\n'.join(solution)+'\nFINISH\n/POST1\nSET,LAST\nSAVE\nFINISH\n/EXIT,NOSAVE\n'
    # Complete periodic heating tables cover the entire simulation.
    flux=np.loadtxt(run/'applied_flux.csv',delimiter=',',skiprows=1)
    if 'alpha' in changes or 'external_epsilon' in changes:
        alpha=changes.get('alpha',.77);eps=changes.get('external_epsilon',.72)
        assert 0<=alpha<=1 and 0<=eps<=1
        flux[:,20:26]=alpha*(flux[:,2:8]+flux[:,8:14])+eps*flux[:,14:20]
        np.savetxt(run/'applied_flux.csv',flux,delimiter=',',header=','.join(['time_s','eclipse']+[f'{k}_panel{i}_W_m2' for k in ['solar','albedo','earth_ir','absorbed'] for i in range(1,7)]),comments='')
    ext={(int(e),int(f)) for e,f,side,area in np.load(run/'mesh.npz')['faces'] if side>0}
    def emission(match):
        e,f,v=int(match[1]),int(match[2]),float(match[3]);key='external_epsilon' if (e,f) in ext else 'internal_epsilon'
        return f'SFE,{e},{f},RDSF,1,{changes.get(key,v):.12g}'
    text=re.sub(r'^SFE,(\d+),(\d+),RDSF,1,([^\n]+)$',emission,text,flags=re.M)
    for prop,mid,value in [('KXX',1,changes.get('panel_k',1.03)),('KXX',3,changes.get('pcb_k',.64))]:
        text=re.sub(rf'^MP,{prop},{mid},.*$',f'MP,{prop},{mid},{value:.12g}',text,flags=re.M)
    scale=changes.get('capacity_scale',1.0);assert scale>0
    text=re.sub(r'^MP,C,(\d+),([^\n]+)$',lambda match:f'MP,C,{match[1]},{pm.PROPS[int(match[1])][1]*scale:.12g}',text,flags=re.M)
    text=re.sub(r'^(?:\*DIM,Q[1-6],TABLE.*|Q[1-6]\(.*)\n','',text,flags=re.M)
    extended=np.vstack([np.c_[flux[:-1,0]+cycle*5580,flux[:-1,1:]] for cycle in range(cycles)]+[np.r_[duration,flux[0,1:]][None,:]])
    tables=[]
    for side in range(1,7):
        tables.append(f'*DIM,Q{side},TABLE,{len(extended)},1,1,TIME')
        for j,row in enumerate(extended,1):tables.extend([f'Q{side}({j},0)={row[0]:.12g}',f'Q{side}({j},1)={row[19+side]:.12g}'])
    text=text.replace('SFE,','\n'.join(tables)+'\nSFE,',1)
    (run/'model.inp').write_text(text,encoding='utf-8')
    (run/'study_provenance.json').write_bytes((json.dumps({'dt_s':dt,'cycles':cycles,'settling_cycles':cycles-2,'duration_s':duration,'recorded_cycles':2,'radiosity_solver':solver,'changes':changes,'input_uncertainty_note':'Scenario ranges are declared engineering assumptions, not measured probability bounds.'},indent=2)+'\n').encode('utf-8'))

def mesh_case(h,features,source,cycles=3):
    generated=pm.prepare(h,.5,16740,'orbital',0,features)
    dest=ROOT/'data/ansys'/f'convergence_h{h:g}_features{features}';shutil.copytree(generated,dest)
    target=np.load(dest/'mesh.npz');initial=interpolate_field(source,target)
    configure(dest,initial,cycles=cycles)
    p=json.loads((dest/'study_provenance.json').read_text());p.update(mesh_h_mm=h,thin_region_divisions=features,seed_source=source.name,seed_method='Piecewise trilinear interpolation within the same occupied part; no interpolation across vacuum')
    (dest/'study_provenance.json').write_bytes((json.dumps(p,indent=2)+'\n').encode('utf-8'))
    return dest

def clone_case(label,source,dt=10,cycles=4,changes=None,use_cache=False):
    dest=ROOT/'data/ansys'/label;dest.mkdir(exist_ok=False)
    for name in ['mesh.npz','mesh_audit.json','applied_flux.csv','model.inp','post.inp']:shutil.copy2(source/name,dest/name)
    # Cache reuse is opt-in and requires a completed equivalence check.
    if use_cache and (source/'orbital.vf').exists():shutil.copy2(source/'orbital.vf',dest/'cached_view_factors.vf')
    initial=seed_state(source)
    # Initial-condition acceleration by a physical mean-radiation balance;
    # this changes only the starting iterate, never the solved parameters.
    if changes and ('alpha' in changes or 'external_epsilon' in changes):
        flux=np.loadtxt(source/'applied_flux.csv',delimiter=',',skiprows=1);new=changes.get('alpha',.77)*(flux[:,2:8]+flux[:,8:14])+changes.get('external_epsilon',.72)*flux[:,14:20]
        prior=json.loads((source/'study_provenance.json').read_text()) if (source/'study_provenance.json').exists() else {}
        prior_epsilon=prior.get('changes',{}).get('external_epsilon',.72)
        ratio=new.sum()/flux[:,20:26].sum()*prior_epsilon/changes.get('external_epsilon',.72)
        initial=initial*ratio**.25
    configure(dest,initial,dt,cycles,changes=changes)
    p=json.loads((dest/'study_provenance.json').read_text());p.update(seed_source=source.name,seed_method='Stored phase-zero field; optical scenarios use declared mean-radiation balance scaling for startup only')
    (dest/'study_provenance.json').write_bytes((json.dumps(p,indent=2)+'\n').encode('utf-8'))
    return dest

def solve(case):
    assert not (case/'orbital.rth').exists(),'Existing results protected'
    nodes=len(np.load(case/'mesh.npz')['nodes']);memory=384 if nodes<8000 else (768 if nodes<16000 else 1280)
    for job,inp,out in [('orbital','model.inp','solver.out'),('export','post.inp','post.out')]:
        subprocess.run([ANSYS,'-b','-smp','-np','2','-m',str(memory),'-db','64','-j',job,'-i',inp,'-o',out],cwd=case,check=True)
        counts=re.findall(r'NUMBER OF ERROR\s+MESSAGES ENCOUNTERED\s*=\s*(\d+)',(case/out).read_text(errors='replace'));assert counts and int(counts[-1])==0,(case,out)
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1'
    subprocess.run([sys.executable,str(ROOT/'analyze_study.py'),str(case)],check=True,env=env)

if __name__=='__main__':
    source=ROOT/'data/ansys/interface_bonded_control'
    if not source.exists():source=ROOT/'data/ansys/captured_convergence_h15_features1'
    medium=mesh_case(10,2,source);solve(medium)
    fine=mesh_case(7.5,3,medium);solve(fine)
