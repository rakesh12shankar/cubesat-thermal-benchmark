"""Controlled zero-conductance interface limits; geometry and radiation unchanged."""
from pathlib import Path
from collections import defaultdict, deque
import json, re, shutil, subprocess, sys,os,argparse
import numpy as np
from prepare_model import ROOT

ANSYS=os.environ.get('ANSYS_EXECUTABLE',r'C:\Program Files\ANSYS Inc\v241\ansys\bin\winx64\ANSYS241.exe')

def prepare(label):
    base=ROOT/'data/ansys/orbital_h15_e0.5_balanced_strict'
    if not base.exists():base=ROOT/'data/ansys/captured_interface_bonded_control'
    dest=ROOT/'data/ansys'/f'interface_{label}'
    dest.mkdir(exist_ok=False)
    m=np.load(base/'mesh.npz'); old=m['nodes']; tags=m['part_names']; conn=m['connectivity']
    if not base.name.startswith('captured_') and (base/'nodal_temperatures.csv').exists():
        with (base/'nodal_temperatures.csv').open() as f: tail=list(deque(f,maxlen=len(old)))
        seed=np.array([[float(v) for v in line.split(',')] for line in tail])
        assert np.all(seed[:,0]==11160) and np.all(seed[:,1]==np.arange(1,len(old)+1))
    else:
        state=np.load(base/'initial_temperatures.npy')
        seed=np.column_stack([np.full(len(old),11160),np.arange(1,len(old)+1),state])
    # Split only nodes affected by the selected interface. At each original
    # node, retain all other pairwise perfect connections, including legitimate
    # alternate paths through frame/bolts at junctions. Connected components
    # implement transitive perfect continuity without introducing new paths.
    incident=defaultdict(set)
    for row,tag in zip(conn,tags):
        for n in row: incident[int(n)].add(str(tag))
    def blocked(a,b):
        seam=a.startswith('panel') and b.startswith('panel')
        pcb=(a=='frame' and b.startswith('pcb')) or (b=='frame' and a.startswith('pcb'))
        return (seam and label in ['seams_off','both_off']) or (pcb and label in ['pcb_frame_off','both_off'])
    mapping={}; xyz=[]; initial=[]; split_count=0; retained_junctions=0
    for n in range(len(old)):
        remaining=set(incident[n]); groups=[]
        while remaining:
            group={min(remaining)};remaining-=group
            while True:
                add={b for b in remaining if any(not blocked(a,b) for a in group)}
                if not add:break
                group|=add;remaining-=add
            groups.append(group)
        split_count+=len(groups)>1
        for group in groups:
            if any(blocked(a,b) for a in group for b in group if a!=b):retained_junctions+=1
            new=len(xyz);xyz.append(old[n]);initial.append(seed[n,2])
            for tag in group:mapping[n,tag]=new
    newconn=np.array([[mapping[int(n),str(tag)] for n in row] for row,tag in zip(conn,tags)])
    xyz=np.array(xyz);initial=np.array(initial)
    # Element shapes, volumes, loads, and radiating faces remain identical.
    assert np.array_equal(xyz[newconn],old[conn])
    np.savez_compressed(dest/'mesh.npz',**{k:(xyz if k=='nodes' else newconn if k=='connectivity' else m[k]) for k in m.files})
    np.save(dest/'initial_temperatures.npy',initial)
    for name in ['mesh_audit.json','applied_flux.csv']:shutil.copy2(base/name,dest/name)
    text=(base/'model.inp').read_text(encoding='utf-8')
    text=re.sub(r'^N,.*\n','',text,flags=re.M)
    node_lines='\n'.join(f'N,{i},{p[0]:.12g},{p[1]:.12g},{p[2]:.12g}' for i,p in enumerate(xyz,1))+'\n'
    text=text.replace('MP,KXX,4,23\n','MP,KXX,4,23\n'+node_lines)
    element_lines={eid:'EN,'+str(eid)+','+','.join(str(n+1) for n in row) for eid,row in enumerate(newconn,1)}
    text=re.sub(r'^EN,(\d+),.*$',lambda match:element_lines[int(match[1])],text,flags=re.M)
    text=re.sub(r'^IC,.*\n','',text,flags=re.M)
    text=text.replace('SPCTEMP,2,2.7','\n'.join(f'IC,{i},TEMP,{t:.12g}' for i,t in enumerate(initial,1))+'\nSPCTEMP,2,2.7')
    (dest/'model.inp').write_text(text,encoding='utf-8')
    post=(base/'post.inp').read_text(encoding='utf-8').replace(str(len(old)),str(len(xyz)))
    (dest/'post.inp').write_text(post,encoding='utf-8')
    provenance={'label':label,'duration_s':11160,'dt_max_s':10,'source_phase_zero':'balanced_strict final time 11160 s',
      'split_original_nodes':split_count,'nodes':len(xyz),'retained_alternate_path_junctions':retained_junctions,
      'interpretation':'Zero direct conductance limiting diagnostic, not recovered author interfaces. Alternate connections at junctions retained; surfaces are coincident and sealed, so no new radiation is introduced.',
      'control':'Same element geometry, properties, exterior loads, radiation facets, numerical settings and mapped initial field. Bonded control uses identical original connectivity.'}
    (dest/'interface_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    audit=json.loads((dest/'mesh_audit.json').read_text(encoding='utf-8'))
    audit['nodes']=len(xyz)
    audit.pop('shared_node_part_connections',None)
    audit['interface_provenance']=provenance
    audit['status']='Element geometry preserved; interface nodes partitioned as documented. Solution not yet audited.'
    (dest/'mesh_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(provenance),flush=True)
    return dest,initial

def run(dest,initial):
    for job,inp,out in [('orbital','model.inp','solver.out'),('export','post.inp','post.out')]:
        subprocess.run([ANSYS,'-b','-smp','-np','2','-j',job,'-i',inp,'-o',out],cwd=dest,check=True)
        counts=re.findall(r'NUMBER OF ERROR\s+MESSAGES ENCOUNTERED\s*=\s*(\d+)',(dest/out).read_text(errors='replace'))
        assert counts and int(counts[-1])==0,(dest,out)
    csv=dest/'nodal_temperatures.csv'; csv.rename(dest/'nodal_temperatures_exported.csv')
    with csv.open('w',encoding='utf-8') as f:
        for i,t in enumerate(initial,1):f.write(f'0,{i},{t:.12g}\n')
        with (dest/'nodal_temperatures_exported.csv').open() as source:shutil.copyfileobj(source,f)
    subprocess.run([sys.executable,str(ROOT/'analyze_results.py'),str(dest)],check=True,stdout=subprocess.DEVNULL)
    a=json.loads((dest/'solution_audit.json').read_text())
    print(json.dumps({'case':dest.name,'periodicity_K':a['last_orbit_periodicity_max_node_K'],'energy_rms_W':a['last_orbit_rms_energy_residual_W']}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--ansys',default=ANSYS)
    ANSYS=parser.parse_args().ansys
    for label in ['bonded_control','seams_off','pcb_frame_off','both_off']:
        dest,initial=prepare(label);run(dest,initial)
