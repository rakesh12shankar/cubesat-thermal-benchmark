"""Stream MAPDL nodal histories; retain one orbit for exact phase comparisons.

Uses actual case material capacities and exterior emissivities, including
perturbed inputs. Does not hold every element's entire history in memory.
"""
from pathlib import Path
import argparse,csv,json,re
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
ROOT=Path(__file__).resolve().parent
FACE={1:[1,0,3,2],2:[0,1,5,4],3:[1,2,6,5],4:[2,3,7,6],5:[3,0,4,7],6:[4,5,6,7]}
CORNER=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
CENTERS={'panel1':[1,50,50],'panel2':[99,50,50],'panel3':[50,1,50],'panel4':[50,99,50],
 'panel5':[50,50,1],'panel6':[50,50,99],'pcb1':[50,50,21],'pcb2':[50,50,38],
 'pcb3':[50,50,62],'pcb4':[50,50,79],'battery':[50,50,43.5]}

def analyze(run):
    m=np.load(run/'mesh.npz');nodes=m['nodes'];conn=m['connectivity'];vol=m['volumes'];tags=m['part_names'];mat=m['materials'];faces=m['faces'];n=len(nodes)
    text=(run/'model.inp').read_text(encoding='utf-8')
    duration=max(float(t) for t in re.findall(r'^TIME,([^\n]+)$',text,re.M));start=duration-5580
    props={}
    for kind,mid,value in re.findall(r'^MP,(DENS|C|KXX),(\d+),([^\n]+)$',text,re.M):props.setdefault(int(mid),{})[kind]=float(value)
    cap=vol*np.array([props[int(mid)]['DENS']*props[int(mid)]['C'] for mid in mat])
    capnod=np.bincount(conn.ravel(),weights=np.repeat(cap/8,8),minlength=n)
    unique=sorted(set(tags));weights=[];partids={}
    for tag in unique:
        mask=tags==tag; partids[tag]=np.unique(conn[mask]);weights.append(np.bincount(conn[mask].ravel(),weights=np.repeat(vol[mask]/8,8),minlength=n)/vol[mask].sum())
    weights=np.array(weights)
    xyz=nodes[conn];low=xyz.min(axis=1);high=xyz.max(axis=1);cw=[]
    for tag,point in CENTERS.items():
        point=np.array(point)/1000;ids=np.flatnonzero((tags==tag)&np.all(point>=low-1e-12,axis=1)&np.all(point<=high+1e-12,axis=1));assert len(ids)
        w=np.zeros(n)
        for e in ids:
            q=np.clip((point-low[e])/(high[e]-low[e]),0,1);np.add.at(w,conn[e],np.prod(np.where(CORNER,q,1-q),axis=1)/len(ids))
        cw.append(w)
    cw=np.array(cw)
    ext=faces[faces[:,2]>0];fc=np.array([conn[int(e)-1][FACE[int(f)]] for e,f,s,a in ext])
    emiss={(int(e),int(f)):float(v) for e,f,v in re.findall(r'^SFE,(\d+),(\d+),RDSF,1,([^\n]+)$',text,re.M)}
    ecoef=np.array([a*emiss[int(e),int(f)]*5.670374419e-8 for e,f,s,a in ext])
    flux=np.loadtxt(run/'applied_flux.csv',delimiter=',',skiprows=1)
    areas=np.array([sum(a for e,f,s,a in ext if s==side) for side in range(1,7)])
    rows=[];centers=[];spreads=[];mins=[];maxs=[];oldphase={};periodic=0.;compared=0;errors=[];through=[];previous=None;snap={}
    extrema={tag:{'minimum_K':float('inf'),'maximum_K':0.,'maximum_simultaneous_spread_K':0.} for tag in unique}
    with (run/'nodal_temperatures.csv').open(encoding='utf-8') as f:
        while True:
            block=[f.readline() for _ in range(n)]
            if not block[0]:break
            assert all(block),'Incomplete nodal export block'
            data=np.fromstring(''.join(block).replace('\n',','),sep=',').reshape(n,3)
            assert np.all(data[:,0]==data[0,0]) and np.array_equal(data[:,1],np.arange(1,n+1))
            t=float(data[0,0]);temp=data[:,2];assert np.isfinite(temp).all() and temp.min()>0
            avg=weights@temp;energy=float(capnod@temp);pin=float(areas@np.array([np.interp(t%5580,flux[:,0],flux[:,20+i]) for i in range(6)]));pout=float(ecoef@((temp[fc].mean(axis=1))**4-2.7**4))
            rows.append([t,*avg,energy,pin,pout]);centers.append([t,*(cw@temp)])
            lo=np.array([temp[partids[tag]].min() for tag in unique]);hi=np.array([temp[partids[tag]].max() for tag in unique]);spread=hi-lo
            mins.append([t,*lo]);maxs.append([t,*hi]);spreads.append([t,*spread])
            if t>=start:
                for j,tag in enumerate(unique):
                    v=extrema[tag];v['minimum_K']=min(v['minimum_K'],float(lo[j]));v['maximum_K']=max(v['maximum_K'],float(hi[j]));v['maximum_simultaneous_spread_K']=max(v['maximum_simultaneous_spread_K'],float(spread[j]))
                if previous is not None:
                    pt,pe=previous;dt=t-pt
                    # Exports may omit warm-up; never treat an unobserved gap as one time step.
                    if 0<dt<=10.000001:errors.append((energy-pe)/dt-(pin-pout));through.append(max(pin,pout))
            if t>=duration-11160-1e-6:
                phase=round(t%5580,6)
                if t>start+1e-6 and phase in oldphase:
                    oldtime,oldtemp=oldphase[phase]
                    if abs(t-oldtime-5580)<1e-5:periodic=max(periodic,float(np.max(abs(temp-oldtemp))));compared+=1
                oldphase[phase]=(t,temp.copy())
            if t>=start and round(t-start,6) in [1720.,3880.]:snap[round(t-start,6)]=temp.copy()
            previous=(t,energy);final=temp.copy()
    d=np.array(rows);assert np.all(np.diff(d[:,0])>0)
    assert abs(d[-1,0]-duration)<1e-6,'Nodal export does not reach the specified final phase-zero time'
    last=d[:,0]>=start-1e-6;phase=d[last,0]-start;avg_last=d[last,1:1+len(unique)]
    maxstep=float(np.diff(d[last,0]).max());expected=int(round(5580/maxstep))
    reference=np.genfromtxt(ROOT/'config/journal_figure11_reference.csv',delimiter=',',names=True)
    comparison={}
    for tag in reference.dtype.names[1:]:
        j=unique.index(tag);diff=np.interp(reference['time_s'],phase,avg_last[:,j])-reference[tag]
        comparison[tag]={'rmse_K':float(np.sqrt(np.mean(diff**2))),'max_abs_difference_K':float(abs(diff).max()),'mean_bias_K':float(diff.mean())}
    for name,header,values in [('part_temperatures.csv',['time_s',*unique,'energy_J','absorbed_power_W','emitted_power_W'],rows),('center_temperatures.csv',['time_s',*CENTERS],centers),('part_spreads.csv',['time_s',*unique],spreads),('part_minimum_temperatures.csv',['time_s',*unique],mins),('part_maximum_temperatures.csv',['time_s',*unique],maxs)]:
        with (run/name).open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(header);w.writerows(values)
    np.savetxt(run/'last_orbit.csv',np.c_[phase,d[last]],delimiter=',',header=','.join(['phase_s','time_s',*unique,'energy_J','absorbed_power_W','emitted_power_W']),comments='')
    np.save(run/'final_temperatures.npy',final)
    if len(snap)==2:np.savez_compressed(run/'phase_fields.npz',nodes=nodes,connectivity=conn,part_names=tags,phase_s=np.array(sorted(snap)),temperatures_K=np.array([snap[p] for p in sorted(snap)]))
    audit={'result_steps':len(d),'final_time_s':duration,'heat_capacity_J_K':float(cap.sum()),'actual_material_properties':props,
     'last_orbit_periodicity_max_node_K':periodic,'phase_comparisons':compared,'expected_phase_comparisons':expected,
     'last_orbit_sampling_max_step_s':maxstep,'periodicity_pass_0p1K':periodic<.1 and compared>=expected,
     'last_orbit_rms_energy_residual_W':float(np.sqrt(np.mean(np.array(errors)**2))),
     'last_orbit_rms_energy_residual_fraction_of_heat_throughput':float(np.sqrt(np.mean(np.array(errors)**2))/np.mean(through)),
     'last_orbit_part_average_range_K':{tag:{'min':float(avg_last[:,j].min()),'max':float(avg_last[:,j].max())} for j,tag in enumerate(unique)},
     'last_orbit_part_nodal_range_K':extrema,'direct_journal_comparison':comparison,
     'pooled_journal_rmse_K':float(np.sqrt(np.mean([v['rmse_K']**2 for v in comparison.values()]))),
     'export_note':'Streaming audit of recorded nodal history. Energy audit excludes unobserved gaps. Periodicity compares matching phases over the final half-open orbit (start,end].',
     'validation_status':'Numerical comparison to journal curves, not experimental validation.'}
    (run/'solution_audit.json').write_bytes((json.dumps(audit,indent=2)+'\n').encode('utf-8'))
    print(json.dumps({k:audit[k] for k in ['last_orbit_periodicity_max_node_K','phase_comparisons','last_orbit_rms_energy_residual_W','pooled_journal_rmse_K']},indent=2))
    return audit

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);a=p.parse_args();analyze(a.case)
