"""Audit exported MAPDL temperatures with the actual mesh and surface loads."""
from pathlib import Path
import argparse,json,csv,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT=Path(__file__).resolve().parent
PROPS={1:(2325,1103),2:(2810,948),3:(2120,975),4:(2247,1110)}
FACE_NODES={1:[1,0,3,2],2:[0,1,5,4],3:[1,2,6,5],4:[2,3,7,6],5:[3,0,4,7],6:[4,5,6,7]}

def analyze(run):
    mesh=np.load(run/'mesh.npz');nodes=mesh['nodes'];conn=mesh['connectivity'];vol=mesh['volumes'];mat=mesh['materials'];tags=mesh['part_names'];faces=mesh['faces']
    rows=np.loadtxt(run/'nodal_temperatures.csv',delimiter=',')
    assert len(rows)%len(nodes)==0
    rows=rows.reshape(-1,len(nodes),3);times=rows[:,0,0];temp=rows[:,:,2]
    assert np.all(rows[:,:,1]==np.arange(1,len(nodes)+1))
    assert np.isfinite(temp).all() and temp.min()>0
    etemp=temp[:,conn].mean(axis=2)
    capacity=vol*np.array([PROPS[m][0]*PROPS[m][1] for m in mat])
    energy=etemp@capacity
    unique=sorted(set(tags));averages={tag:np.average(etemp[:,tags==tag],axis=1,weights=vol[tags==tag]) for tag in unique}
    centers={'panel1':[1,50,50],'panel2':[99,50,50],'panel3':[50,1,50],'panel4':[50,99,50],
             'panel5':[50,50,1],'panel6':[50,50,99],'pcb1':[50,50,21],'pcb2':[50,50,38],
             'pcb3':[50,50,62],'pcb4':[50,50,79],'battery':[50,50,43.5]}
    center_series={}
    xyz=nodes[conn];low=xyz.min(axis=1);high=xyz.max(axis=1)
    corner=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    for tag,point in centers.items():
        point=np.array(point)/1000
        candidates=np.flatnonzero((tags==tag)&np.all(point>=low-1e-12,axis=1)&np.all(point<=high+1e-12,axis=1))
        assert len(candidates),tag
        sampled=[]
        for e in candidates:
            q=np.clip((point-low[e])/(high[e]-low[e]),0,1)
            weight=np.prod(np.where(corner,q,1-q),axis=1)
            sampled.append(temp[:,conn[e]]@weight)
        center_series[tag]=np.mean(sampled,axis=0)
    with (run/'center_temperatures.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s',*centers]);w.writerows([t,*[center_series[k][j] for k in centers]] for j,t in enumerate(times))
    power_out=np.zeros(len(times))
    for eid,face,side,area in faces:
        if side>0:
            ft=temp[:,conn[int(eid)-1][FACE_NODES[int(face)]]].mean(axis=1)
            power_out+=area*.72*5.670374419e-8*(ft**4-2.7**4)
    model_text=(run/'model.inp').read_text(encoding='utf-8') if (run/'model.inp').exists() else ''
    mode='orbital' if '*DIM,Q1,TABLE' in model_text else 'cooldown'
    surface_emissivities=[float(v) for v in re.findall(r'SFE,\d+,\d+,RDSF,1,([^\s]+)',model_text)]
    epsilon_half=any(abs(v-.5)<1e-9 for v in surface_emissivities)
    if mode=='orbital':
        flux_path=run/'applied_flux.csv'
        if not flux_path.exists():flux_path=ROOT/'config/orbital_flux_reconstruction.csv'
        flux=np.loadtxt(flux_path,delimiter=',',skiprows=1)
        # exact six external areas are 0.01 m2 each
        power_in=sum(.01*np.interp(times%5580,flux[:,0],flux[:,20+i]) for i in range(6))
    else:power_in=np.zeros(len(times))
    # Backward-Euler time integration should close this discrete global balance.
    denergy=np.diff(energy);step=np.diff(times);rhs=(power_in[1:]-power_out[1:])*step
    residual=denergy-rhs
    denom=np.maximum(np.maximum(np.abs(denergy),np.abs(rhs)),1e-8)
    audit={'result_steps':len(times),'final_time_s':float(times[-1]),'temperature_min_K':float(temp.min()),'temperature_max_K':float(temp.max()),
           'heat_capacity_J_K':float(capacity.sum()),'maximum_discrete_energy_residual_J':float(np.max(abs(residual))) if len(residual) else None,
           'maximum_relative_energy_residual':float(np.max(abs(residual)/denom)) if len(residual) else None,
           'rms_energy_residual_W':float(np.sqrt(np.mean((residual/step)**2))) if len(residual) else None,
           'final_part_average_K':{k:float(v[-1]) for k,v in averages.items()},
           'validation_status':'Solver and conservation audit of a reconstructed model; paper agreement is evaluated separately.'}
    audit['center_point_coordinates_mm']=centers
    if (run/'initial_temperatures.npy').exists():
        initial=np.load(run/'initial_temperatures.npy')
        key='seeded_cycle_end_minus_start_max_node_K' if times[-1]<=5580 else 'seeded_run_end_minus_start_max_node_K'
        audit[key]=float(np.max(abs(temp[-1]-initial)))
        audit['initial_state_note']='Imported stored nodal temperatures; equivalent phase-zero origin is documented in the case provenance.'
    if times[-1]>=2*5580:
        mask=times>=times[-1]-5580
        audit['last_orbit_part_average_range_K']={k:{'min':float(v[mask].min()),'max':float(v[mask].max())} for k,v in averages.items()}
        audit['last_orbit_part_nodal_range_K']={}
        for tag in unique:
            ids=np.unique(conn[tags==tag])
            local_temperature=temp[mask][:,ids]
            audit['last_orbit_part_nodal_range_K'][tag]={
                'minimum_K':float(local_temperature.min()),
                'maximum_K':float(local_temperature.max()),
                'maximum_simultaneous_spread_K':float(np.ptp(local_temperature,axis=1).max())}
        final_steps=times[1:]>=times[-1]-5580
        rms=float(np.sqrt(np.mean((residual[final_steps]/step[final_steps])**2)))
        scale=float(np.mean(np.maximum(power_in[1:][final_steps],power_out[1:][final_steps])))
        audit['last_orbit_rms_energy_residual_W']=rms
        audit['last_orbit_rms_energy_residual_fraction_of_heat_throughput']=rms/scale
        audit['energy_residual_note']='Relative step error is ill-conditioned when net stored-energy change approaches zero; use RMS power residual and heat-throughput normalization.'
        audit['last_orbit_periodicity_max_node_K']=float(np.max(abs(temp[mask]-np.stack([np.interp(times[mask]-5580,times,temp[:,i]) for i in range(len(nodes))],axis=1))))
        audit['periodicity_pass_1K']=audit['last_orbit_periodicity_max_node_K']<1
    if mode=='orbital' and epsilon_half and times[-1]>=5580:
        ref=np.genfromtxt(ROOT/'config/thesis_figure48_vector_reference.csv',delimiter=',',names=True)
        local=times-(times[-1]-5580);mask=local>=0
        comparison={}
        fig,ax=plt.subplots(1,2,figsize=(12,4.5))
        for tag in ref.dtype.names[1:]:
            sampled=np.interp(local[mask],ref['time_s'],ref[tag])
            diff=averages[tag][mask]-sampled
            comparison[tag]={'rmse_K':float(np.sqrt(np.mean(diff**2))),
                             'max_abs_difference_K':float(np.max(abs(diff))),
                             'mean_bias_K':float(np.mean(diff))}
            a=ax[0 if tag.startswith('panel') else 1]
            line=a.plot(local[mask],averages[tag][mask],label=tag)[0]
            a.plot(ref['time_s'],ref[tag],color=line.get_color(),ls='--',alpha=.65)
        for a in ax:a.set(xlabel='Last-orbit time (s)',ylabel='Volume-average temperature (K)');a.grid(alpha=.2);a.legend(fontsize=8)
        fig.suptitle('Solid: ANSYS reconstruction; dashed: supporting thesis Figure 48 (orbit definitions differ)')
        fig.tight_layout();fig.savefig(run/'thesis_comparison.png',dpi=170);plt.close(fig)
        audit['supporting_thesis_comparison']=comparison
        audit['comparison_limitation']='Supporting thesis, not journal validation; 5580 s reconstruction versus approximately 5600 s plotted reference, geometric/contact assumptions, and view-factor method differences.'
        journal_path=ROOT/'config/journal_figure11_reference.csv'
        if journal_path.exists():
            journal=np.genfromtxt(journal_path,delimiter=',',names=True)
            direct={}
            for tag in journal.dtype.names[1:]:
                diff=np.interp(journal['time_s'],local,averages[tag])-journal[tag]
                direct[tag]={'rmse_K':float(np.sqrt(np.mean(diff**2))),
                             'max_abs_difference_K':float(np.max(abs(diff))),
                             'mean_bias_K':float(diff.mean())}
            audit['direct_journal_comparison']=direct
            audit['journal_comparison_note']='Figure 11 volume averages, epsilon 0.5; original time axis, no fitting. Source and geometry ambiguities remain. This is numerical comparison, not experimental validation.'
    (run/'solution_audit.json').write_text(json.dumps(audit,indent=2))
    with (run/'part_temperatures.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s',*unique,'energy_J','absorbed_power_W','emitted_power_W'])
        for j,t in enumerate(times):w.writerow([t,*[averages[k][j] for k in unique],energy[j],power_in[j],power_out[j]])
    fig,ax=plt.subplots(1,2,figsize=(12,4.5))
    for tag,series in averages.items():
        if tag.startswith('panel'):ax[0].plot(times,series,label=tag)
        elif tag.startswith('pcb') or tag=='battery':ax[1].plot(times,series,label=tag)
    for a in ax:a.set(xlabel='Time (s)',ylabel='Volume-average temperature (K)');a.grid(alpha=.2);a.legend(fontsize=8)
    fig.suptitle('ANSYS solid conduction and radiation - '+run.name+' (reconstructed geometry/loads)')
    fig.tight_layout();fig.savefig(run/'temperature_history.png',dpi=170);plt.close(fig)
    # Internal solids only, with cell faces consolidated by geometry mesh.
    fig=plt.figure(figsize=(7,6));ax=fig.add_subplot(111,projection='3d')
    polygons=[];colors=[]
    color={'frame':'#536779','bolts':'#657c9c','battery':'#b84a7c','pcb1':'#36b6b0','pcb2':'#36b6b0','pcb3':'#36b6b0','pcb4':'#36b6b0'}
    for eid,face,side,area in faces:
        e=int(eid)-1;tag=tags[e]
        if tag.startswith('panel'):continue
        polygons.append(nodes[conn[e][FACE_NODES[int(face)]]]*1000);colors.append(color[tag])
    ax.add_collection3d(Poly3DCollection(polygons,facecolors=colors,edgecolors='none',alpha=.95))
    ax.set(xlim=(0,100),ylim=(0,100),zlim=(0,100),xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)')
    ax.set_box_aspect([1,1,1]);ax.view_init(elev=24,azim=-55);ax.set_title('CubeSat 2021 reconstruction - exterior panels hidden')
    fig.tight_layout();fig.savefig(ROOT/'data/cad/geometry_preview.png',dpi=180);plt.close(fig)
    print(json.dumps(audit,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');a=p.parse_args();analyze(Path(a.run))
