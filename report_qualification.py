"""Regenerate spatial/time and assumption-based uncertainty summaries from saved CSVs."""
from pathlib import Path
import argparse,json
from functools import lru_cache
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
ROOT=Path(__file__).resolve().parent
PARTS=['battery','bolts','frame','panel1','panel2','panel3','panel4','panel5','panel6','pcb1','pcb2','pcb3','pcb4']

@lru_cache(maxsize=64)
def load(case):
    audit=json.loads((case/'solution_audit.json').read_text())
    d=np.genfromtxt(case/'last_orbit.csv',delimiter=',',names=True)
    assert np.isfinite(np.column_stack([d[p] for p in PARTS])).all()
    return audit,d

def resolution_label(case,family,mesh=None):
    p=json.loads((case/'study_provenance.json').read_text())
    if family=='temporal':return f'Δt = {p["dt_s"]:g} s'
    if family=='angular':return f'HRES = {p.get("changes",{}).get("hemicube_resolution",20)}'
    mesh=mesh or json.loads((case/'mesh_audit.json').read_text())
    h=mesh.get('maximum_background_spacing_mm',p.get('mesh_h_mm'))
    thin=mesh.get('small_feature_divisions',p.get('thin_region_divisions',1))
    return f'h = {h:g} mm; thin = {thin} ({mesh["solid_elements"]:,} elements)'

def pcb_fields(case,out):
    path=case/'phase_fields.npz'
    if not path.exists():return False
    m=np.load(path);xyz=m['nodes'];conn=m['connectivity'];tags=m['part_names'];fields=m['temperatures_K'];phases=m['phase_s']
    fig,axs=plt.subplots(len(phases),4,figsize=(12,6.5),layout='constrained')
    vals=[];cells=[];topnodes=[]
    for pcb in ['pcb1','pcb2','pcb3','pcb4']:
        elements=conn[tags==pcb];top=xyz[elements,2].max();elements=elements[np.isclose(xyz[elements,2].max(axis=1),top,rtol=0,atol=1e-10)]
        nodes=elements[:,[4,5,6,7]];assert np.allclose(xyz[nodes,2],top)
        cells.append(xyz[nodes,:2]*1000);topnodes.append(nodes)
        vals.extend([temp[nodes].mean(axis=1) for temp in fields])
    vmin=min(v.min() for v in vals);vmax=max(v.max() for v in vals)
    for row,(phase,temp) in enumerate(zip(phases,fields)):
        for col,(polygons,nodes) in enumerate(zip(cells,topnodes)):
            ax=axs[row,col];colors=temp[nodes].mean(axis=1);collection=PolyCollection(polygons,array=colors,cmap='inferno',clim=(vmin,vmax),edgecolors='none');ax.add_collection(collection)
            spread=float(np.ptp(temp[np.unique(nodes)]));ax.set(xlim=(5,95),ylim=(5,95),aspect='equal',xlabel='x (mm)',ylabel='y (mm)',title=f'PCB {col+1}, phase {phase:g} s\nTop-node spread {spread:.1f} K')
            if col==1:ax.plot([20,80,80,20,20],[20,20,80,80,20],'--',color='cyan',lw=.8)
    fig.colorbar(collection,ax=axs.ravel().tolist(),label='Face-average temperature (K)',shrink=.8)
    fig.suptitle(f'{case.name}: bonded PCB top surfaces\nCyan outline: battery footprint')
    fig.savefig(out/'refined_pcb_fields.png',dpi=180);plt.close(fig);return True

def delta(a,b):
    aa,ad=load(a);ba,bd=load(b)
    phase=np.unique(np.r_[ad['phase_s'],bd['phase_s']])
    diff={p:float(np.max(abs(np.interp(phase,bd['phase_s'],bd[p])-np.interp(phase,ad['phase_s'],ad[p])))) for p in PARTS}
    result={'cases':[a.name,b.name],'maximum_component_average_curve_difference_K':max(diff.values()),'by_part_K':diff,
        'battery_nodal_minimum_difference_K':abs(ba['last_orbit_part_nodal_range_K']['battery']['minimum_K']-aa['last_orbit_part_nodal_range_K']['battery']['minimum_K']),
        'battery_nodal_maximum_difference_K':abs(ba['last_orbit_part_nodal_range_K']['battery']['maximum_K']-aa['last_orbit_part_nodal_range_K']['battery']['maximum_K']),
        'pcb_maximum_spread_differences_K':{p:abs(ba['last_orbit_part_nodal_range_K'][p]['maximum_simultaneous_spread_K']-aa['last_orbit_part_nodal_range_K'][p]['maximum_simultaneous_spread_K']) for p in ['pcb1','pcb2','pcb3','pcb4']}}
    result['analyst_screening_targets_pass']=max(diff.values())<=.25 and max(result['battery_nodal_minimum_difference_K'],result['battery_nodal_maximum_difference_K'])<=.25 and max(result['pcb_maximum_spread_differences_K'].values())<=.5
    return result

def main(cases,out,selection):
    out.mkdir(parents=True,exist_ok=True)
    spec=json.loads(selection.read_text());plan=json.loads((ROOT/'config/input_uncertainty_plan.json').read_text())
    convergence={};lines=['# Numerical qualification and input uncertainty','',
        'This study assesses numerical resolution and assumed-input sensitivity of the reconstructed conduction–radiation model. Source-paper curve comparison is numerical benchmarking, not experimental validation. Input intervals are declared engineering assumptions; their probability distributions were not measured.','']
    families=['spatial','temporal']+(['angular'] if 'angular' in spec else [])
    fig,axs=plt.subplots(1,len(families),figsize=(5.5*len(families),4.2))
    detail,dx=plt.subplots(2,len(families),figsize=(5.5*len(families),7))
    for col,family in enumerate(families):
        names=spec[family];records=[];pairs=[]
        for name in names:
            case=cases/name;a,d=load(case);mesh=json.loads((case/'mesh_audit.json').read_text())
            records.append({'case':name,'nodes':mesh['nodes'],'elements':mesh['solid_elements'],'periodicity_K':a['last_orbit_periodicity_max_node_K'],'periodicity_pass':a['periodicity_pass_0p1K'],'energy_residual_W':a['last_orbit_rms_energy_residual_W'],'journal_rmse_K':a['pooled_journal_rmse_K']})
            axs[col].plot(d['phase_s'],d['battery'],label=resolution_label(case,family,mesh))
        for a,b in zip(names,names[1:]):pairs.append(delta(cases/a,cases/b))
        convergence[family]={'cases':records,'adjacent_comparisons':pairs}
        if family=='temporal' and len(names)==3:
            histories=[load(cases/name)[1] for name in names];phase=histories[0]['phase_s'];orders={}
            for part in PARTS:
                values=[np.interp(phase,d['phase_s'],d[part]) for d in histories]
                d01=float(np.sqrt(np.trapezoid((values[1]-values[0])**2,phase)/5580));d12=float(np.sqrt(np.trapezoid((values[2]-values[1])**2,phase)/5580))
                orders[part]={'coarse_medium_rms_difference_K':d01,'medium_fine_rms_difference_K':d12,'observed_order_for_halved_steps':float(np.log2(d01/d12)) if d01>0 and d12>1e-10 else None}
            convergence[family]['observed_order_diagnostics']=orders
        lines += [f'## {family.capitalize()} resolution','', '| Case | Elements | Periodicity (K) | Journal RMSE (K) |','|---|---:|---:|---:|']
        lines += [f'| {v["case"]} | {v["elements"]} | {v["periodicity_K"]:.5f} | {v["journal_rmse_K"]:.4f} |' for v in records]
        lines += ['', '| Adjacent cases | Max average-curve change (K) | Max battery nodal-extreme change (K) | Max PCB spread change (K) |','|---|---:|---:|---:|']
        lines += [f'| {v["cases"][0]} → {v["cases"][1]} | {v["maximum_component_average_curve_difference_K"]:.4f} | {max(v["battery_nodal_minimum_difference_K"],v["battery_nodal_maximum_difference_K"]):.4f} | {max(v["pcb_maximum_spread_differences_K"].values()):.4f} |' for v in pairs]
        lines += ['', f'Final adjacent-pair screening targets: {"PASS" if pairs[-1]["analyst_screening_targets_pass"] else "NOT MET"}. Targets: 0.25 K for average curves and battery nodal extrema; 0.5 K for PCB maximum spreads.','']
        axs[col].set(xlabel='Orbital phase (s)',ylabel='Battery average temperature (K)',title=family.capitalize());axs[col].grid(alpha=.25);axs[col].legend(fontsize=8)
        finest_a,finest=load(cases/names[-1]);bars=[]
        for name in names:
            a,d=load(cases/name)
            if name!=names[-1]:dx[0,col].plot(finest['phase_s'],np.interp(finest['phase_s'],d['phase_s'],d['battery'])-finest['battery'],label=resolution_label(cases/name,family))
            bars.append([a['last_orbit_part_nodal_range_K'][p]['maximum_simultaneous_spread_K'] for p in ['pcb1','pcb2','pcb3','pcb4']])
        dx[0,col].set(title=family.capitalize(),xlabel='Orbital phase (s)',ylabel='Battery average minus finest tested (K)');dx[0,col].grid(alpha=.25);dx[0,col].legend(fontsize=8)
        x=np.arange(4);width=.8/len(names)
        for j,(name,values) in enumerate(zip(names,bars)):dx[1,col].bar(x+(j-(len(names)-1)/2)*width,values,width,label=resolution_label(cases/name,family))
        dx[1,col].set_ylim(0,max(35,float(np.max(bars))*1.28))
        dx[1,col].set(xticks=x,xticklabels=['PCB1','PCB2','PCB3','PCB4'],ylabel='Maximum simultaneous nodal spread (K)');dx[1,col].legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'numerical_convergence.png',dpi=180);plt.close(fig)
    detail.tight_layout();detail.savefig(out/'numerical_resolution_differences.png',dpi=180);plt.close(detail)
    has_fields=pcb_fields(cases/spec['spatial'][-1],out)
    lines += ['The spatial sequence refines background spacing and thin-region divisions together. It is anisotropic; no uniform-ratio Richardson extrapolation or GCI is claimed. Reported differences are resolution diagnostics, not a bound on every local temperature.','']
    lines += ['Solver warning counts are retained in each packaged execution summary. MAPDL warns about radiation temperature offsets and the simultaneous presence of view-factor scaling and space-temperature commands. Temperatures here are absolute kelvin with zero offset. The input audit confirms that closure/reciprocity scaling applies only to the closed internal enclosure 1, while the 2.7 K space temperature applies only to the open exterior enclosure 2; the exterior is not forced closed.','']
    lines += ['Interface-gradient conclusions remain conditional on the original contact-study meshes. Matched refinement of the zero-direct-conductance variants is needed to establish contact-effect independence from numerical resolution.','']
    if 'angular' in spec:lines += ['The angular study changes the hemicube resolution on the same solid mesh and time step. ANSYS documents [HEMIOPT](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/ans_cmd/Hlp_C_HEMIOPT.html) as controlling view-factor accuracy. It tests radiation quadrature separately from solid refinement; closure and energy conservation alone cannot establish view-factor accuracy.','']
    if 'observed_order_diagnostics' in convergence['temporal']:
        lines += ['Time-step RMS differences also provide observed-order diagnostics for the 10→5→2.5 s sequence in the JSON report. These are indicators of approach to an asymptotic regime, not certified error bounds; residual periodicity and phase interpolation can affect small differences.','']
    lines += ['The temporal study holds the 10 s tabulated orbital forcing fixed and lets ANSYS interpolate it at smaller solver steps. It isolates time integration; it does not establish convergence of the forcing time sampling or finite-Earth-disk quadrature. Captured starting states use continuation, and each selected final orbit must independently pass the periodicity check.','']
    finest_a,_=load(cases/spec['spatial'][-1])
    lines += ['## Finest bonded case: source-curve discrepancies','','These comparisons use digitized component-average curves at their recorded orbital phases, without shifting the time axis or fitting properties. Numerical-resolution agreement does not remove geometry, forcing or source-model differences.','','| Part | RMSE (K) | Maximum absolute difference (K) | Mean bias (K) |','|---|---:|---:|---:|']
    lines += [f'| {p} | {v["rmse_K"]:.3f} | {v["max_abs_difference_K"]:.3f} | {v["mean_bias_K"]:.3f} |' for p,v in finest_a['direct_journal_comparison'].items()]
    lines += ['','## Finest bonded case: local temperatures','','Extrema cover the final recorded orbit. Maximum nodal spread is the largest simultaneous within-part difference, rather than the difference between extrema at unrelated phases.','','| Part | Average minimum (K) | Average maximum (K) | Nodal minimum (K) | Nodal maximum (K) | Maximum nodal spread (K) |','|---|---:|---:|---:|---:|---:|']
    for p in ['battery','pcb1','pcb2','pcb3','pcb4']:
        av=finest_a['last_orbit_part_average_range_K'][p];nv=finest_a['last_orbit_part_nodal_range_K'][p]
        lines += [f'| {p} | {av["min"]:.3f} | {av["max"]:.3f} | {nv["minimum_K"]:.3f} | {nv["maximum_K"]:.3f} | {nv["maximum_simultaneous_spread_K"]:.3f} |']
    lines += ['']
    base_a,base=load(cases/spec['uncertainty_control']);phase=base['phase_s'];factors=list(plan['factors']);sens=[];curv=[];endpoint=[]
    for factor in factors:
        lo_a,lo=load(cases/spec['uncertainty'][factor]['low']);hi_a,hi=load(cases/spec['uncertainty'][factor]['high'])
        low=np.column_stack([np.interp(phase,lo['phase_s'],lo[p]) for p in PARTS]);high=np.column_stack([np.interp(phase,hi['phase_s'],hi[p]) for p in PARTS]);nom=np.column_stack([base[p] for p in PARTS])
        s=(high-low)/2;c=(high+low)/2-nom;sens.append(s);curv.append(c)
        metrics={}
        for p in ['battery','pcb1','pcb2','pcb3','pcb4']:
            metrics[p]={}
            for metric in ['minimum_K','maximum_K','maximum_simultaneous_spread_K']:
                lv=lo_a['last_orbit_part_nodal_range_K'][p][metric];hv=hi_a['last_orbit_part_nodal_range_K'][p][metric];nv=base_a['last_orbit_part_nodal_range_K'][p][metric]
                metrics[p][metric]={'nominal':nv,'low':lv,'high':hv,'half_range':(hv-lv)/2,'midpoint_curvature':(hv+lv)/2-nv}
        endpoint.append({'factor':factor,'low_case':spec['uncertainty'][factor]['low'],'high_case':spec['uncertainty'][factor]['high'],'maximum_half_range_average_curve_K':float(abs(s).max()),'maximum_midpoint_curvature_K':float(abs(c).max()),'maximum_battery_half_range_K':float(abs(s[:,0]).max()),'nodal_metrics':metrics,'low_periodicity_K':lo_a['last_orbit_periodicity_max_node_K'],'high_periodicity_K':hi_a['last_orbit_periodicity_max_node_K'],'periodicity_pass':lo_a['periodicity_pass_0p1K'] and hi_a['periodicity_pass_0p1K']})
    sens=np.array(sens);draw=np.random.default_rng(plan['reproducibility_seed']).uniform(-1,1,(plan['propagation_draws'],len(factors)))
    reference=np.genfromtxt(ROOT/'config/journal_figure11_reference.csv',delimiter=',',names=True)
    bands=[];trajectory_extrema={};variance_shares={};coverage={};fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    for j,p in enumerate(PARTS):
        # One component at a time limits memory; each draw is shared across all phases.
        values=draw@sens[:,:,j]+base[p][None,:];q=np.quantile(values,[.05,.5,.95],axis=0)
        trajectory_extrema[p]={'minimum_K_p05_p50_p95':np.quantile(values.min(axis=1),[.05,.5,.95]).tolist(),'maximum_K_p05_p50_p95':np.quantile(values.max(axis=1),[.05,.5,.95]).tolist()}
        v=np.trapezoid(sens[:,:,j]**2,phase,axis=1)/(3*5580)
        variance_shares[p]={factor:float(value/v.sum()) if v.sum()>0 else 0 for factor,value in zip(factors,v)}
        if p in reference.dtype.names:
            lower=np.interp(reference['time_s'],phase,q[0]);upper=np.interp(reference['time_s'],phase,q[2]);actual=reference[p]
            coverage[p]={'fraction_of_digitized_reference_samples_in_band':float(np.mean((actual>=lower)&(actual<=upper))),'maximum_departure_outside_band_K':float(max(np.max(lower-actual),np.max(actual-upper),0))}
        bands.append(np.c_[phase,q.T]);
        if p in ['battery','pcb2']:
            ax=axs[0 if p=='battery' else 1];ax.fill_between(phase,q[0],q[2],alpha=.3,label='5–95% affine scenario band');ax.plot(phase,base[p],label='Nominal');ax.plot(reference['time_s'],reference[p],'--',color='black',lw=1.1,label='2021 digitized curve');ax.set(title=p,xlabel='Orbital phase (s)',ylabel='Average temperature (K)');ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.suptitle('Assumed independent uniform ranges: original mesh, 10 s steps, HRES 20',fontsize=11)
    fig.tight_layout();fig.savefig(out/'input_uncertainty_bands.png',dpi=180);plt.close(fig)
    for p,d in zip(PARTS,bands):np.savetxt(out/f'uncertainty_band_{p}.csv',d,delimiter=',',header='phase_s,p05_K,p50_K,p95_K',comments='')
    recorded=out/'uncertainty_bands_recorded.npz'
    if recorded.exists():
        anchor=np.load(recorded)
        for p,d in zip(PARTS,bands):assert np.allclose(d,anchor[p],rtol=0,atol=1e-8),f'Uncertainty-band regeneration mismatch: {p}'
    uncertainty={'assumptions':plan,'control_case':spec['uncertainty_control'],'endpoint_checks':endpoint,'phase_integrated_affine_variance_shares':variance_shares,'affine_trajectory_extrema_quantiles':trajectory_extrema,'band_note':'Conditional pointwise local affine propagation with independent assumed uniform ranges. These are scenario quantiles, not measured confidence intervals or full nonlinear joint ANSYS results. The same sampled inputs apply at every phase of each trajectory.','numerical_resolution_note':spec.get('numerical_resolution_note','See spatial and temporal comparisons separately.')}
    joint=[]
    for item in spec.get('joint_checks',[]):
        a,d=load(cases/item['case']);z=np.array([item['normalized_inputs'][f] for f in factors]);pred=np.column_stack([base[p] for p in PARTS])+np.einsum('f,ftj->tj',z,sens)
        error={p:float(abs(np.interp(phase,d['phase_s'],d[p])-pred[:,j]).max()) for j,p in enumerate(PARTS)}
        joint.append({**item,'maximum_average_curve_affine_error_K':max(error.values()),'by_part_K':error,'periodicity_K':a['last_orbit_periodicity_max_node_K'],'periodicity_pass':a['periodicity_pass_0p1K']})
    uncertainty['nonlinear_joint_checks']=joint
    uncertainty['conditional_source_band_comparison']=coverage
    uncertainty['source_band_comparison_note']='Coverage is a descriptive fraction of digitized source-curve samples under analyst-chosen bands, not a probabilistic model-validation test. Geometry, forcing interpretation and model-form discrepancies remain outside the input ranges.'
    resolution={}
    for part,band in zip(PARTS,bands):
        observed={}
        for family in families:
            values=[]
            for name in spec[family]:
                _,history=load(cases/name);times=np.unique(np.r_[base['phase_s'],history['phase_s']])
                values.append(float(abs(np.interp(times,history['phase_s'],history[part])-np.interp(times,base['phase_s'],base[part])).max()))
            observed[family]=max(values)
        resolution[part]={'largest_observed_nominal_curve_shifts_K':observed,'sum_of_largest_observed_shifts_K':sum(observed.values()),'largest_conditional_band_half_width_K':float(np.max((band[:,3]-band[:,1])/2)),'interpretation':'Separate nominal-resolution diagnostics; their sum is not a certified error bound or a probability interval. Resolution/input interactions were not exhaustively tested.'}
    uncertainty['separate_nominal_resolution_diagnostics']=resolution
    labels={'alpha':'Solar/albedo α','external_epsilon':'External ε','internal_epsilon':'Internal ε','capacity_scale':'Specific heat','panel_k':'Panel k','pcb_k':'PCB k'}
    ordered=sorted(endpoint,key=lambda v:v['maximum_battery_half_range_K']);y=np.arange(len(ordered));fig,axes=plt.subplots(1,2,figsize=(10.5,4.2),layout='constrained')
    axes[0].barh(y,[v['maximum_battery_half_range_K'] for v in ordered],color='#397ca8')
    axes[1].barh(y,[max(abs(v['nodal_metrics'][p]['maximum_simultaneous_spread_K']['half_range']) for p in ['pcb1','pcb2','pcb3','pcb4']) for v in ordered],color='#c87536')
    for ax in axes:ax.set(yticks=y,yticklabels=[labels[v['factor']] for v in ordered]);ax.grid(axis='x',alpha=.25);ax.set_axisbelow(True)
    axes[0].set(xlabel='Battery average-curve half-range (K)',title='Component average sensitivity');axes[1].set(xlabel='Maximum PCB spread half-range (K)',title='Local spread sensitivity')
    fig.suptitle('Assumed input ranges: original mesh, 10 s time step, hemicube resolution 20');fig.savefig(out/'input_sensitivity_ranking.png',dpi=180);plt.close(fig)
    lines += ['## Input uncertainty','',uncertainty['band_note'],'',f'{plan["propagation_draws"]:,} reproducible surrogate evaluations use seed {plan["reproducibility_seed"]}; these are not thousands of ANSYS solves. The full solver supplies the twelve paired endpoints and the joint checks.','',uncertainty['numerical_resolution_note'],'', '| Input | Largest average-curve half-range (K) | Battery half-range (K) | Largest midpoint curvature (K) |','|---|---:|---:|---:|']
    lines += [f'| {v["factor"]} | {v["maximum_half_range_average_curve_K"]:.3f} | {v["maximum_battery_half_range_K"]:.3f} | {v["maximum_midpoint_curvature_K"]:.3f} |' for v in sorted(endpoint,key=lambda x:-x['maximum_half_range_average_curve_K'])]
    lines += ['', '| Input | Phase-integrated battery affine variance share |','|---|---:|']
    lines += [f'| {f} | {100*v:.2f}% |' for f,v in sorted(variance_shares['battery'].items(),key=lambda item:-item[1])]
    lines += ['', '| Part | Sum of largest observed nominal resolution shifts (K) | Largest conditional band half-width (K) |','|---|---:|---:|']
    lines += [f'| {p} | {v["sum_of_largest_observed_shifts_K"]:.3f} | {v["largest_conditional_band_half_width_K"]:.3f} |' for p,v in resolution.items()]
    lines += ['', 'The resolution column adds the largest nominal changes from the separate mesh, time-step and angular sequences. It is a diagnostic comparison of magnitudes, not a certified combined error bound. It is kept separate from the assumed input probabilities.','']
    lines += ['', 'Small sensitivities and curvature values can be comparable with residual periodicity or numerical-resolution effects. Their precise ranking is not established by the input bands. Each endpoint audit retains its actual periodicity defect.','']
    if joint:
        lines += ['',plan.get('joint_initial_state_note',''),'']
        lines += ['', '| Joint case | Largest affine prediction error (K) | Battery error (K) |','|---|---:|---:|']
        lines += [f'| {v["case"]} | {v["maximum_average_curve_affine_error_K"]:.3f} | {v["by_part_K"]["battery"]:.3f} |' for v in joint]
        lines += ['','These two joint checks quantify approximation error at specific combined-input corners. They do not establish an error bound throughout the six-dimensional input domain or certify the affine quantile estimates.','']
    lines += ['', '## Conditional comparison with source curves','',uncertainty['source_band_comparison_note'],'','| Part | Digitized samples inside assumed 5–95% band | Largest out-of-band departure (K) |','|---|---:|---:|']
    lines += [f'| {p} | {100*v["fraction_of_digitized_reference_samples_in_band"]:.1f}% | {v["maximum_departure_outside_band_K"]:.3f} |' for p,v in coverage.items()]
    lines += ['','Curvature measures the endpoint midpoint relative to the nominal response. Significant curvature limits the affine approximation. The study excludes uncertain geometry, finite contact conductance, attitude/ephemeris uncertainty, temperature-dependent properties, and source-solver model-form differences. No input was tuned to the journal curves.','', '![Numerical convergence](numerical_convergence.png)','', '![Resolution differences and local spreads](numerical_resolution_differences.png)','', '![Conditional input bands](input_uncertainty_bands.png)','', '![Average and local sensitivities](input_sensitivity_ranking.png)','']
    if has_fields:lines += ['![Refined PCB temperatures](refined_pcb_fields.png)','','The maps show face-average colors on PCB top surfaces at two recorded phases. Reported spreads use the surface nodes at that phase; the numerical tables separately retain whole-component extrema and spreads over the entire orbit. White areas lie outside the PCB material, including frame and bolt cutouts.','']
    result={'selection':spec,'convergence':convergence,'input_uncertainty':uncertainty}
    (out/'numerical_uncertainty.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
    (out/'Numerical_Uncertainty_Report.md').write_bytes(('\n'.join(lines)).encode())
    anchor=out/'numerical_uncertainty_recorded.json'
    if anchor.exists():
        def same(a,b):
            if isinstance(a,dict):return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
            if isinstance(a,list):return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
            if isinstance(a,(int,float)) and not isinstance(a,bool):return abs(a-b)<1e-8
            return a==b
        assert same(result,json.loads(anchor.read_text())),'Numerical-summary regeneration mismatch'
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cases',type=Path,default=ROOT/'results/cases');p.add_argument('--out',type=Path,default=ROOT/'results');p.add_argument('--selection',type=Path,required=True);a=p.parse_args();main(a.cases,a.out,a.selection)
