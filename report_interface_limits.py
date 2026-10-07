"""Summarize interface diagnostics without fitting curves or input properties."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
LABELS=['bonded_control','seams_off','pcb_frame_off','both_off']

def report():
    records={}; tables={}
    for label in LABELS:
        folder=ROOT/'data/ansys'/f'interface_{label}'
        settled=ROOT/'data/ansys'/f'interface_{label}_settled'
        if (settled/'solution_audit.json').exists():folder=settled
        qualified=ROOT/'data/ansys'/f'interface_{label}_qualified'
        if (qualified/'solution_audit.json').exists():folder=qualified
        if not folder.exists():folder=ROOT/'results/cases'/f'interface_{label}'
        a=json.loads((folder/'solution_audit.json').read_text(encoding='utf-8'))
        table_path=folder/'part_temperatures.csv'
        if not table_path.exists():table_path=folder/'last_orbit.csv'
        table=np.genfromtxt(table_path,delimiter=',',names=True)
        phase=table['time_s']-(table['time_s'][-1]-5580); mask=phase>=0
        tables[label]=(phase[mask],table[mask])
        comparison=a['direct_journal_comparison']
        records[label]={
            'pooled_journal_rmse_K':float(np.sqrt(np.mean([v['rmse_K']**2 for v in comparison.values()]))),
            'maximum_journal_temperature_difference_K':max(v['max_abs_difference_K'] for v in comparison.values()),
            'periodicity_max_node_K':a['last_orbit_periodicity_max_node_K'],
            'settled_to_0p1K':a['last_orbit_periodicity_max_node_K']<.1,
            'energy_rms_W':a['last_orbit_rms_energy_residual_W'],
            'battery_average_min_K':a['last_orbit_part_average_range_K']['battery']['min'],
            'battery_average_max_K':a['last_orbit_part_average_range_K']['battery']['max'],
            'battery_nodal_range':a.get('last_orbit_part_nodal_range_K',{}).get('battery',{}),
            'component_nodal_ranges':a.get('last_orbit_part_nodal_range_K',{}),
            'per_component_journal_comparison':comparison,
            'interface_provenance':json.loads((folder/'interface_provenance.json').read_text(encoding='utf-8'))}
        initial_path=ROOT/'data/ansys'/f'interface_{label}'/'solution_audit.json'
        if not initial_path.exists():initial_path=folder/'initial_screening_audit.json'
        if initial_path.exists():records[label]['initial_screening_periodicity_K']=json.loads(initial_path.read_text(encoding='utf-8'))['last_orbit_periodicity_max_node_K']
    bp,base=tables['bonded_control']
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for label,(phase,table) in tables.items():
        changes={tag:float(np.max(abs(table[tag]-np.interp(phase,bp,base[tag])))) for tag in records[label]['per_component_journal_comparison']}
        records[label]['maximum_component_average_change_from_control_K']=max(changes.values())
        records[label]['component_average_changes_from_control_K']=changes
        axes[0].plot(phase,table['battery'],label=label.replace('_',' '))
        axes[1].plot(phase,table['panel1'],label=label.replace('_',' '))
    ref=np.genfromtxt(ROOT/'config/journal_figure11_reference.csv',delimiter=',',names=True)
    for ax,tag in zip(axes,['battery','panel1']):
        ax.plot(ref['time_s'],ref[tag],'k--',label='2021 paper')
        ax.set(xlabel='Orbital phase (s)',ylabel='Volume-average temperature (K)',title=tag)
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    output=ROOT/'results' if (ROOT/'results').exists() else ROOT/'results'
    fig.tight_layout();fig.savefig(output/'interface_limits.png',dpi=170);plt.close(fig)
    result={'study':'Direct-interface zero-conductance limits; unchanged geometry and radiation',
      'acceptance_note':'0.1 K nodal periodicity is a study target, not an externally mandated tolerance. Cases above it remain transient screening results.',
      'limitations':'No finite contact-conductance calibration or experimental validation. Coincident sealed interfaces, unchanged radiation surfaces, and alternate paths through frame at junctions remain. This cannot recover author CAD or isolate view-factor algorithm differences. Nodal extremes and spatial spreads are coarse-mesh diagnostics; formal spatial and through-thickness convergence remains necessary.',
      'cases':records}
    (output/'interface_limits.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf-8'))
    lines=['# Direct-interface limiting-case study','',result['study'],
      '', 'Each initial diagnostic starts from the same mapped phase-zero temperature field and solves two complete 10 s orbits with closure/reciprocity and strict heat convergence. Unsettled variants are extended using six 60 s warm-up orbits followed by two 10 s orbits, then two additional 10 s orbits to reduce the remaining drift. Only the final resolved orbit is compared. Material properties, volumes, radiating faces and absorbed loads are fixed. No properties or time shifts are fitted.',
      '', '| Case | Journal RMSE (K) | Maximum curve error (K) | Cycle change (K) | Battery min/max (K) | Max change from control (K) |',
      '|---|---:|---:|---:|---:|---:|']
    for label,r in records.items():
        lines.append(f"| {label} | {r['pooled_journal_rmse_K']:.3f} | {r['maximum_journal_temperature_difference_K']:.3f} | {r['periodicity_max_node_K']:.4f} | {r['battery_average_min_K']:.2f} / {r['battery_average_max_K']:.2f} | {r['maximum_component_average_change_from_control_K']:.3f} |")
    lines+=['', '| Case | Coldest battery node (K) | Hottest battery node (K) | Maximum simultaneous battery spread (K) |', '|---|---:|---:|---:|']
    for label,r in records.items():
        b=r['battery_nodal_range']
        if b:lines.append(f"| {label} | {b['minimum_K']:.3f} | {b['maximum_K']:.3f} | {b['maximum_simultaneous_spread_K']:.3f} |")
    lines+=['', 'Maximum simultaneous nodal temperature spread on each PCB:', '', '| PCB | Bonded control (K) | Seams off (K) | PCB/frame off (K) | Both off (K) |', '|---|---:|---:|---:|---:|']
    for tag in ['pcb1','pcb2','pcb3','pcb4']:
        spreads=[records[label]['component_nodal_ranges'].get(tag,{}).get('maximum_simultaneous_spread_K',float('nan')) for label in LABELS]
        lines.append('| '+tag+' | '+' | '.join(f'{v:.3f}' for v in spreads)+' |')
    pcbs=['pcb1','pcb2','pcb3','pcb4']; x=np.arange(4)
    fig,ax=plt.subplots(figsize=(7.5,4.5))
    for shift,label,color in [(-.18,'bonded_control','#2474ab'),(.18,'pcb_frame_off','#d66a24')]:
        values=[records[label]['component_nodal_ranges'][tag]['maximum_simultaneous_spread_K'] for tag in pcbs]
        bars=ax.bar(x+shift,values,.36,label=label.replace('_',' '),color=color)
        ax.bar_label(bars,fmt='%.1f',padding=3,fontsize=9)
    ax.set(xticks=x,xticklabels=[tag.upper() for tag in pcbs],ylabel='Maximum simultaneous temperature spread (K)',ylim=(0,38),title='PCB spatial temperature sensitivity: coarse-mesh diagnostic')
    ax.legend(loc='upper right',fontsize=9);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    fig.tight_layout();fig.savefig(output/'interface_spreads.png',dpi=170);plt.close(fig)
    spread_change=max(abs(records['bonded_control']['component_nodal_ranges'][tag]['maximum_simultaneous_spread_K']-records['pcb_frame_off']['component_nodal_ranges'][tag]['maximum_simultaneous_spread_K']) for tag in pcbs)
    lines+=['',f"Removing frame-to-PCB contact changes main-component volume-average temperatures by at most {records['pcb_frame_off']['maximum_component_average_change_from_control_K']:.3f} K, while changing a PCB's maximum simultaneous nodal spread by as much as {spread_change:.3f} K. This demonstrates why volume-average agreement alone does not establish local-temperature accuracy. The spatial result requires mesh qualification.", '', '![PCB spatial spreads](interface_spreads.png)']
    best=min(records,key=lambda label:records[label]['pooled_journal_rmse_K'])
    lines+=['',f"The smallest pooled RMSE among these cases is {records[best]['pooled_journal_rmse_K']:.3f} K ({best}), compared with {records['bonded_control']['pooled_journal_rmse_K']:.3f} K for the bonded control."]
    if all(r['settled_to_0p1K'] for r in records.values()):
        lines+=['', 'All four final comparisons meet the prespecified 0.1 K maximum nodal cycle-change target. Global final-orbit RMS power residuals are below '+f"{max(r['energy_rms_W'] for r in records.values()):.3g} W."]
    else:lines+=['', '**Some cases exceed the 0.1 K cycle-change target. Their error comparisons are provisional and must not be presented as settled contact effects.**']
    lines+=['','The pooled RMSE weights the eleven component curves equally. Battery values are volume averages, not hottest battery-node temperatures.', '', result['acceptance_note'], '', result['limitations'], '', 'Panel seam removal deletes only direct panel-to-panel links. Frame-mediated continuity at geometric junctions is retained. PCB/frame removal leaves the bolt paths and battery/PCB contact intact.', '', 'Source rationale: the paper describes the supports as conductive paths to the PCBs and battery attachment to PCB2; it does not fully specify the reconstructed PCB/frame interface locations. Its lumped model explicitly omits direct panel-to-panel conduction, but that does not establish the same assumption for the finite-volume model. These limits therefore test reconstruction ambiguity rather than assert a source correction.', '', '![Interface comparison](interface_limits.png)', '', 'Interpret these runs as hypothesis tests and bounding diagnostics. A smaller paper error does not prove a more accurate physical contact arrangement. Formal mesh/time convergence and source-supported finite contact ranges remain separate studies.']
    (output/'Interface_Limits_Report.md').write_bytes(('\n'.join(lines)+'\n').encode('utf-8'))
    print(json.dumps({label:{k:v for k,v in r.items() if not isinstance(v,dict)} for label,r in records.items()},indent=2))

if __name__=='__main__':report()
