"""Check public files, references, numerical exports and script syntax."""
from pathlib import Path
import ast,json,re,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent

def main():
    issues=[]
    for p in ROOT.rglob('*'):
        if not p.is_file() or '.git' in p.parts or '__pycache__' in p.parts:continue
        if p.name.startswith('~$'):issues.append(f'Office/CAD lock file: {p.relative_to(ROOT)}')
        if p.suffix.lower() in {'.pdf','.rth','.db','.rdb','.esav','.osav','.vf','.full'}:issues.append(f'Unwanted archive file: {p.relative_to(ROOT)}')
        if p.stat().st_size>10_000_000:issues.append(f'Oversized file: {p.relative_to(ROOT)}')
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
        if p.suffix in {'.md','.json','.csv','.py','.ps1','.txt','.cff'}:
            s=p.read_text(encoding='utf-8')
            if re.search(r'C:[/\\]Users[/\\]',s,re.I):issues.append(f'Machine-specific user path: {p.relative_to(ROOT)}')
        if p.suffix=='.md':
            for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
                if link.startswith(('http:','https:','#')):continue
                if not (p.parent/link).exists():issues.append(f'Broken local link in {p.relative_to(ROOT)}: {link}')
    for p in (ROOT/'results/cases').glob('*/last_orbit.csv'):
        d=np.genfromtxt(p,delimiter=',',names=True)
        if not np.all(np.diff(d['phase_s'])>0):issues.append(f'Nonmonotonic phase: {p.parent.name}')
        if d['phase_s'][-1]!=5580 or d['phase_s'][0]>10:issues.append(f'Incomplete final orbit: {p.parent.name}')
        for k in d.dtype.names:
            if not np.isfinite(d[k]).all():issues.append(f'Nonfinite {k}: {p.parent.name}')
    recorded=json.loads((ROOT/'results/journal_validation.json').read_text(encoding='utf-8'))['temperature_comparisons']
    computed=json.loads((ROOT/'results/recomputed_temperature_metrics.json').read_text(encoding='utf-8'))
    for label,case in [('Coarse mesh, orbit 8','orbital_h15_e0.5_qbase'),('Fine mesh, orbit 7','orbital_h10_e0.5_warm')]:
        for tag,v in recorded[label].items():
            if abs(v['rmse']-computed[case][tag]['rmse_K'])>1e-8:issues.append(f'Metric regeneration mismatch: {case}/{tag}')
    verify=json.loads((ROOT/'results/analytical_cooling_verification.json').read_text(encoding='utf-8'))
    if not verify['pass']:issues.append('Analytical verification did not pass.')
    internal=ROOT/'results/closed_enclosure_verification.json'
    if internal.exists() and not json.loads(internal.read_text(encoding='utf-8'))['pass']:issues.append('Closed internal-enclosure check did not pass.')
    manifest=json.loads((ROOT/'SHA256SUMS.json').read_text(encoding='utf-8'))
    # Image encoding/fonts and floating-point formatting can vary by platform.
    # Their numerical source tables and recorded metrics remain hash-checked.
    regenerated={'results/numerical_convergence.png','results/numerical_resolution_differences.png','results/refined_pcb_fields.png','results/input_uncertainty_bands.png','results/input_sensitivity_ranking.png','results/numerical_uncertainty.json','results/interface_spreads.png','results/interface_limits.png','results/temperature_comparison.png','results/emissivity_effect.png','results/load_interpretation.png','results/improved_solver_comparison.png','results/temperature_fields.png','results/recomputed_temperature_metrics.json'}
    for rel,digest in manifest.items():
        if rel in regenerated or (rel.startswith('results/uncertainty_band_') and rel.endswith('.csv')):continue
        if not (ROOT/rel).exists() or hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=digest:issues.append(f'Changed/missing packaged file: {rel}')
    selection=json.loads((ROOT/'config/qualification_case_selection.json').read_text())
    study_cases=list(dict.fromkeys(selection['spatial']+selection['temporal']+selection.get('angular',[])+[v for pair in selection['uncertainty'].values() for v in pair.values()]+[v['case'] for v in selection.get('joint_checks',[])]))
    uncertainty_cases={v for pair in selection['uncertainty'].values() for v in pair.values()}|{v['case'] for v in selection.get('joint_checks',[])}
    target=json.loads((ROOT/'config/input_uncertainty_plan.json').read_text())['periodicity_target_K']
    for name in study_cases:
        case=ROOT/'results/cases'/name
        d=np.genfromtxt(case/'last_orbit.csv',delimiter=',',names=True)
        residual=np.diff(d['energy_J'])/np.diff(d['phase_s'])-(d['absorbed_power_W'][1:]-d['emitted_power_W'][1:])
        if np.sqrt(np.mean(residual**2))>5e-4:issues.append(f'Study energy-balance check failed: {name}')
        audit=json.loads((case/'solution_audit.json').read_text())
        if abs(d['time_s'][-1]-audit['final_time_s'])>1e-6:issues.append(f'Incomplete final orbit: {name}')
        if audit['phase_comparisons']<audit['expected_phase_comparisons']:issues.append(f'Incomplete matching-phase coverage: {name}')
        if not audit['periodicity_pass_0p1K']:issues.append(f'Study periodicity check failed: {name}')
        if name in uncertainty_cases and audit['last_orbit_periodicity_max_node_K']>=target:issues.append(f'Input-scenario periodicity check failed: {name}')
    if issues:raise SystemExit('\n'.join(issues))
    print(f'Release checks passed: {len(manifest)} manifest files; scripts, links, source exclusions and regenerated metrics verified.')

if __name__=='__main__':main()
