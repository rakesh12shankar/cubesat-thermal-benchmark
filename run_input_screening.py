"""Run declared paired input cases; do not interpret their ranges as measurements."""
from concurrent.futures import ThreadPoolExecutor
import json
import shutil,time
from qualification_study import ROOT,clone_case,solve

def scenario(item):
    name,spec,side=item
    value=spec['nominal']+(-1 if side=='low' else 1)*spec['half_width']
    case=ROOT/'data/ansys'/f'uncertainty_{name}_{side}'
    if (case/'solution_audit.json').exists():
        print(f'Already analyzed: {case.name}',flush=True);return
    if case.exists():
        if (case/'orbital.rth').exists():
            raise RuntimeError(f'Existing unfinished result requires inspection: {case}')
        if (case/'solver.out').exists():shutil.copy2(case/'solver.out',case/f'solver_failed_memory_{time.time_ns()}.out')
    else:
        source=ROOT/'data/ansys/interface_bonded_control'
        if not source.exists():source=ROOT/'data/ansys/captured_convergence_h15_features1'
        case=clone_case(case.name,source,dt=10,cycles=3,changes={name:value})
    print(f'Running {case.name}: {name}={value}',flush=True);solve(case)
    audit=json.loads((case/'solution_audit.json').read_text())
    print(json.dumps({'case':case.name,'periodicity_K':audit['last_orbit_periodicity_max_node_K'],'pass':audit['periodicity_pass_0p1K']}),flush=True)

if __name__=='__main__':
    plan=json.loads((ROOT/'config/input_uncertainty_plan.json').read_text())
    tasks=[(name,spec,side) for name,spec in plan['factors'].items() for side in ['low','high']]
    with ThreadPoolExecutor(max_workers=1) as pool:list(pool.map(scenario,tasks))
