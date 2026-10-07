"""Verify conservation and equilibration in an isolated gray radiation enclosure."""
from pathlib import Path
import argparse,json,re,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parent

def main(ansys):
    run=ROOT/'data/ansys/verification_closed_enclosure_strict'
    if run.exists():raise FileExistsError(run)
    run.mkdir()
    # Six independent 100x100x2 mm walls, with their inward faces enclosing
    # a 100 mm cube. Duplicate corner nodes prevent inter-wall conduction.
    walls=[([-.002,0,0],[0,.1,.1],3),([.1,0,0],[.102,.1,.1],5),
           ([0,-.002,0],[.1,0,.1],4),([0,.1,0],[.1,.102,.1],2),
           ([0,0,-.002],[.1,.1,0],6),([0,0,.1],[.1,.1,.102],1)]
    corners=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    lines=['/BATCH','/CLEAR,START','/PREP7','ET,1,SOLID278','TOFFST,0','STEF,5.670374419E-8','MP,DENS,1,1000','MP,C,1,1000','MP,KXX,1,1000']
    for wall,(lo,hi,face) in enumerate(walls):
        xyz=np.array(lo)+(np.array(hi)-lo)*corners;ids=np.arange(wall*8+1,wall*8+9)
        for nid,p in zip(ids,xyz):lines.append(f'N,{nid},'+','.join(f'{x:.12g}' for x in p))
        lines.extend([f'EN,{wall+1},'+','.join(map(str,ids)),f'SFE,{wall+1},{face},RDSF,1,.5',f'SFE,{wall+1},{face},RDSF,2,1'])
    lines+=['FINISH','/AUX12','VFSM,DEFINE,1,2,1000,1E-8','FINISH','/SOLU','ANTYPE,TRANS','TRNOPT,FULL','THOPT,FULL','TUNIF,300']
    for nid in range(1,9):lines.append(f'IC,{nid},TEMP,400')
    lines+=['CNVTOL,HEAT,,1E-8,2,1E-6','HEMIOPT,40,1E-6','RADOPT,,1E-5,2,1000,1E-6,1','AUTOTS,OFF','DELTIM,1','OUTRES,ALL,ALL','TIME,500','SOLVE','FINISH','/POST1','SET,LAST','*GET,NSETS,ACTIVE,0,SET,NSET','*DIM,TEMPALL,ARRAY,48','*DIM,IDS,ARRAY,48','*VFILL,IDS(1),RAMP,1,1','*CFOPEN,temperatures,csv','*DO,I,1,NSETS','SET,,,,,,,I','*GET,TT,ACTIVE,0,SET,TIME','*VGET,TEMPALL(1),NODE,1,TEMP','*VWRITE,TT,IDS(1),TEMPALL(1)',"(F12.6,',',F5.0,',',E20.12)",'*ENDDO','*CFCLOS','FINISH','/EXIT,NOSAVE']
    (run/'model.inp').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    subprocess.run([ansys,'-b','-smp','-np','2','-j','enclosure','-i','model.inp','-o','solver.out'],cwd=run,check=True)
    counts=re.findall(r'NUMBER OF ERROR\s+MESSAGES ENCOUNTERED\s*=\s*(\d+)',(run/'solver.out').read_text(errors='replace'));assert counts and int(counts[-1])==0
    rows=np.loadtxt(run/'temperatures.csv',delimiter=',').reshape(-1,48,3)
    temps=rows[:,:,2].reshape(-1,6,8).mean(axis=2);means=temps.mean(axis=1)
    initial=np.array([400,300,300,300,300,300]);exact_mean=initial.mean()
    deviation=float(abs(means-exact_mean).max());variance_ratio=float(np.var(temps[-1])/np.var(initial))
    result={'description':'Six separate equal-capacity walls, each 20 J/K; epsilon .5; closed 100 mm cubic enclosure; initial wall temperatures [400,300,300,300,300,300] K. Outer faces adiabatic; no conduction between walls. Exact total stored energy and mean temperature are constant. This verifies gray internal radiation conservation/equilibration, not detailed CubeSat view factors.','exact_capacity_weighted_mean_K':float(exact_mean),'maximum_mean_temperature_drift_K':deviation,'maximum_total_energy_drift_J':deviation*120,'final_wall_average_temperatures_K':temps[-1].tolist(),'final_to_initial_wall_temperature_variance_ratio':variance_ratio,'solver_errors':0,'pass':bool(deviation<.005 and variance_ratio<1 and temps.min()>=299.99 and temps.max()<=400.01)}
    (ROOT/'results/closed_enclosure_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if not result['pass']:raise RuntimeError('Closed-enclosure conservation check failed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--ansys',required=True);a=p.parse_args();main(a.ansys)
