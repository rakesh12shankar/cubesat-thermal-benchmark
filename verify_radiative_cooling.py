"""Check MAPDL against the analytical cooling of a uniform radiating cube."""
from pathlib import Path
import argparse, json, re, subprocess
import numpy as np
ROOT=Path(__file__).resolve().parent
SIGMA=5.670374419e-8

def main(ansys):
    summary={}
    for dt in [10,1]:
        run=ROOT/'data/ansys'/f'verification_cooling_dt{dt}'
        if run.exists():raise FileExistsError(f'Existing verification protected: {run}')
        run.mkdir()
        lines=['/BATCH','/CLEAR,START','/PREP7','ET,1,SOLID278','TOFFST,0',f'STEF,{SIGMA}',
               'MP,DENS,1,1000','MP,C,1,1000','MP,KXX,1,1000000']
        xyz=[[0,0,0],[.1,0,0],[.1,.1,0],[0,.1,0],[0,0,.1],[.1,0,.1],[.1,.1,.1],[0,.1,.1]]
        for i,p in enumerate(xyz,1):lines.append('N,'+str(i)+','+','.join(map(str,p)))
        lines.append('EN,1,1,2,3,4,5,6,7,8')
        for face in range(1,7):lines.extend([f'SFE,1,{face},RDSF,1,.72',f'SFE,1,{face},RDSF,2,1'])
        lines+=['FINISH','/SOLU','ANTYPE,TRANS','TRNOPT,FULL','THOPT,FULL','TUNIF,300','SPCTEMP,1,0',
                'HEMIOPT,20,1E-6','RADOPT,,1E-5,2,1000,1E-6,1','AUTOTS,OFF',f'DELTIM,{dt}',
                'OUTRES,ALL,ALL','TIME,1000','SOLVE','FINISH','/POST1','SET,LAST',
                '*GET,NSETS,ACTIVE,0,SET,NSET','*CFOPEN,cooling,csv','*DO,I,1,NSETS','SET,,,,,,,I',
                '*GET,TT,ACTIVE,0,SET,TIME','*GET,TEMP1,NODE,1,TEMP','*VWRITE,TT,TEMP1',
                "(F12.6,',',E20.12)",'*ENDDO','*CFCLOS','FINISH','/EXIT,NOSAVE']
        (run/'model.inp').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        subprocess.run([ansys,'-b','-smp','-np','2','-j','cooling','-i','model.inp','-o','solver.out'],cwd=run,check=True)
        log=(run/'solver.out').read_text(errors='replace');counts=re.findall(r'NUMBER OF ERROR\s+MESSAGES ENCOUNTERED\s*=\s*(\d+)',log)
        assert counts and int(counts[-1])==0,run
        data=np.loadtxt(run/'cooling.csv',delimiter=',')
        exact=(300**-3+3*.72*SIGMA*.06*data[:,0]/1000)**(-1/3)
        error=data[:,1]-exact
        summary[str(dt)]={'maximum_step_s':dt,'maximum_temperature_error_K':float(abs(error).max()),'final_MAPDL_K':float(data[-1,1]),'final_analytical_K':float(exact[-1]),'solver_errors':0}
    summary['description']='Uniform 0.1 m cube; rho=1000 kg/m3, c=1000 J/(kg K), k=1e6 W/(m K), epsilon=.72, all six faces radiate to a 0 K sink; no applied heat. Equal nodal temperatures follow from symmetric loads. T(t)=[T0^-3+3 epsilon sigma A t/(rho c V)]^-1/3. This verifies the external radiation implementation and transient integration, not the CubeSat geometry or internal view factors.'
    summary['pass']=summary['1']['maximum_temperature_error_K']<.02 and summary['1']['maximum_temperature_error_K']<summary['10']['maximum_temperature_error_K']
    (ROOT/'results/analytical_cooling_verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    if not summary['pass']:raise RuntimeError('Analytical cooling check failed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--ansys',required=True);a=p.parse_args();main(a.ansys)
