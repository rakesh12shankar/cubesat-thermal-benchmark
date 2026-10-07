"""Reproduce a case without overwriting existing solver results."""
from pathlib import Path
import argparse, subprocess, sys, os, shutil
from prepare_model import prepare

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--ansys',default=os.environ.get('ANSYS_EXECUTABLE') or shutil.which('ansys') or r'C:\Program Files\ANSYS Inc\v241\ansys\bin\winx64\ANSYS241.exe')
    p.add_argument('--epsilon',type=float,default=.5)
    p.add_argument('--h',type=float,default=15)
    p.add_argument('--duration',type=float,default=27900)
    p.add_argument('--warm-orbits',type=int,default=3)
    a=p.parse_args()
    if not Path(a.ansys).is_file():
        raise SystemExit('Set --ansys to a licensed MAPDL executable, or set ANSYS_EXECUTABLE.')
    root=Path(__file__).resolve().parent
    expected=root/'data/ansys'/f'orbital_h{a.h:g}_e{a.epsilon:g}_warm'
    if a.warm_orbits==0: expected=root/'data/ansys'/f'orbital_h{a.h:g}_e{a.epsilon:g}'
    if (expected/'orbital.rth').exists():
        raise SystemExit(f'Existing results protected: {expected}')
    run=prepare(a.h,a.epsilon,a.duration,'orbital',a.warm_orbits)
    for job,inp,out in [('orbital','model.inp','solver.out'),('export','post.inp','post.out')]:
        subprocess.run([a.ansys,'-b','-smp','-np','2','-j',job,'-i',inp,'-o',out],cwd=run,check=True)
        log=(run/out).read_text(errors='replace')
        if 'NUMBER OF ERROR MESSAGES ENCOUNTERED=          0' not in log:
            # ANSYS formatting varies: confirm the final error count with a regex.
            import re
            counts=re.findall(r'NUMBER OF ERROR\s+MESSAGES ENCOUNTERED\s*=\s*(\d+)',log)
            if not counts or int(counts[-1]):
                raise RuntimeError(f'Inspect solver log: {run/out}')
    subprocess.run([sys.executable,str(root/'analyze_results.py'),str(run)],check=True)

if __name__=='__main__': main()

