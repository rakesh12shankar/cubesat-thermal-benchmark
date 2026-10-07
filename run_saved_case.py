"""Run a captured APDL input in its folder, preserving existing results."""
from pathlib import Path
import argparse, re, subprocess, sys, shutil
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('--ansys',required=True);a=p.parse_args()
    run=a.case.resolve()
    if not (run/'model.inp').is_file() or not (run/'post.inp').is_file():raise SystemExit('Case must contain model.inp and post.inp.')
    if (run/'orbital.rth').exists() or (run/'orbital.lock').exists():raise SystemExit(f'Existing results protected: {run}')
    if not Path(a.ansys).is_file():raise SystemExit('MAPDL executable not found.')
    for job,inp,out in [('orbital','model.inp','solver.out'),('export','post.inp','post.out')]:
        subprocess.run([a.ansys,'-b','-smp','-np','2','-m',str(384 if len(np.load(run/'mesh.npz')['nodes'])<8000 else 1280),'-db','64','-j',job,'-i',inp,'-o',out],cwd=run,check=True)
        counts=re.findall(r'NUMBER OF ERROR\s+MESSAGES ENCOUNTERED\s*=\s*(\d+)',(run/out).read_text(errors='replace'))
        if not counts or int(counts[-1]):raise RuntimeError(f'Inspect {run/out}')
    if (run/'initial_temperatures.npy').exists() and not (run/'study_provenance.json').exists():
        path=run/'nodal_temperatures.csv'
        with path.open() as f:first_time=float(f.readline().split(',')[0])
        if first_time>0:
            initial=np.load(run/'initial_temperatures.npy')
            archived=run/'nodal_temperatures_exported.csv'
            if archived.exists():raise RuntimeError('Existing nodal export archive protected.')
            path.rename(archived)
            with path.open('w',encoding='utf-8') as f:
                for i,v in enumerate(initial,1):f.write(f'0,{i},{v:.12g}\n')
                with archived.open() as oldfile:shutil.copyfileobj(oldfile,f)
    subprocess.run([sys.executable,str(Path(__file__).with_name('analyze_study.py' if (run/'study_provenance.json').exists() else 'analyze_results.py')),str(run)],check=True)

if __name__=='__main__':main()
