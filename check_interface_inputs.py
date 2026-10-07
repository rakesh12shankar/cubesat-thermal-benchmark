"""Check captured interface cases against their mesh and fixed physical inputs."""
from pathlib import Path
import json,re
import numpy as np
ROOT=Path(__file__).resolve().parent

def main():
    lab=ROOT/'results'
    def case(label):
        for name in [f'interface_{label}_qualified',f'interface_{label}_settled',
                     f'interface_{label}',f'captured_interface_{label}']:
            folder=ROOT/'data/ansys'/name
            if (folder/'model.inp').is_file():return folder
        raise FileNotFoundError(f'No generated or captured interface case for {label}')
    base=case('bonded_control'); original=np.load(base/'mesh.npz')
    fixed=lambda s:[line for line in s.splitlines() if line.startswith(('MP,','SFE,','VFSM,','CNVTOL,','RADOPT,','HEMIOPT,','SPCTEMP,'))]
    reference=fixed((base/'model.inp').read_text(encoding='utf-8')); checks={}
    for label in ['bonded_control','seams_off','pcb_frame_off','both_off']:
        p=case(label);m=np.load(p/'mesh.npz');s=(p/'model.inp').read_text(encoding='utf-8')
        assert np.array_equal(m['nodes'][m['connectivity']],original['nodes'][original['connectivity']])
        assert fixed(s)==reference
        assert (p/'applied_flux.csv').read_bytes()==(base/'applied_flux.csv').read_bytes()
        lines=[line.split(',') for line in s.splitlines() if line.startswith('EN,')]
        assert np.array_equal(np.array([[int(n)-1 for n in line[2:]] for line in lines]),m['connectivity'])
        node_lines=[line.split(',') for line in s.splitlines() if line.startswith('N,')]
        assert [int(line[1]) for line in node_lines]==list(range(1,len(m['nodes'])+1))
        assert np.allclose(np.array([[float(v) for v in line[2:]] for line in node_lines]),m['nodes'],rtol=0,atol=1e-14)
        ics=re.findall(r'^IC,(\d+),TEMP,([^\n]+)$',s,re.M)
        assert [int(n) for n,t in ics]==list(range(1,len(m['nodes'])+1))
        assert np.allclose([float(t) for n,t in ics],np.load(p/'initial_temperatures.npy'),rtol=0,atol=1e-8)
        duration=max(float(t) for t in re.findall(r'^TIME,([^\n]+)$',s,re.M))
        flux=np.loadtxt(p/'applied_flux.csv',delimiter=',',skiprows=1)
        for side in range(1,7):
            times=[float(t) for t in re.findall(rf'^Q{side}\(\d+,0\)=([^\n]+)$',s,re.M)]
            loads=[float(t) for t in re.findall(rf'^Q{side}\(\d+,1\)=([^\n]+)$',s,re.M)]
            assert len(times)==len(loads) and np.all(np.diff(times)>0) and max(times)>=duration
            expected=np.interp(np.array(times)%5580,flux[:,0],flux[:,19+side])
            assert np.allclose(loads,expected,rtol=0,atol=1e-8)
        checks[label]={'element_geometry_unchanged':True,'materials_radiation_settings_unchanged':True,
          'input_matches_mesh':True,'all_nodes_initialized_from_captured_state':True,'periodic_loads_cover_entire_run':True}
    output=lab if lab.exists() else ROOT/'results'
    (output/'interface_input_checks.json').write_bytes((json.dumps(checks,indent=2)+'\n').encode('utf-8'))
    print('Four interface input checks passed: geometry, physics, connectivity, initialization and full load history.')

if __name__=='__main__':main()
