"""Extend a contact diagnostic with six coarse warm-up and two resolved orbits."""
from collections import deque
import json,re,shutil,sys
import numpy as np
from study_interface_limits import ROOT,run

def settle(label,resolved_only=False):
    base=ROOT/'data/ansys'/f"interface_{label}{'_settled' if resolved_only else ''}"
    dest=ROOT/'data/ansys'/f"interface_{label}{'_qualified' if resolved_only else '_settled'}";dest.mkdir(exist_ok=False)
    mesh=np.load(base/'mesh.npz');n=len(mesh['nodes'])
    with (base/'nodal_temperatures.csv').open() as f:tail=list(deque(f,maxlen=n))
    rows=np.array([[float(v) for v in line.split(',')] for line in tail])
    assert np.all(rows[:,0]==(44640 if resolved_only else 11160)) and np.all(rows[:,1]==np.arange(1,n+1))
    initial=rows[:,2];np.save(dest/'initial_temperatures.npy',initial)
    for name in ['mesh.npz','mesh_audit.json','applied_flux.csv']:shutil.copy2(base/name,dest/name)
    text=(base/'model.inp').read_text(encoding='utf-8')
    text=re.sub(r'^IC,.*$',lambda m:'',text,flags=re.M)
    text=text.replace('SPCTEMP,2,2.7','\n'.join(f'IC,{i},TEMP,{t:.12g}' for i,t in enumerate(initial,1))+'\nSPCTEMP,2,2.7')
    if resolved_only:
        old='DELTIM,60,1,60\nTIME,33480\nSOLVE\nDELTIM,10,1,10\nTIME,44640\nSOLVE'
        assert old in text
        text=text.replace(old,'DELTIM,10,1,10\nTIME,11160\nSOLVE')
        (dest/'model.inp').write_text(text,encoding='utf-8');shutil.copy2(base/'post.inp',dest/'post.inp')
        p=json.loads((base/'interface_provenance.json').read_text(encoding='utf-8'))
        p.update(duration_s=11160,initial_state_source=base.name+' final phase-zero state',extension='Two additional fully resolved 10 s orbits after the six warm-up and two resolved-orbit extension; comparison uses final orbit only.')
        (dest/'interface_provenance.json').write_text(json.dumps(p,indent=2),encoding='utf-8')
        print('Resolving '+label,flush=True);run(dest,initial);return
    text=re.sub(r'^(?:\*DIM,Q[1-6],TABLE.*|Q[1-6]\(.*)\n','',text,flags=re.M)
    flux=np.loadtxt(base/'applied_flux.csv',delimiter=',',skiprows=1)
    extended=np.vstack([np.c_[flux[:-1,0]+cycle*5580,flux[:-1,1:]] for cycle in range(8)]+[np.r_[44640,flux[0,1:]][None,:]])
    table=[]
    for side in range(1,7):
        table.append(f'*DIM,Q{side},TABLE,{len(extended)},1,1,TIME')
        for j,row in enumerate(extended,1):table.extend([f'Q{side}({j},0)={row[0]:.12g}',f'Q{side}({j},1)={row[19+side]:.12g}'])
    text=text.replace('SFE,','\n'.join(table)+'\nSFE,',1)
    assert 'TIME,11160\nSOLVE' in text
    text=text.replace('TIME,11160\nSOLVE','DELTIM,60,1,60\nTIME,33480\nSOLVE\nDELTIM,10,1,10\nTIME,44640\nSOLVE')
    (dest/'model.inp').write_text(text,encoding='utf-8');shutil.copy2(base/'post.inp',dest/'post.inp')
    p=json.loads((base/'interface_provenance.json').read_text(encoding='utf-8'))
    p.update(duration_s=44640,initial_state_source=base.name+' final phase-zero state',extension='Six 60 s warm-up orbits followed by two 10 s orbits; comparison uses final 10 s orbit only.')
    (dest/'interface_provenance.json').write_text(json.dumps(p,indent=2),encoding='utf-8')
    print('Extending '+label,flush=True);run(dest,initial)

if __name__=='__main__':settle(sys.argv[1],resolved_only='--resolved-only' in sys.argv)
