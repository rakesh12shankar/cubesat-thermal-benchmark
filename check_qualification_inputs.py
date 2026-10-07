"""Audit captured study physics and reference-invariant mesh quantities."""
from pathlib import Path
import argparse,json,re
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
ROOT=Path(__file__).resolve().parent
NOMINAL={1:(2325,1103,1.03),2:(2810,948,140),3:(2120,975,.64),4:(2247,1110,23)}

def check(case):
    m=np.load(case/'mesh.npz');text=(case/'model.inp').read_text();p=json.loads((case/'study_provenance.json').read_text());changes=p.get('changes',{})
    props={}
    hres=changes.get('hemicube_resolution',20)
    assert re.search(rf'^HEMIOPT,{hres},1E-6$',text,re.M)
    assert re.findall(r'^VFSM,DEFINE,(\d+),([^\n]+)$',text,re.M)==[('1','2,1000,1E-8')]
    assert re.findall(r'^SPCTEMP,(\d+),([^\n]+)$',text,re.M)==[('2','2.7')]
    for kind,mid,value in re.findall(r'^MP,(DENS|C|KXX),(\d+),([^\n]+)$',text,re.M):props.setdefault(int(mid),{})[kind]=float(value)
    for mid,(rho,cp,k) in NOMINAL.items():
        assert props[mid]['DENS']==rho
        assert abs(props[mid]['C']-cp*changes.get('capacity_scale',1))<1e-8
        expected=changes.get('panel_k',k) if mid==1 else (changes.get('pcb_k',k) if mid==3 else k)
        assert abs(props[mid]['KXX']-expected)<1e-10
    vol=m['volumes'];assert abs(vol.sum()-.000247936)<1e-14
    mass=sum(vol[m['materials']==mid].sum()*v[0] for mid,v in NOMINAL.items());assert abs(mass-.57848424)<1e-9
    ext={(int(e),int(f)) for e,f,s,a in m['faces'] if s>0}
    for side in range(1,7):assert abs(sum(a for e,f,s,a in m['faces'] if s==side)-.01)<1e-12
    surfaces=re.findall(r'^SFE,(\d+),(\d+),RDSF,1,([^\n]+)$',text,re.M)
    expected_faces={(int(e),int(f)) for e,f,s,a in m['faces']}
    assert len(surfaces)==len(expected_faces) and {(int(e),int(f)) for e,f,v in surfaces}==expected_faces
    enclosures={(int(e),int(f)):int(v) for e,f,v in re.findall(r'^SFE,(\d+),(\d+),RDSF,2,(\d+)$',text,re.M)}
    assert enclosures=={key:2 if key in ext else 1 for key in expected_faces}
    heating={(int(e),int(f)):int(v) for e,f,v in re.findall(r'^SFE,(\d+),(\d+),HFLUX,1,%Q([1-6])%$',text,re.M)}
    assert heating=={(int(e),int(f)):int(s) for e,f,s,a in m['faces'] if s>0}
    for e,f,v in surfaces:
        expected=changes.get('external_epsilon',.72) if (int(e),int(f)) in ext else changes.get('internal_epsilon',.5)
        assert abs(float(v)-expected)<1e-12
    flux=np.loadtxt(case/'applied_flux.csv',delimiter=',',skiprows=1)
    expected=changes.get('alpha',.77)*(flux[:,2:8]+flux[:,8:14])+changes.get('external_epsilon',.72)*flux[:,14:20]
    assert np.max(abs(flux[:,20:26]-expected))<1e-7
    duration=max(float(v) for v in re.findall(r'^TIME,([^\n]+)$',text,re.M))
    for side in range(1,7):
        entries=re.findall(rf'^Q{side}\((\d+),([01])\)=([^\n]+)$',text,re.M)
        times={int(i):float(v) for i,col,v in entries if col=='0'}
        values={int(i):float(v) for i,col,v in entries if col=='1'}
        assert times.keys()==values.keys() and times,'Incomplete APDL heating table'
        ids=sorted(times);t=np.array([times[i] for i in ids]);q=np.array([values[i] for i in ids])
        assert ids==list(range(1,len(ids)+1)) and np.all(np.diff(t)>0)
        assert t[0]==0 and t[-1]>=duration
        actual=np.interp(t%5580,flux[:,0],flux[:,19+side])
        assert np.max(abs(q-actual))<1e-7,f'APDL forcing mismatch on exterior side {side}'
    assert not re.search(r'^(?:SF|SFE|SFA|SFL),[^\n]*\bCONV\b',text,re.M)
    assert 'cached_view_factors' not in text,'Public qualification inputs must recompute radiation'
    init=np.load(case/'initial_temperatures.npy');ics=re.findall(r'^IC,(\d+),TEMP,([^\n]+)$',text,re.M)
    assert len(ics)==len(init)==len(m['nodes'])
    assert np.array_equal([int(i) for i,t in ics],np.arange(1,len(init)+1))
    assert np.max(abs(np.array([float(t) for i,t in ics])-init))<1e-8
    return {'case':case.name,'nodes':len(init),'elements':len(vol),'mass_kg':mass,'volume_m3':float(vol.sum()),'pass':True}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('cases',nargs='+',type=Path);a=p.parse_args();print(json.dumps([check(c) for c in a.cases],indent=2))
