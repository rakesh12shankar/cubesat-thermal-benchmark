"""Build a conformal solid-only mesh and reproducible ANSYS MAPDL inputs.

The same analytic boxes and cuts define the SolidWorks model. Coordinates in
the geometry routines are mm; MAPDL receives SI metres and Kelvin.
"""
from pathlib import Path
import argparse
import csv
import json
import math
from collections import defaultdict
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np

ROOT = Path(__file__).resolve().parent
PROPS = {1:(2325,1103,1.03),2:(2810,948,140),3:(2120,975,.64),4:(2247,1110,23)}
PCB_Z = [20,37,61,78]

def classify(x,y,z):
    if z<2: return 'panel5',1
    if z>98: return 'panel6',1
    if x<2: return 'panel1',1
    if x>98: return 'panel2',1
    if y<2: return 'panel3',1
    if y>98: return 'panel4',1
    ex=x<7 or x>93; ey=y<7 or y>93; ez=z<7 or z>93
    if sum([ex,ey,ez])>=2: return 'frame',2
    if (10<x<15 or 85<x<90) and (10<y<15 or 85<y<90): return 'bolts',2
    if 5<x<95 and 5<y<95:
        for i,lo in enumerate(PCB_Z):
            if lo<z<lo+2: return f'pcb{i+1}',3
    if 20<x<80 and 20<y<80 and 39<z<48: return 'battery',4
    return None

def subdivide(breaks,h,feature_divisions=1):
    out=[]
    for a,b in zip(breaks[:-1],breaks[1:]):
        count=max(1,math.ceil((b-a)/h),feature_divisions if b-a<=2.000001 else 1)
        out.extend(np.linspace(a,b,count+1)[:-1])
    return np.array(out+[breaks[-1]])

def make_mesh(h,feature_divisions=1):
    if not isinstance(feature_divisions,int) or feature_divisions<1:raise ValueError('Feature divisions must be a positive integer.')
    xy=subdivide([0,2,5,7,10,15,20,50,80,85,90,93,95,98,100],h,feature_divisions)
    zz=subdivide([0,2,7,20,22,37,39,48,61,63,78,80,93,98,100],h,feature_divisions)
    axes=[xy,xy,zz];cells={};nodes={};elems=[];parts=defaultdict(list)
    corners=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    for k in range(len(zz)-1):
        for j in range(len(xy)-1):
            for i in range(len(xy)-1):
                center=[(a[q]+a[q+1])/2 for a,q in zip(axes,[i,j,k])]
                tag=classify(*center)
                if tag is None: continue
                conn=[]
                for di,dj,dk in corners:
                    key=(i+di,j+dj,k+dk)
                    if key not in nodes:nodes[key]=len(nodes)+1
                    conn.append(nodes[key])
                vol=np.prod([a[q+1]-a[q] for a,q in zip(axes,[i,j,k])])*1e-9
                eid=len(elems)+1
                elems.append((eid,tag[0],tag[1],conn,vol,center))
                cells[i,j,k]=eid;parts[tag[0]].append(eid)
    # SOLID278 standard IJKL/MNOP: bottom=1, -Y=2,+X=3,+Y=4,-X=5,top=6.
    directions=[((0,0,-1),1,2,-1),((0,-1,0),2,1,-1),((1,0,0),3,0,1),
                ((0,1,0),4,1,1),((-1,0,0),5,0,-1),((0,0,1),6,2,1)]
    faces=[]
    for key,eid in cells.items():
        for delta,face,axis,sign in directions:
            neighbor=tuple(a+b for a,b in zip(key,delta))
            if neighbor in cells:continue
            coord=axes[axis][key[axis]+(sign>0)]
            external=coord==0 or coord==100
            area=np.prod([a[q+1]-a[q] for n,(a,q) in enumerate(zip(axes,key)) if n!=axis])*1e-6
            side={ (0,-1):1,(0,1):2,(1,-1):3,(1,1):4,(2,-1):5,(2,1):6}[axis,sign] if external else 0
            faces.append((eid,face,side,area))
    xyz=np.zeros((len(nodes),3))
    for key,nid in nodes.items():xyz[nid-1]=[axes[n][q]/1000 for n,q in enumerate(key)]
    return xyz,elems,parts,faces

def flux_history(period=5580,dt=10):
    """Declared beta=0 circular-orbit load reconstruction, not digitized data.

    Phase fixes ingress to the paper's t=1720 s. Earth view factors are
    integrated over the finite apparent Earth disk with Gaussian quadrature.
    """
    re=6371.;r=re+431.;theta=math.asin(re/r);n=2*math.pi/period
    gauss,weights=np.polynomial.legendre.leggauss(80)
    mu=(gauss+1)/2*(1-math.cos(theta))+math.cos(theta)
    wm=weights/2*(1-math.cos(theta))
    az=np.arange(240)*2*math.pi/240
    direction=np.stack(np.broadcast_arrays(mu[:,None],np.sqrt(1-mu[:,None]**2)*np.cos(az),np.sqrt(1-mu[:,None]**2)*np.sin(az)),axis=-1)
    normals=np.array([[-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]])
    earth=np.array([np.sum(np.maximum(direction@v,0)*wm[:,None])*(2*math.pi/240)/math.pi for v in normals])
    # Nadir is +X; Earth-to-spacecraft is -X. Ingress occurs when
    # sun-to-spacecraft angle from +X reaches theta.
    start=theta+n*1720
    times=np.unique(np.r_[np.arange(0,period+dt,dt),period])
    rows=[]
    for t in times:
        phase=start-n*t
        sun=np.array([math.cos(phase),math.sin(phase),0])
        eclipsed=math.acos(np.clip(sun[0],-1,1))<theta
        solar=1367*np.maximum(normals@sun,0)*(not eclipsed)
        albedo=1367*.3*earth*max(-sun[0],0)
        infrared=237*earth
        absorbed=.77*(solar+albedo)+.72*infrared
        rows.append([t,int(eclipsed),*solar,*albedo,*infrared,*absorbed])
    return np.array(rows),earth

def prepare(h,epsilon,duration,mode,warm_orbits=0,feature_divisions=1):
    if not (math.isfinite(h) and h>0 and math.isfinite(epsilon) and 0<=epsilon<=1):
        raise ValueError('Mesh size must be positive and emissivity must be between 0 and 1.')
    if not (math.isfinite(duration) and duration>0 and warm_orbits>=0):
        raise ValueError('Duration must be positive and warm-up orbits nonnegative.')
    if warm_orbits and duration<=warm_orbits*5580:
        raise ValueError('Duration must include a resolved cycle after warm-up.')
    suffix=f'{mode}_h{h:g}_e{epsilon:g}'+('_warm' if warm_orbits else '')
    if feature_divisions!=1:suffix+=f'_features{feature_divisions}'
    run=ROOT/'data/ansys'/suffix
    if run.exists() and any(run.iterdir()):
        raise FileExistsError(f'Existing case protected: {run}')
    run.mkdir(exist_ok=True)
    xyz,elems,parts,faces=make_mesh(h,feature_divisions)
    groups={k:[e for e in elems if e[1]==k] for k in parts}
    volumes={k:sum(e[4] for e in es) for k,es in groups.items()}
    expected={'panel5':2e-5,'panel6':2e-5,'panel1':1.92e-5,'panel2':1.92e-5,
              'panel3':1.8432e-5,'panel4':1.8432e-5,'frame':(.000025*.096*4+.000025*.086*8),
              'bolts':4*.005*.005*.096,'battery':.06*.06*.009,
              **{f'pcb{i}':(.09*.09-4*.002*.002-4*.005*.005)*.002 for i in range(1,5)}}
    for k,v in expected.items():assert abs(volumes[k]-v)<1e-12,(k,volumes[k],v)
    areas={f'panel{i}':sum(f[3] for f in faces if f[2]==i) for i in range(1,7)}
    assert all(abs(v-.01)<1e-12 for v in areas.values())
    # Every body is connected through shared nodes. Check graph, no vacuum cells.
    adjacency=defaultdict(set)
    for _,tag,_,conn,_,_ in elems:
        for nid in conn:adjacency[nid].add(tag)
    graph=defaultdict(set)
    for tags in adjacency.values():
        for a in tags:graph[a].update(tags-{a})
    seen={'battery'}
    while True:
        new=seen|set().union(*(graph[g] for g in seen))
        if new==seen:break
        seen=new
    assert seen==set(parts)
    report={'nodes':len(xyz),'solid_elements':len(elems),'internal_radiation_facets':sum(f[2]==0 for f in faces),
            'small_feature_divisions':feature_divisions,'maximum_background_spacing_mm':h,
            'external_facets':sum(f[2]>0 for f in faces),'volumes_m3':volumes,'external_area_m2':areas,
            'mass_kg':sum(e[4]*PROPS[e[2]][0] for e in elems),
            'shared_node_part_connections':{k:sorted(v) for k,v in graph.items()},
            'relative_heat_convergence_tolerance':1e-8,
            'internal_view_factor_adjustment':'closure and reciprocity, VFSM option 2' if epsilon>0 else 'no internal radiation',
            'status':'geometry partition verified; solution not yet audited'}
    (run/'mesh_audit.json').write_text(json.dumps(report,indent=2))
    np.savez_compressed(run/'mesh.npz',nodes=xyz,connectivity=np.array([e[3] for e in elems])-1,
                        volumes=np.array([e[4] for e in elems]),materials=np.array([e[2] for e in elems]),
                        part_names=np.array([e[1] for e in elems]),faces=np.array(faces))
    flux,earth=flux_history()
    np.savetxt(run/'applied_flux.csv',flux,delimiter=',',
               header=','.join(['time_s','eclipse']+[f'{k}_panel{i}_W_m2' for k in ['solar','albedo','earth_ir','absorbed'] for i in range(1,7)]),comments='')
    lines=['/BATCH','/CLEAR,START','/TITLE,CubeSat 2021 conduction and radiation','/PREP7','ET,1,SOLID278','TOFFST,0','STEF,5.670374419E-8']
    for mat,(rho,cp,k) in PROPS.items():lines.extend([f'MP,DENS,{mat},{rho}',f'MP,C,{mat},{cp}',f'MP,KXX,{mat},{k}'])
    for i,p in enumerate(xyz,1):lines.append(f'N,{i},{p[0]:.12g},{p[1]:.12g},{p[2]:.12g}')
    last=0
    for eid,tag,mat,conn,_,_ in elems:
        if mat!=last:lines.append(f'MAT,{mat}');last=mat
        lines.append('EN,'+str(eid)+','+','.join(map(str,conn)))
    for tag,ids in parts.items():
        lines.append('ESEL,NONE')
        # Individual IDs avoids assumptions about ordering between material parts.
        for eid in ids:lines.append(f'ESEL,A,ELEM,,{eid}')
        lines.extend([f'CM,{tag},ELEM','ALLSEL,ALL'])
    if mode=='orbital':
        # Extend a complete, continuous periodic table across the requested run.
        flux=np.vstack([np.c_[flux[:-1,0]+cycle*5580,flux[:-1,1:]] for cycle in range(math.ceil(duration/5580))]
                       + [np.r_[duration,flux[0,1:]][None,:]])
        for side in range(1,7):
            lines.append(f'*DIM,Q{side},TABLE,{len(flux)},1,1,TIME')
            for j,row in enumerate(flux,1):lines.extend([f'Q{side}({j},0)={row[0]:.12g}',f'Q{side}({j},1)={row[19+side]:.12g}'])
    for eid,face,side,area in faces:
        if side:
            lines.extend([f'SFE,{eid},{face},RDSF,1,.72',f'SFE,{eid},{face},RDSF,2,2'])
            q=f'%Q{side}%' if mode=='orbital' else '0'
            lines.append(f'SFE,{eid},{face},HFLUX,1,{q}')
        elif epsilon>0:lines.extend([f'SFE,{eid},{face},RDSF,1,{epsilon}',f'SFE,{eid},{face},RDSF,2,1'])
    lines.extend(['ALLSEL,ALL','FINISH'])
    if epsilon>0:lines.extend(['/AUX12','VFSM,DEFINE,1,2,1000,1E-8','FINISH'])
    lines.extend(['/SOLU','ANTYPE,TRANS','TRNOPT,FULL','THOPT,FULL','CNVTOL,HEAT,,1E-8,2,1E-6','TUNIF,300',
                  'SPCTEMP,2,2.7','HEMIOPT,20,1E-6','RADOPT,,1E-4,2,1000,1E-5,1','AUTOTS,ON',
                  'DELTIM,10,1,10','OUTRES,ALL,ALL','NEQIT,100','RESCONTROL,DEFINE,ALL,LAST','SAVE'])
    if warm_orbits:
        assert duration>warm_orbits*5580
        lines.extend(['DELTIM,60,1,60',f'TIME,{warm_orbits*5580}','SOLVE','DELTIM,10,1,10'])
    lines.extend([f'TIME,{duration}','SOLVE','FINISH',
                  '/POST1','SET,LAST','SAVE','FINISH','/EXIT,NOSAVE'])
    (run/'model.inp').write_text('\n'.join(lines)+'\n')
    post=['/BATCH','/POST1',f'FILE,{"cooldown" if mode=="cooldown" else "orbital"},RTH',
          'SET,LAST','*GET,NSETS,ACTIVE,0,SET,NSET',f'*DIM,NTEMP,ARRAY,{len(xyz)}',
          f'*DIM,NIDS,ARRAY,{len(xyz)}','*VFILL,NIDS(1),RAMP,1,1',
          '*CFOPEN,nodal_temperatures,csv','*DO,ISTEP,1,NSETS','SET,,,,,,,ISTEP',
          '*GET,TTIME,ACTIVE,0,SET,TIME','*VGET,NTEMP(1),NODE,1,TEMP',
          '*VWRITE,TTIME,NIDS(1),NTEMP(1)',"(F14.6,',',F10.0,',',E20.12)",
          '*ENDDO','*CFCLOS','FINISH','/EXIT,NOSAVE']
    # Resume DB is essential: result-file selection alone does not create nodes.
    post.insert(1,f'RESUME,{"cooldown" if mode=="cooldown" else "orbital"},DB')
    (run/'post.inp').write_text('\n'.join(post)+'\n')
    print(json.dumps({'run':str(run),'audit':report,'earth_view_factors':earth.tolist()},indent=2))
    return run

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--h',type=float,default=15)
    parser.add_argument('--epsilon',type=float,default=.5);parser.add_argument('--duration',type=float,default=100)
    parser.add_argument('--mode',choices=['cooldown','orbital'],default='cooldown')
    parser.add_argument('--warm-orbits',type=int,default=0)
    parser.add_argument('--feature-divisions',type=int,default=1)
    a=parser.parse_args();prepare(a.h,a.epsilon,a.duration,a.mode,a.warm_orbits,a.feature_divisions)
