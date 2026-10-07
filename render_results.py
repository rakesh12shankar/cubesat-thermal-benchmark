"""Rebuild public figures and metrics from compact, included output tables."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'

def main():
    ref=np.genfromtxt(ROOT/'config/journal_figure11_reference.csv',delimiter=',',names=True)
    t=ref['time_s'];metrics={}
    files=sorted((RESULTS/'cases').glob('*/last_orbit.csv'))
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
    for row,group in enumerate([list(ref.dtype.names[1:7]),list(ref.dtype.names[7:])]):
        for col,name in enumerate(['orbital_h15_e0.5_qbase','orbital_h10_e0.5_warm']):
            data=np.genfromtxt(RESULTS/'cases'/name/'last_orbit.csv',delimiter=',',names=True)
            ax=axes[row,col]
            for j,tag in enumerate(group):
                color=f'C{j}';ax.plot(t,ref[tag],color=color,label=tag,lw=1.4)
                ax.plot(t,np.interp(t,data['phase_s'],data[tag]),color=color,ls='--',lw=1.4)
            ax.set(title='Coarse mesh, orbit 8' if col==0 else 'Fine mesh, orbit 7',ylabel='Volume-average temperature (K)',xlabel='Orbital phase (s)');ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=2)
    fig.suptitle('Internal emissivity 0.5: journal solid; spectral-absorption ANSYS baseline dashed')
    fig.tight_layout();fig.savefig(RESULTS/'temperature_comparison.png',dpi=160);plt.close(fig)
    for file in files:
        if 'e0.5' not in file.parent.name:continue # Figure 11 is only emissivity 0.5.
        data=np.genfromtxt(file,delimiter=',',names=True);case={}
        for tag in ref.dtype.names[1:]:
            delta=np.interp(t,data['phase_s'],data[tag])-ref[tag]
            case[tag]={'rmse_K':float(np.sqrt(np.mean(delta**2))),'maximum_absolute_error_K':float(abs(delta).max()),'mean_bias_K':float(delta.mean())}
        metrics[file.parent.name]=case
    (RESULTS/'recomputed_temperature_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for eps,name in [(0,'orbital_h15_e0_warm'),(.5,'orbital_h15_e0.5_qbase'),(1,'orbital_h15_e1_warm')]:
        data=np.genfromtxt(RESULTS/'cases'/name/'last_orbit.csv',delimiter=',',names=True)
        for ax,part in zip(axes,['panel4','pcb2']):ax.plot(data['phase_s'],data[part],label=f'Internal emissivity {eps:g}');ax.set(title=part,xlabel='Orbital phase (s)',ylabel='Volume-average temperature (K)');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Reconstructed model: effect of internal radiation (spectral-absorption baseline)')
    fig.tight_layout();fig.savefig(RESULTS/'emissivity_effect.png',dpi=160);plt.close(fig)
    literal=RESULTS/'cases/orbital_h15_e0.5_equation_literal/last_orbit.csv'
    if literal.exists():
        data=np.genfromtxt(literal,delimiter=',',names=True)
        base=np.genfromtxt(RESULTS/'cases/orbital_h15_e0.5_qbase/last_orbit.csv',delimiter=',',names=True)
        fig,axes=plt.subplots(1,3,figsize=(14,4.5))
        for ax,tag in zip(axes,['panel4','pcb2','battery']):
            ax.plot(t,ref[tag],label='Journal');ax.plot(base['phase_s'],base[tag],'--',label='Spectral absorption');ax.plot(data['phase_s'],data[tag],':',label='Eq. 20–24 literal');ax.set(title=tag,xlabel='Orbital phase (s)',ylabel='Temperature (K)');ax.grid(alpha=.2);ax.legend(fontsize=8)
        fig.suptitle('Boundary-condition interpretation study; no fitted parameters')
        fig.tight_layout();fig.savefig(RESULTS/'load_interpretation.png',dpi=160);plt.close(fig)
    strictfile=RESULTS/'cases/orbital_h15_e0.5_balanced_strict/last_orbit.csv'
    if strictfile.exists():
        strict=np.genfromtxt(strictfile,delimiter=',',names=True)
        base=np.genfromtxt(RESULTS/'cases/orbital_h15_e0.5_qbase/last_orbit.csv',delimiter=',',names=True)
        fig,axes=plt.subplots(1,3,figsize=(14,4.5))
        for ax,tag in zip(axes,['panel4','pcb2','battery']):
            ax.plot(t,ref[tag],label='Journal',lw=1.5);ax.plot(base['phase_s'],base[tag],'--',label='Original settings',lw=1.4);ax.plot(strict['phase_s'],strict[tag],':',label='Balanced / strict',lw=1.8)
            ax.set(title=tag,xlabel='Orbital phase (s)',ylabel='Volume-average temperature (K)');ax.grid(alpha=.2);ax.legend(fontsize=8)
        fig.suptitle('Conservation settings: physical inputs unchanged; residual paper mismatch retained')
        fig.tight_layout();fig.savefig(RESULTS/'improved_solver_comparison.png',dpi=160);plt.close(fig)
    fields=RESULTS/'final_orbit_fields.npz'
    if fields.exists():
        d=np.load(fields);xyz=d['nodes'];conn=d['connectivity'];face_nodes={1:[1,0,3,2],2:[0,1,5,4],3:[1,2,6,5],4:[2,3,7,6],5:[3,0,4,7],6:[4,5,6,7]}
        fig=plt.figure(figsize=(11,10.5));norm=matplotlib.colors.Normalize(220,360);cmap=plt.get_cmap('inferno')
        for row in range(2):
            for j,phase in enumerate(d['phase_s']):
                ax=fig.add_subplot(2,2,row*2+j+1,projection='3d');polys=[];colors=[]
                for eid,face,side,area in d['faces']:
                    e=int(eid)-1;tag=str(d['part_names'][e])
                    if row==0 and side==0:continue
                    if row==1 and tag in {'panel2','panel3','panel6'}:continue
                    ids=conn[e][face_nodes[int(face)]];polys.append(xyz[ids]*1000);colors.append(cmap(norm(d['temperature_K'][j,ids].mean())))
                ax.add_collection3d(Poly3DCollection(polys,facecolors=colors,edgecolors=(0,0,0,.12),linewidths=.12))
                ax.set(xlim=(0,100),ylim=(0,100),zlim=(0,100),xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',title=f"{'Exterior' if row==0 else 'Cutaway'}: orbital phase {phase:.0f} s")
                ax.set_box_aspect((1,1,1));ax.view_init(elev=23,azim=45 if row==0 else -55)
                ax.set_xticks([0,25,50,75]);ax.set_yticks([25,50,75,100]);ax.set_zticks([0,25,50,75,100])
        fig.suptitle('Actual ANSYS nodal fields: balanced / strict case, internal emissivity 0.5')
        fig.subplots_adjust(left=.02,right=.86,bottom=.06,top=.91,hspace=.20,wspace=.02)
        cax=fig.add_axes([.9,.2,.022,.6]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='Facet-average temperature (K)')
        fig.savefig(RESULTS/'temperature_fields.png',dpi=160);plt.close(fig)
    print(f'Recomputed journal metrics for {len(metrics)} emissivity-0.5 cases and public figures.')

if __name__=='__main__':main()
