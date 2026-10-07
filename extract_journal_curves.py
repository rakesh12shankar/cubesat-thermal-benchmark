"""Extract Figures 9/11 from the verified user-supplied journal PDF.

The PDF is not redistributed. Vector path indices are specific to its hash.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pymupdf
ROOT=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--pdf',type=Path,required=True);p.add_argument('--output',type=Path,default=Path('extracted_reference'));a=p.parse_args()
    digest=hashlib.sha256(a.pdf.read_bytes()).hexdigest()
    expected=json.loads((ROOT/'config/source_hashes.json').read_text(encoding='utf-8'))['Morsch_2021_Journal.pdf']
    if digest!=expected:raise SystemExit('Different PDF bytes: inspect figures and recalibrate vector paths before using this extractor.')
    out=a.output.resolve()
    if out.exists() and any(out.iterdir()):raise SystemExit(f'Existing extraction protected: {out}')
    out.mkdir(parents=True,exist_ok=True);doc=pymupdf.open(a.pdf);t=np.arange(10,5581,10)
    def curve(page,index,x0,x5000,y0,ytop,v0,vtop):
        points=[]
        for item in doc[page].get_drawings()[index]['items']:
            if item[0]=='l':points.extend([tuple(item[1]),tuple(item[2])])
        b=np.array(points);b=b[np.argsort(b[:,0],kind='stable')]
        x,idx=np.unique(b[:,0],return_index=True);y=b[idx,1]
        tt=(x-x0)*5000/(x5000-x0);v=v0+(y-y0)*(vtop-v0)/(ytop-y0)
        assert tt.min()<11 and tt.max()>5579
        return np.interp(t,tt,v)
    temps={f'panel{i}':curve(8,106+i,106.005,274.286,206.756,57.800,230,350) for i in range(1,7)}
    for tag,index in [('pcb4',249),('pcb3',250),('pcb2',251),('pcb1',252),('battery',253)]:temps[tag]=curve(8,index,318.930,487.210,206.756,57.800,260,285)
    flux={f'panel{i}':curve(6,141+i,339.959,527.325,418.8959,258.072,0,1600) for i in range(1,7)}
    for name,data in [('journal_figure11_reference.csv',temps),('journal_figure9_reference.csv',flux)]:
        np.savetxt(out/name,np.column_stack([t,*data.values()]),delimiter=',',header=','.join(['time_s',*data]),comments='')
    print(f'Extracted source curves from verified PDF into {out}')

if __name__=='__main__':main()
