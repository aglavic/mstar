#!/usr/bin/env python3

from subprocess import call
from scipy.optimize import minimize
from numpy import arange

def run_mcstas(p):
    call('rm -r testrun'.split())
    call('./compile_if_needed.sh')
    call(('./mstar.out -d testrun -n 1e6 selene_b=%s selene_xi=%s'%tuple(p)).split())
    
    res=eval_mcsim('testrun/mccode.sim')
    FOM=float(res['mon3']['values'].split()[0])
    return FOM
    

def eval_mcsim(fname):
    txt=open(fname, 'r').read()
    
    output={}
    nexts=txt.find('begin data')
    while True:
        start=nexts+txt[nexts:].find('\n')
        end=start+txt[start:].find('\nend data')
        block=txt[start:end].strip()
        out=dict([map(str.strip, li.strip().split(':', 1)) for li in block.splitlines()])
        
        if not 'component' in out:
            continue
        
        output[out['component']]=out
        search=txt[end:].find('begin data')
        if search==-1:
            break
        nexts=end+search
    return output

if __name__=='__main__':
    results=[]
    
    out_p=[0.03, 0.7]
    out_I=run_mcstas(out_p)
    
    scale=1.0
    
    for i in range(15):
        did_break=False
        for db in [-0.005*scale, 0.005*scale]:
            p0=[out_p[0]+db, out_p[1]]
            res=run_mcstas(p0)
            if res>out_I:
                print('\n\nUpdate to:', p0, res, '\n')
                out_p=p0
                out_I=res
                did_break=True
                break
        if not did_break:
            scale*=0.5
            if scale<0.1:
                break
    print('\n\nFound Result:', out_p)
    #minimize(run_mcstas, p0, bounds=[(0.005, 0.1),(0.3, 0.85)])
    