#!/usr/bin/env python3
"""Launch each example in its own fresh MATLAB process (standard library)."""
import argparse,json,os,subprocess
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matlab',default='matlab',help='MATLAB executable or absolute path')
    p.add_argument('--intlab-root',default=os.environ.get('INTLAB_ROOT'),required=not bool(os.environ.get('INTLAB_ROOT')))
    p.add_argument('--case',choices=['num5','num7','pbvp','all'],default='all')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--client-threads',type=int,default=2)
    p.add_argument('--smoke',action='store_true',help='GE: one Phase2/3 interval only; PBVP: full proof')
    a=p.parse_args()
    if not 1<=a.workers<=10 or min(a.threads,a.client_threads)<1:p.error('Invalid process/thread counts')
    root=Path(__file__).resolve().parent
    for case in (['num5','num7','pbvp'] if a.case=='all' else [a.case]):
        out=a.output.resolve()/case
        if out.exists() and any(out.iterdir()):p.error(f'Output is not empty: {out}')
        out.mkdir(parents=True,exist_ok=True)
        env=os.environ.copy()
        env.update(INTLAB_ROOT=str(Path(a.intlab_root).resolve()),AI_EXAMPLE=case,AI_OUTPUT=str(out),AI_REPOSITORY=str(root))
        env['MATLABPATH']=str(root/'runtime')+os.pathsep+str(root)+os.pathsep+env.get('MATLABPATH','')
        # Paths are environment values, never interpolated as MATLAB source.
        expression=("addpath(getenv('AI_REPOSITORY'));run_case(getenv('AI_EXAMPLE'),"
                    "getenv('INTLAB_ROOT'),getenv('AI_OUTPUT'),struct("
                    f"'Workers',{a.workers},'Threads',{a.threads},'ClientThreads',{a.client_threads},'Smoke',{int(a.smoke)}))")
        subprocess.run([a.matlab,'-batch',expression],cwd=root,env=env,check=True)
        if not (out/'SUCCESS').is_file():raise RuntimeError(f'No successful certificate: {out}')
        print(json.dumps(json.loads((out/'summary.json').read_text()),indent=2))

if __name__=='__main__':main()
