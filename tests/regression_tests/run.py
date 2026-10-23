#!/usr/bin/env python3
"""Build fresh baseline/candidate checkouts and run the same end-to-end cases."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys

import numpy as np
from compare import close, force, lightcone, numeric, particles, TOLERANCES
from baseline import load_baseline

HERE = Path(__file__).resolve().parent
CASES = ('gr', 'horndeski', 'kmouflage', 'lightcone')
EXAMPLES = HERE.parents[1]/'examples'/'small_simulations'
sys.path.insert(0, str(EXAMPLES))
from run_simulations import command, make_args, frontend, lua


def git(repo, *args):
    return subprocess.check_output(['git','-C',str(repo),*args],text=True).strip()


def clone(repo, ref, dest):
    # Clone through Git rather than copying the developer's working files or objects.
    command(['git','clone','--no-local','--no-checkout',repo,dest],dest.parent,dest.parent/'clone.log')
    try:
        resolved = git(dest,'rev-parse','--verify','--quiet',ref+'^{commit}')
    except subprocess.CalledProcessError:
        resolved = git(dest,'rev-parse','--verify',f'refs/remotes/origin/{ref}^{{commit}}')
    command(['git','checkout','--detach',resolved],dest,dest.parent/'checkout.log')
    command(['git','submodule','update','--init','--recursive'],dest,dest.parent/'dependencies.log')
    return dict(commit=git(dest,'rev-parse','HEAD'),submodules=git(dest,'submodule','status','--recursive'))


def build(source, args, baseline=False):
    backend,values = make_args(source,args,baseline)
    command(values,source,source.parent/'build.log',timeout=1200)
    # The probe links the application's real compiled numerical objects.
    probe_make = source.parent/'force-probe.mk'
    probe_make.write_text('force_probe: $(filter-out Main.o,$(OBJS))\n'
        f'\t$(CC) $(OPTIONS) $(INC) -Isrc {HERE / "force_probe.cpp"} -o $@ $^ $(LIB) $(LINK)\n')
    command(values+['-f','Makefile','-f',probe_make,'force_probe'],source,source.parent/'probe-build.log',timeout=1200)
    for name in ('nbody','force_probe'):
        if not (backend/name).is_file():
            raise RuntimeError(f'Missing build output: {backend/name}')
    return backend


def run_cases(source, backend, work, args):
    tables = {model:frontend(source,work,model) for model in ('horndeski','kmouflage')}
    for ranks in args.ranks:
        for case in args.cases:
            folder = work/f'{case}-np{ranks}'
            folder.mkdir()
            config = lua(source,folder,case,tables.get(case))
            launch = ['mpiexec','-n',ranks]
            command(launch+[backend/'nbody',config],source,folder/'run.log')
            log = (folder/'run.log').read_text()
            if 'Doing output 0 at redshift z = 0' not in log:
                raise ValueError(f'{case}: did not output the final epoch')
            for screened in (False,True) if case in ('horndeski','kmouflage') else (False,):
                probe_config = folder/f'force-{screened}.lua'
                text = config.read_text()+f'\ngravity_model_screening = {str(screened).lower()}\n'
                # Stress the gradient screening with a manufactured, finite coefficient.
                if screened and case=='kmouflage':
                    table = numeric(tables[case]['force']).copy()
                    table[:,1] = 1e4
                    path = folder/'screening-stress.txt'
                    np.savetxt(path,table)
                    text += f'HiCOLA_preforce_filename = {json.dumps(str(path))}\n'
                probe_config.write_text(text)
                command(launch+[backend/'force_probe',probe_config,folder/f'force-{screened}'],source,folder/f'force-{screened}.log')


def compare_runs(baseline, candidate, args):
    results = {}
    for model in ('horndeski','kmouflage'):
        for name in (('expansion','force','full') if model=='horndeski' else ('expansion','force')):
            ref = next((baseline/'frontend'/model).glob(f'{model}_*_{name}.txt'))
            got = candidate/'frontend'/model/ref.name
            results[f'{model}/frontend/{name}'] = close(numeric(got),numeric(ref),'tables')
    for ranks in args.ranks:
        for case in args.cases:
            key = f'{case}-np{ranks}'
            ref,got = baseline/key/'output',candidate/key/'output'
            reference_files = sorted(p.relative_to(ref) for p in ref.rglob('*.txt') if not p.relative_to(ref).parts[0].startswith('lightcone_'))
            actual_files = sorted(p.relative_to(got) for p in got.rglob('*.txt') if not p.relative_to(got).parts[0].startswith('lightcone_'))
            if actual_files != reference_files or not reference_files:
                raise ValueError(f'{key}: output file inventory differs or is empty')
            if not list(ref.glob('gravitymodel_*')) or not list(ref.rglob('pofk*txt')):
                raise ValueError(f'{key}: missing growth or power spectrum')
            for path in reference_files:
                kind = 'power' if 'pofk' in path.name else 'tables'
                results[f'{key}/{path}'] = close(numeric(got/path),numeric(ref/path),kind)
            ref_parts = sorted(ref.rglob('gadget_z0.000.*'))
            got_parts = sorted(got.rglob('gadget_z0.000.*'))
            if len(ref_parts)!=ranks or len(got_parts)!=ranks:
                raise ValueError(f'{key}: missing snapshot rank files')
            ids,p,v = particles(ref_parts)
            got_ids,got_p,got_v = particles(got_parts)
            np.testing.assert_array_equal(got_ids,ids)
            results[f'{key}/positions'] = close(got_p,p,'positions')
            results[f'{key}/velocities'] = close(got_v,v,'velocities')
            for screened in (False,True) if case in ('horndeski','kmouflage') else (False,):
                a = sorted((baseline/key).glob(f'force-{screened}.*.txt'))
                b = sorted((candidate/key).glob(f'force-{screened}.*.txt'))
                if len(a)!=ranks or len(b)!=ranks:
                    raise ValueError(f'{key}: missing force rank files')
                a,b = force(a),force(b)
                np.testing.assert_array_equal(a[:,:3],b[:,:3])
                results[f'{key}/force-{screened}'] = close(b[:,3:],a[:,3:],'force')
            if case=='lightcone':
                results[f'{key}/lightcone'] = close(lightcone(sorted(got.rglob('*ascii_step*'))),
                    lightcone(sorted(ref.rglob('*ascii_step*'))),'lightcone')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository',default=str(Path.cwd()))
    parser.add_argument('--candidate-ref',default='HEAD')
    parser.add_argument('--work',type=Path,required=True,help='New directory for checkouts, logs and artifacts')
    parser.add_argument('--prefix',type=Path,default=os.environ.get('CONDA_PREFIX'))
    parser.add_argument('--backend',default='HiCOLA/FML/COLASolver',help='Candidate backend directory relative to checkout')
    parser.add_argument('--fml-include',default='HiCOLA',help='Candidate directory containing FML')
    parser.add_argument('--make-arg',action='append',default=[])
    parser.add_argument('--ranks',type=int,nargs='+',default=[1,2],choices=(1,2))
    parser.add_argument('--cases',nargs='+',choices=CASES,default=list(CASES))
    args = parser.parse_args()
    if args.prefix is None:
        parser.error('Supply --prefix or set CONDA_PREFIX to the installed dependency prefix')
    args.work = args.work.resolve()
    args.work.mkdir(parents=True,exist_ok=False)
    os.environ.update(OMP_NUM_THREADS='1',MPLBACKEND='Agg')
    os.environ['PATH'] = str(args.prefix/'bin')+os.pathsep+os.environ['PATH']
    os.environ.pop('PYTHONPATH',None)
    baseline = load_baseline()
    manifest = dict(baseline_specification=baseline, candidate_ref=args.candidate_ref,
        host=socket.gethostname(),python=sys.version,cases=args.cases,ranks=args.ranks,
        tolerances=TOLERANCES,status='running',
        fixture_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted([HERE/'baseline.json', *HERE.glob('*.cpp'), *EXAMPLES.glob('*.ini'), *EXAMPLES.glob('*.lua'), EXAMPLES/'run_simulations.py'])})
    try:
        command([sys.executable,'-m','pip','freeze'],args.work,args.work/'environment.log')
        command(['mpicxx','--version'],args.work,args.work/'compiler.log')
        command(['mpiexec','--version'],args.work,args.work/'mpi.log')
        for label,ref in (('baseline',baseline['commit']),('candidate',args.candidate_ref)):
            work = args.work/label
            work.mkdir()
            source = work/'source'
            manifest[label] = clone(args.repository,ref,source)
            backend = build(source,args,baseline=label=='baseline')
            run_cases(source,backend,work,args)
        manifest['comparisons'] = compare_runs(args.work/'baseline',args.work/'candidate',args)
        manifest['status'] = 'passed'
    except Exception as error:
        manifest.update(status='failed',error=str(error))
        raise
    finally:
        (args.work/'report.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'PASS: {len(manifest["comparisons"])} comparisons; report {args.work/"report.json"}')


if __name__=='__main__':
    main()
