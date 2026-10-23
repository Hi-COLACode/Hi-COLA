#!/usr/bin/env python3
"""Run small GR, Horndeski, K-mouflage and lightcone simulations."""
import argparse
import json
import os
from pathlib import Path
import shlex
import signal
import socket
import subprocess
import sys
import time

import numpy as np
from configobj import ConfigObj

HERE = Path(__file__).resolve().parent
CASES = ('gr', 'horndeski', 'kmouflage', 'lightcone')


def command(argv, cwd, log, timeout=900):
    """Keep commands and output after failures; kill MPI children on timeout."""
    argv = [str(x) for x in argv]
    log.parent.mkdir(parents=True,exist_ok=True)
    start = time.monotonic()
    with log.open('w') as stream:
        with subprocess.Popen(argv,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT,
                              start_new_session=True) as process:
            info = dict(command=argv,cwd=str(cwd),host=socket.gethostname(),pid=process.pid,
                        log=str(log),stop=f'kill -TERM -- -{process.pid}')
            print(json.dumps(info),flush=True)
            log.with_suffix('.command.json').write_text(json.dumps(info,indent=2)+'\n')
            try:
                code = process.wait(timeout=timeout)
            except (subprocess.TimeoutExpired,KeyboardInterrupt):
                os.killpg(process.pid,signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL)
                    process.wait()
                raise
    info.update(returncode=code,elapsed_seconds=time.monotonic()-start)
    log.with_suffix('.command.json').write_text(json.dumps(info,indent=2)+'\n')
    if code:
        raise RuntimeError(f'Exit {code}: {shlex.join(argv)}; see {log}')


def make_args(source, args, baseline=False):
    prefix = args.prefix
    backend = source/('HiCOLA/FML/COLASolver' if baseline else args.backend)
    include = source/('HiCOLA' if baseline else args.fml_include)
    values = ['make','-C',backend,'-j2',f'FML_INCLUDE={include}','USE_MPI=true','USE_OMP=false']
    for library in ('FFTW','GSL','LUA'):
        values += [f'{library}_INCLUDE={prefix}/include',f'{library}_LIB={prefix}/lib']
    values += args.make_arg
    return backend,values


def frontend(source, work, model):
    out = work/'frontend'/model
    out.mkdir(parents=True)
    if model == 'kmouflage':
        # Analytic LCDM expansion and nonzero coupling exercise gradient screening.
        a = np.geomspace(1e-4,1.,200)
        r, m = 4.28e-5/.7307**2, (.02196+.1274)/.7307**2
        E2 = r/a**4+m/a**3+1-r-m
        expansion = np.column_stack((a,np.sqrt(E2),-.5*(4*r/a**4+3*m/a**3)/E2))
        data = {'expansion':expansion,'force':np.column_stack((a,np.full_like(a,1e4),np.full_like(a,.02)))}
        paths = {}
        for name,values in data.items():
            paths[name] = out/f'{model}_example_{name}.txt'
            np.savetxt(paths[name],values)
        return paths
    theory = ConfigObj(str(HERE/f'{model}.ini'))
    theory['output_directory'] = str(out)
    theory['model_name'] = model
    theory.filename = str(out/'model.ini')
    theory.write()
    config = ConfigObj(str(HERE/'numerical.ini'))
    config.update(dict(Npoints='400',max_redshift='9999',suppression_flag='False',forwards_flag='False'))
    config.filename = str(out/'numerical.ini')
    config.write()
    module = 'generate_simulation_input'
    command([sys.executable,'-m',f'HiCOLA.Frontend.{module}',theory.filename,config.filename],source,out/'run.log')
    paths = {}
    for name in ('expansion','force','full'):
        matches = list(out.glob(f'{model}_*_{name}.txt'))
        if len(matches) != 1:
            raise RuntimeError(f'Expected one {model} {name} table; got {matches}')
        data = np.loadtxt(matches[0], ndmin=2)
        if len(data) != 400 or not np.isfinite(data).all() or not (np.all(np.diff(data[:,0])>0) or np.all(np.diff(data[:,0])<0)):
            raise ValueError(f'Invalid frontend scale-factor grid: {matches[0]}')
        paths[name] = matches[0]
    return paths


def lua(source, folder, case, tables):
    """Combine shared settings with the parameters for this simulation."""
    text = (HERE/'backend_minimal.lua').read_text()
    # Append model-specific settings and paths after the shared fixture.
    values = dict(all_parameters_must_be_in_file=False,simulation_name=case,
        simulation_boxsize=256., particle_Npart_1D=16, particle_allocation_factor=4.,
        ic_nmesh=16, force_nmesh=32, pofk_nmesh=32, output_particles=True,
        output_folder=str(folder/'output'),output_fileformat='GADGET',
        ic_random_seed=1234,ic_fix_amplitude=True,ic_initial_redshift=49.,
        force_linear_massive_neutrinos=False,cosmology_model='LCDM',gravity_model='GR',
        cosmology_h=.7307,cosmology_Omegab=.02196/.7307**2,
        cosmology_OmegaCDM=.1274/.7307**2,cosmology_OmegaMNu=0.,
        ic_input_filename=str(source/'HiCOLA/Parameter_files/Backend_parameter_files/LCDM_1306dot3219_matterpower_z0.000.dat'),
        lightcone=case=='lightcone',plc_pos_observer=[.5,.5,.5],
        plc_z_init=.05,plc_z_finish=0.,plc_use_replicas=False,plc_boundary_mpch=0.,
        plc_output_asciifile=True,plc_output_gadgetfile=False,plc_output_in_batches=False,
        plc_make_onion_density_maps=False,fof=False,pofk=True,
        output_redshifts=[0.],timestep_nsteps=[5])
    if case in ('horndeski','kmouflage'):
        values.update(cosmology_model='HiCOLA',gravity_model='HiCOLA' if case=='horndeski' else 'kmouflage',
            HiCOLA_expansion_filename=str(tables['expansion']),HiCOLA_preforce_filename=str(tables['force']),
            gravity_model_screening=True,gravity_model_HiCOLA_smoothing_filter='gaussian',
            gravity_model_HiCOLA_smoothing_scale_over_boxsize=3./256.,
            gravity_model_screening_enforce_largescale_linear=True,gravity_model_screening_linear_scale_hmpc=.1)
    def literal(value):
        if isinstance(value,bool):
            return str(value).lower()
        if isinstance(value,list):
            return '{'+','.join(map(literal,value))+'}'
        return json.dumps(value)
    text += '\n-- Simulation-specific settings.\n'
    text += '\n'.join(f'{key} = {literal(value)}' for key,value in values.items())+'\n'
    path = folder/'run.lua'
    path.write_text(text)
    return path



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=HERE.parents[1])
    parser.add_argument('--work', type=Path, required=True, help='New directory for logs, inputs and outputs')
    parser.add_argument('--prefix', type=Path, default=os.environ.get('CONDA_PREFIX', '/usr'))
    parser.add_argument('--backend', default='HiCOLA/FML/COLASolver')
    parser.add_argument('--fml-include', default='HiCOLA')
    parser.add_argument('--make-arg', action='append', default=[])
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    args = parser.parse_args()
    if args.ranks < 1:
        parser.error('--ranks must be positive')
    source = args.source.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    os.environ.setdefault('MPLBACKEND', 'Agg')
    backend, values = make_args(source, args)
    command(values, source, work/'build.log', timeout=1200)
    tables = {model: frontend(source, work, model)
              for model in ('horndeski', 'kmouflage') if model in args.cases}
    for case in args.cases:
        folder = work/case
        folder.mkdir()
        config = lua(source, folder, case, tables.get(case))
        command(['mpiexec', '-n', args.ranks, backend/'nbody', config], source, folder/'run.log')
        if 'Doing output 0 at redshift z = 0' not in (folder/'run.log').read_text():
            raise RuntimeError(f'{case}: missing final output')
    print(f'Completed {len(args.cases)} simulations; outputs: {work}')


if __name__ == '__main__':
    main()
