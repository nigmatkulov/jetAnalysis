"""Install and consume the shared libraries from a separate CMake project."""
import argparse
import importlib.util
from pathlib import Path
import subprocess
import tempfile
p=argparse.ArgumentParser();p.add_argument('--build',type=Path,required=True);a=p.parse_args()
project=Path(__file__).resolve().parents[2];build=a.build.resolve()
def run(command,cwd):
    result=subprocess.run([str(x) for x in command],cwd=cwd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'{command}\n{result.stdout}')
    if 'Error in <TCling::LoadPCM>' in result.stdout:raise RuntimeError(result.stdout)
    return result.stdout
with tempfile.TemporaryDirectory(prefix='jetanalysis-consumer-') as temporary:
    base=Path(temporary);prefix=base/'install';consumer=base/'consumer'
    run(['cmake','--install',build,'--prefix',prefix],base)
    run(['cmake','-S',project/'framework/examples','-B',consumer,f'-DCMAKE_PREFIX_PATH={prefix}','-DCMAKE_BUILD_TYPE=Release'],base)
    run(['cmake','--build',consumer,'-j','2'],base)
    run([prefix/'bin/jetAnalysis','--help'],base)
    spec=importlib.util.spec_from_file_location('workflow',project/'framework/tests/test_workflow.py')
    workflow=importlib.util.module_from_spec(spec);spec.loader.exec_module(workflow)
    source=base/'HiForest_fixture.root';workflow.forest(source)
    counts=run([consumer/'track_inventory',source],base)
    assert '8 events, 16 tracks' in counts,counts
    run([consumer/'dijet_task',source,base/'dijet.root',project/'aux_files/pPb_8160'],base)
    run([build/'processForestSimple',source,base/'reference.root',1,0,540,0,2,0],project)
    run([build/'compare_root_outputs',base/'reference.root',base/'dijet.root'],base)
print('PASS: installed core-only task, installed DiJet task, runtime library/PCM lookup, exact reference comparison')
