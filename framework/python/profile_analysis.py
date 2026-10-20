#!/usr/bin/env python3
"""Run the framework with stage timings; optionally sample its CPU stacks on macOS."""
import argparse
import json
from pathlib import Path
import subprocess
import time
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--executable",type=Path,required=True);p.add_argument("--input",required=True)
p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--sample",action="store_true")
p.add_argument("--analysis-args",nargs="+",default=["1","0","540","0","2","0"])
a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
command=[str(a.executable.resolve()),a.input,str(out/"profile.root"),*a.analysis_args,"--profile-stages"]
with (out/"analysis.log").open("w") as log:
    child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
    (out/"pid").write_text(str(child.pid))
    print(f"Profiling PID {child.pid}",flush=True)
    sample_status=None
    if a.sample:
        time.sleep(1)
        sample=subprocess.run(["sample",str(child.pid),"3","1","-file",str(out/"sample.txt")],text=True,capture_output=True)
        (out/"sample.log").write_text(sample.stdout+sample.stderr);sample_status=sample.returncode
    status=child.wait()
(out/"status.json").write_text(json.dumps(dict(command=command,analysis_status=status,sample_status=sample_status),indent=2)+"\n")
if status:raise SystemExit(status)
print((out/"analysis.log").read_text().split("STAGE",1)[-1])
if sample_status:raise SystemExit(sample_status)
