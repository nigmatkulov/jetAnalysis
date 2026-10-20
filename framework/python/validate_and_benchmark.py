#!/usr/bin/env python3
"""Run sequential reference/framework comparisons and record per-child resources."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time


def measured(command, log):
    started=time.perf_counter()
    with log.open("w") as stream:
        child=subprocess.Popen([str(x) for x in command], stdout=stream, stderr=subprocess.STDOUT)
        _, status, usage=os.wait4(child.pid,0)
        child.returncode=os.waitstatus_to_exitcode(status)
    if child.returncode:
        raise RuntimeError(f"Command failed ({child.returncode}): {command}; see {log}")
    text=log.read_text()
    match=re.search(r"Processed (\d+) events, of which (\d+) passed selection",text)
    elapsed=time.perf_counter()-started
    return dict(command=[str(x) for x in command],log=str(log),elapsed_seconds=elapsed,
                user_seconds=usage.ru_utime,system_seconds=usage.ru_stime,
                peak_rss_bytes=int(usage.ru_maxrss*(1 if sys.platform=="darwin" else 1024)),
                processed_events=int(match[1]) if match else None,
                selected_events=int(match[2]) if match else None,
                events_per_second=int(match[1])/elapsed if match else None)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--repeat",type=int,default=1)
    p.add_argument("--input",required=True)
    p.add_argument("--analysis-args",nargs="+",default=["1","0","540","0","2","0"])
    args=p.parse_args()
    if args.repeat<1:p.error("--repeat must be positive")
    build=args.build.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    report=dict(input=args.input,analysis_args=args.analysis_args,platform=sys.platform,
                note="Sequential runs; alternating order after first pair; filesystem cache and host load are uncontrolled.",runs=[])
    for i in range(args.repeat):
        results={}
        order=["reference","framework"] if i%2==0 else ["framework","reference"]
        for name in order:
            exe=build/("processForestSimple" if name=="reference" else "jetAnalysis")
            results[name]=measured([exe,args.input,out/f"{name}-{i}.root",*args.analysis_args],out/f"{name}-{i}.log")
        comparison=subprocess.run([str(build/"compare_root_outputs"),str(out/f"reference-{i}.root"),str(out/f"framework-{i}.root")],text=True,capture_output=True)
        (out/f"comparison-{i}.log").write_text(comparison.stdout+comparison.stderr)
        if comparison.returncode:raise RuntimeError(f"Numerical discrepancy: {comparison.stdout}{comparison.stderr}")
        for field in ("processed_events","selected_events"):
            if results["reference"][field]!=results["framework"][field]:raise RuntimeError(f"Different {field}")
        results["comparison"]=comparison.stdout.strip()
        report["runs"].append(results)
        (out/"report.json").write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps(results,indent=2),flush=True)
if __name__=="__main__":main()
