#!/usr/bin/env python3
"""Copy the same leading entries of each required forest tree for bounded validation.

This changes pPb invocation normalization; compare only runs on this same fixture.
"""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from hist_analysis.python.notebook_setup import load_root
ROOT=load_root()
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("input");p.add_argument("output",type=Path);p.add_argument("--entries",type=int,default=2000)
p.add_argument("--jet-tree",default="ak4PFJetAnalyzer/t")
a=p.parse_args()
if a.entries<1:p.error("--entries must be positive")
if a.output.exists():p.error("Output exists; choose another path")
source=ROOT.TFile.Open(a.input)
if not source or source.IsZombie():raise OSError(a.input)
paths=["hiEvtAnalyzer/HiTree","hltanalysis/HltTree","skimanalysis/HltTree",a.jet_tree]
trees=[source.Get(path) for path in paths]
if not all(trees):raise ValueError("Missing forest tree")
counts={t.GetEntries() for t in trees}
if len(counts)!=1:raise ValueError("Unaligned input trees")
n=min(a.entries,counts.pop())
output=ROOT.TFile(str(a.output),"CREATE")
for path,t in zip(paths,trees):
    directory,name=path.rsplit("/",1)
    target=output.GetDirectory(directory) or output.mkdir(directory)
    target.cd();copy=t.CloneTree(n);copy.Write()
output.Close();source.Close()
print(f"Copied {n} aligned entries to {a.output}; fixture normalization differs from the full sample")
