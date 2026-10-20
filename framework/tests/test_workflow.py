"""End-to-end checks using deliberately small forests with known edge cases."""
import argparse
from array import array
from pathlib import Path
import subprocess
import sys
import tempfile

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT))
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()


def forest(path, *, bad_run=False, rho="rho", count=8):
    f = ROOT.TFile(str(path), "RECREATE")
    buffers, trees = {}, []
    def tree(directory, name):
        f.mkdir(directory).cd()
        t = ROOT.TTree(name, name)
        trees.append(t)
        return t
    def scalar(t, name, code, value=0):
        pycode = {"F": "f", "D": "d", "I": "i", "i": "I", "l": "Q", "O": "b"}[code]
        b = array(pycode, [value]); buffers[name] = b
        t.Branch(name, b, f"{name}/{code}")
        return b
    def variable(t, name, code, n, values):
        b = array("f" if code == "F" else "b" if code == "O" else "i", values)
        buffers[name] = b; t.Branch(name, b, f"{name}[{n}]/{code}")
        return b
    hlt = tree("hltanalysis", "HltTree")
    for name in ("HLT_PAAK4PFJet60_Eta5p1_v4", "HLT_PAAK4PFJet80_Eta5p1_v3", "HLT_PAAK4PFJet100_Eta5p1_v3"):
        scalar(hlt, name, "I", 1)
    skim = tree("skimanalysis", "HltTree")
    for name in ("pBeamScrapingFilter", "pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "phfCoincFilter", "pVertexFilterCutdz1p0", "pVertexFilterCutGplus", "pVertexFilterCutVtx1"):
        scalar(skim, name, "I", 1)
    events = tree("hiEvtAnalyzer", "HiTree")
    scalar(events, "run", "I" if bad_run else "i", 285480)
    scalar(events, "lumi", "i", 7)
    event = scalar(events, "evt", "l", (1 << 63) + 17)
    scalar(events, "pthat", "F", 550.)
    vz = scalar(events, "vz", "F")
    weight = scalar(events, "weight", "F", 1.)
    if rho: scalar(events, rho, "F", 15.)
    jets = tree("ak4PFJetAnalyzer", "t")
    nref = scalar(jets, "nref", "I", 3)
    ngen = scalar(jets, "ngen", "I", 3)
    values = {"rawpt": [150, 110, 40], "jteta": [0.3, -0.2, 2.5], "jtphi": [0, 3.14, 1.5],
              "trackMax": [25, 20, 5], "jtPfNHF": [0.1]*3, "jtPfNEF": [0.1]*3,
              "jtPfCHF": [0.6]*3, "jtPfMUF": [0.05]*3, "jtPfCEF": [0.05]*3,
              "refpt": [155, 115, 42], "refeta": [0.3, -0.2, 2.5], "refphi": [0, 3.14, 1.5]}
    for name, v in values.items(): variable(jets, name, "F", "nref", v)
    for name in ("jtPfCHM", "jtPfCEM", "jtPfNHM", "jtPfNEM", "jtPfMUM"):
        variable(jets, name, "I", "nref", [5, 5, 5])
    for name, v in {"genpt": [155,115,42], "geneta": [0.3,-0.2,2.5], "genphi": [0,3.14,1.5]}.items():
        variable(jets, name, "F", "ngen", v)
    tracks = tree("ppTrack", "trackTree")
    scalar(tracks, "nTrk", "I", 2)
    for name,v in {"trkPt":[1.5,2.5],"trkEta":[0.2,-0.4],"trkPhi":[1.0,2.0]}.items(): variable(tracks,name,"F","nTrk",v)
    variable(tracks,"trkCharge","I","nTrk",[1,-1]);variable(tracks,"highPurity","O","nTrk",[1,0])
    particles = tree("HiGenParticleAna", "hi")
    vectors = []
    for name,t,vals in (("pt","float",[2.0,3.0]),("eta","float",[0.1,-0.2]),("phi","float",[1.0,2.0]),("chg","int",[1,-1]),("pdg","int",[211,-211]),("sube","int",[0,1])):
        v=ROOT.std.vector(t)()
        for x in vals:v.push_back(x)
        particles.Branch(name,v);vectors.append(v)
    for i in range(count):
        event[0]=(1 << 63)+17+i
        vz[0] = 16. if i==1 else -2.
        weight[0] = [1.,0.,-0.5,2.][i%4]
        nref[0] = 1 if i==2 else 3
        ngen[0] = 1 if i==3 else 3
        buffers["refpt"][0] = -1 if i==4 else 155
        buffers["refeta"][0] = 2.0 if i==5 else 0.3
        buffers["jtPfNHF"][0] = 1.0 if i==6 else 0.1
        buffers["pBeamScrapingFilter"][0] = 0 if i==7 else 1
        for t in trees: t.Fill()
    f.Write()
    for name in ("ak5PFJetAnalyzer","akCs4PFJetAnalyzer"):
        f.mkdir(name).cd();copied=jets.CloneTree(-1,"fast");copied.Write()
    f.Close()


def main():
    p = argparse.ArgumentParser();p.add_argument("--build",type=Path,required=True);args=p.parse_args()
    build=args.build.resolve()
    def run(command,ok=True):
        result=subprocess.run([str(x) for x in command],cwd=PROJECT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if (result.returncode==0)!=ok:raise AssertionError(f"{command}\n{result.stdout}")
        return result.stdout
    with tempfile.TemporaryDirectory(prefix="jetanalysis-workflow-") as d:
        d=Path(d);a=d/"HiForest_a.root";b=d/"HiForest_b.root";bad=d/"HiForest_bad.root"
        forest(a);forest(b);forest(bad,bad_run=True)
        listing=d/"inputs.list";listing.write_text(f"# input order\n{a}\n\n{b}\n")
        # Entry-by-entry equality, including IDs above 2^63 and file transitions.
        run([build/"validate_forest_reader",listing,1,"pPb"])
        cases=[(1,0,540,0,2,0,"pPb",0,0),(1,1,540,0,2,1,"pPb",1,1),
               (1,0,540,0,1,2,"pPb",2,0),(2,1,0,0,2,0,"pp",0,0),
               (0,0,0,1,2,0,"pPb",0,0),(0,1,0,2,3,1,"pPb",0,0),
               (0,1,0,3,0,2,"pp",0,0)]
        for i,options in enumerate(cases):
            ref=d/f"ref{i}.root";new=d/f"new{i}.root"
            reference=run([build/"processForestSimple",listing,ref,*options])
            candidate=run([build/"jetAnalysis",listing,new,*options])
            run([build/"compare_root_outputs",ref,new])
            # Catch event-count differences even if they happen to leave histograms unchanged.
            import re
            counts=r"Processed (\d+) events, of which (\d+) passed selection"
            assert re.search(counts,reference).groups()==re.search(counts,candidate).groups()
        run([build/"test_dijet_composition",listing,d/"task_a.root",d/"task_b.root"])
        run([build/"compare_root_outputs",d/"ref1.root",d/"task_a.root"])
        run([build/"compare_root_outputs",d/"ref1.root",d/"task_b.root"])
        out=run([build/"jetAnalysis",listing,d/"unused.root",1,0,540,0,2,0,"--reader-only","--read-tracks","--read-particles"])
        assert "tracks=32 particles=32" in out,out
        out=run([build/"jetAnalysis",listing,d/"unused.root",1,0,540,0,2,0,"--reader-only","--no-reco-jets","--no-gen-jets"])
        assert "reco=0 gen=0" in out,out
        run([build/"jetAnalysis",listing,d/"unused.root",1,0,540,0,2,0,"--no-gen-jets"],False)
        run([build/"jetAnalysis",a,d/"unused.root",1,0,540,0,2,0,"--jet-collection","AK5PF"],False)
        run([build/"jetAnalysis",a,d/"alternate.root",1,0,540,0,2,0,"--reader-only","--jet-collection","AK5PF"])
        run([build/"jetAnalysis",a,d/"alternate_cs.root",1,0,540,0,2,0,"--reader-only","--jet-collection","AK4CSPF"])
        run([build/"jetAnalysis",listing,d/"profile.root",*cases[0],"--profile-stages"])
        run([build/"compare_root_outputs",d/"ref0.root",d/"profile.root"])
        alternate=d/"HiForest_alternate.root";forest(alternate,rho="fixedGridRhoFastjetAll")
        run([build/"validate_forest_reader",alternate,1,"pPb"])
        listing.write_text(f"{a}\n{alternate}\n")
        run([build/"jetAnalysis",listing,d/"unused.root",1,0,540,0,2,0,"--reader-only"],False)
        listing.write_text(f"{a}\n{bad}\n")
        run([build/"jetAnalysis",listing,d/"unused.root",1,0,540,0,2,0,"--reader-only"],False)
        # Comparator must detect an actual discrepancy.
        f=ROOT.TFile(str(d/"new0.root"),"UPDATE");h=f.Get("hVzRaw");h.SetBinContent(1,123);h.Write("",ROOT.TObject.kOverwrite);f.Close()
        run([build/"compare_root_outputs",d/"ref0.root",d/"new0.root"],False)
        # Category exclusion must remove precisely the requested prefix.
        run([build/"jetAnalysis",a,d/"filtered.root",1,0,540,0,2,0,"--exclude-output","hVz"])
        f=ROOT.TFile(str(d/"filtered.root"));assert not any(k.GetName().startswith("hVz") for k in f.GetListOfKeys());assert f.Get("hPtHat");f.Close()
    print("PASS: seven exact physics comparisons, reader transitions, optional collections, failure policies, and output selection")
if __name__=="__main__":main()
