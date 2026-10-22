#!/usr/bin/env python3
"""Conservative ROOT integrity check, shared by workers and recovery scans."""
import json
import shutil
import subprocess
import tempfile
import sys
from pathlib import Path


def validate(path):
    path = Path(path)
    try:
        if not path.is_file() or path.stat().st_size == 0:
            return False, 'missing or empty output'
    except OSError as error:
        return False, 'validation unavailable: ' + str(error)
    try:
        import ROOT
    except ImportError:
        if not shutil.which('root'):
            return False, 'validation unavailable: PyROOT or root is required'
        # The ROOT CLI uses the same library checks when PyROOT is unavailable.
        with tempfile.TemporaryDirectory(prefix='jet-output-check-') as directory:
            macro = Path(directory) / 'check.C'
            macro.write_text("""#include <TFile.h>
#include <TKey.h>
#include <TSystem.h>
void check() {
 TFile *f = TFile::Open(PATH);
 if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered)) gSystem->Exit(75);
 if (!f->Get("hVzRaw")) gSystem->Exit(75);
 TIter next(f->GetListOfKeys());
 while (TKey *key = (TKey*)next()) { if (!key->ReadObj()) gSystem->Exit(75); }
 f->Close(); gSystem->Exit(0);
}
""".replace('PATH', json.dumps(str(path))))
            try:
                result = subprocess.run(['root', '-l', '-b', '-q', str(macro)],
                                        capture_output=True, timeout=120)
            except subprocess.TimeoutExpired:
                return False, 'validation unavailable: ROOT check timed out'
            return result.returncode == 0, ('valid ROOT output' if result.returncode == 0
                                            else ('invalid ROOT output' if result.returncode == 75
                                                  else 'validation unavailable: ROOT runtime failed'))
    try:
        f = ROOT.TFile.Open(str(path))
    except OSError:
        return False, 'unreadable ROOT file'
    if not f or f.IsZombie():
        return False, 'unreadable ROOT file'
    try:
        if f.TestBit(ROOT.TFile.kRecovered):
            return False, 'recovered ROOT file'
        if not f.Get('hVzRaw'):
            return False, 'missing required hVzRaw'
        for key in f.GetListOfKeys():
            if not key.ReadObj():
                return False, 'unreadable ROOT object: ' + key.GetName()
        return True, 'valid ROOT output'
    except (OSError, ReferenceError):
        return False, 'unreadable ROOT object'
    finally:
        f.Close()


if __name__ == '__main__':
    ok, reason = validate(sys.argv[1])
    if not ok:
        print(reason, file=sys.stderr)
    sys.exit(0 if ok else 75)
