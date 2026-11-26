"""ROOT canvas names must not invalidate plots retained for later archiving."""
from pathlib import Path
import subprocess
import sys
import unittest


class CanvasLifetimeTest(unittest.TestCase):
    def test_repeated_plots_survive_and_archive_in_subprocess(self):
        # Isolate native lifetime failures so a regression cannot kill the suite.
        code = '''
from pathlib import Path
import tempfile
from hist_analysis.python.notebook_setup import load_root
ROOT = load_root()
from hist_analysis.python.plotting import draw_overlay, draw_closure
from hist_analysis.python.notebook_output import NotebookRootOutput
ROOT.TH1.AddDirectory(False)
canvases = []
for renderer in (draw_overlay, draw_closure):
    for explicit in (None, 'same_logical_name'):
        for index in range(2):
            hist = ROOT.TH1D(f'h_{len(canvases)}', '', 2, 0, 2)
            hist.SetDirectory(0)
            hist.Fill(.5, 2)
            kwargs = dict(title='', x_title='eta', y_title='yield', canvas_name=explicit)
            if renderer is draw_closure:
                kwargs['nominal'] = 'Gen'
            result = renderer({'Gen': hist}, **kwargs)
            canvases.append(result[0] if isinstance(result, tuple) else result)
assert len({canvas.GetName() for canvas in canvases}) == 8
assert ROOT.gROOT.GetListOfCanvases().GetSize() == 8
with tempfile.TemporaryDirectory() as directory:
    archive = NotebookRootOutput(Path(directory), 'canvas_lifetime')
    archive.finish({'canvases': canvases})
    output = ROOT.TFile.Open(str(archive.path))
    assert output.Get('notebookConfiguration')
    assert len(archive.manifest) >= 9
    output.Close()
for canvas in canvases:
    canvas.Update()
print('PASS all retained canvases survive and archive')
'''
        result = subprocess.run([sys.executable, '-c', code], cwd=Path(__file__).resolve().parents[1],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, (result.stdout + result.stderr)[-6000:])
        self.assertNotIn('Deleting canvas with same name', result.stderr)
        self.assertIn('PASS all retained canvases survive and archive', result.stdout)


if __name__ == '__main__':
    unittest.main()
