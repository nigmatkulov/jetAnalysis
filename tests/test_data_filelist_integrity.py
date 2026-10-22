import contextlib
import io
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

PROCESSING = Path(__file__).resolve().parents[1] / 'processing'
sys.path.insert(0, str(PROCESSING))
import check_data_filelists as checker


class FileListIntegrityTest(unittest.TestCase):
    def test_non_utf8_root_diagnostics_preserve_result(self):
        real_run = subprocess.run
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'a.root'
            path.write_text('fixture')
            for code in (0, 12):
                def emit_diagnostics(command, **kwargs):
                    return real_run([
                        sys.executable, '-c',
                        'import sys; '
                        'sys.stdout.buffer.write(b"diagnostic \\xf5\\nROOT_CHECK_RESULT=' + str(code) + '\\n"); '
                        'sys.stderr.buffer.write(b"ROOT error \\xf5\\n"); '
                        'sys.exit(' + str(code) + ')',
                    ], **kwargs)
                with self.subTest(code=code), patch.object(checker.subprocess, 'run', side_effect=emit_diagnostics):
                    issue = checker.check_root_file(path)
                    if code == 0:
                        self.assertIsNone(issue)
                    else:
                        self.assertEqual(issue.category, 'CORRUPT')

    def test_worker_exception_is_reported_and_scan_continues(self):
        paths = [Path('bad.root'), Path('good.root')]
        def validate(path, *args):
            if path == paths[0]:
                raise UnicodeDecodeError('utf-8', b'\xf5', 0, 1, 'invalid start byte')
            return None
        with patch.object(checker, 'check_root_file', side_effect=validate):
            progress = checker.ProgressBar(2)
            good, missing, issues = checker.check_integrity_paths(paths, 2, progress)
        self.assertEqual(good, (paths[1],))
        self.assertEqual(missing, ())
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].path, paths[0])
        self.assertEqual(issues[0].category, 'CHECK_ERROR')
        self.assertIn('UnicodeDecodeError', issues[0].detail)
        self.assertEqual(progress.completed, 2)

    def test_url_and_metadata_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            filelist = Path(directory)/'files.txt'
            filelist.write_text('# comment\n\nroot://host//data/file.root 12\n/tmp/a.root metadata\n')
            self.assertEqual(checker.read_root_paths(filelist), ['root://host//data/file.root', Path('/tmp/a.root')])

    def test_timeout_and_root_crash_are_not_corruption(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'a.root'
            path.write_text('fixture')
            with patch.object(checker.subprocess, 'run', side_effect=subprocess.TimeoutExpired('root', 1)):
                self.assertEqual(checker.check_root_file(path, timeout=1).category, 'TIMEOUT')
            with patch.object(checker.subprocess, 'run', return_value=subprocess.CompletedProcess([], -11, '', 'crash')):
                self.assertEqual(checker.check_root_file(path).category, 'CHECK_ERROR')
            with patch.object(checker.subprocess, 'run', return_value=subprocess.CompletedProcess([], 11, 'ROOT_CHECK_RESULT=11\n', '')):
                self.assertEqual(checker.check_root_file(path).category, 'RECOVERED')

    def test_parallel_results_preserve_order(self):
        paths = [Path('a.root'), Path('b.root'), Path('c.root')]
        def validate(path, *args):
            return None if path == paths[1] else checker.FileProblem(path, 'CORRUPT', 'bad object')
        with patch.object(checker, 'check_root_file', side_effect=validate):
            progress = checker.ProgressBar(3)
            good, missing, issues = checker.check_integrity_paths(paths, 3, progress)
        self.assertEqual(good, (paths[1],))
        self.assertEqual(missing, ())
        self.assertEqual([issue.path for issue in issues], [paths[0], paths[2]])
        self.assertEqual(progress.completed, 3)

    def test_availability_only_lists_problems_without_root(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            good, missing, filelist = base/'good.root', base/'missing.root', base/'files.txt'
            good.write_text('availability only does not inspect contents')
            filelist.write_text(str(good) + '\n' + str(missing) + '\n')
            output = io.StringIO()
            with patch.object(sys, 'argv', ['check']), \
                 patch.object(checker, 'expected_filelists', return_value=[checker.FileListSpec('fixture', filelist)]), \
                 patch.object(checker.shutil, 'which', return_value=None), contextlib.redirect_stdout(output):
                self.assertEqual(checker.main(), 1)
            self.assertIn('availability only', output.getvalue())
            self.assertIn('MISSING: ' + str(missing), output.getvalue())

    def test_missing_root_runtime(self):
        with patch.object(sys, 'argv', ['check', 'a.root']), patch.object(checker.shutil, 'which', return_value=None):
            self.assertEqual(checker.main(), 2)

    @unittest.skipUnless(shutil.which('root'), 'ROOT CLI unavailable')
    def test_real_files_and_final_problem_report(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            good = base/'good.root'
            macro = base/'make.C'
            macro.write_text('#include <TFile.h>\n#include <TH1D.h>\n#include <TTree.h>\n'
                             'void make(){TFile f("' + str(good) + '","RECREATE");'
                             'f.mkdir("nested")->cd(); TH1D h("h","",10,0,10); h.Fill(1); h.Write();'
                             'TTree t("events",""); int n=0; t.Branch("n",&n);'
                             'for(n=0;n<100;n++) t.Fill(); t.Write(); f.Close();}\n')
            subprocess.run(['root','-l','-b','-q',str(macro)], check=True, capture_output=True)
            self.assertIsNone(checker.check_root_file(good))
            self.assertIsNone(checker.check_root_file(good, deep=True))
            zombie = base/'zombie.root'
            zombie.write_text('not ROOT')
            self.assertEqual(checker.check_root_file(zombie).category, 'ZOMBIE')
            empty = base/'empty.root'
            empty.touch()
            self.assertEqual(checker.check_root_file(empty).category, 'EMPTY')
            truncated = base/'truncated.root'
            content = good.read_bytes()
            truncated.write_bytes(content[:len(content)//2])
            self.assertIsNotNone(checker.check_root_file(truncated, deep=True))
            missing = base/'missing.root'
            output = io.StringIO()
            with patch.object(sys, 'argv', ['check', str(good), str(zombie), str(empty), str(missing)]), contextlib.redirect_stdout(output):
                self.assertEqual(checker.main(), 1)
            final = output.getvalue().split('Problem files:')[1]
            for path in (zombie, empty, missing):
                self.assertIn(str(path), final)
            self.assertNotIn(str(good), final)
            with patch.object(sys, 'argv', ['check', str(good), '--deep']), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(checker.main(), 0)


if __name__ == '__main__':
    unittest.main()
