import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

PROCESSING = Path(__file__).resolve().parents[1] / 'processing'
sys.path.insert(0, str(PROCESSING))
import retry_failed_io_jobs as recovery
import validate_job_output as output_check


def event(code, body='', job='10.0'):
    return f'{code} ({job}.000) 10/21 01:00:00 {body}\n'


class RecoveryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.row = recovery.QueueRow('sample_job0001', self.base/'input.list',
                                     self.base/'out.root', '')
        self.row.input_list.write_text('a.root\nb.root\n')

    def classify(self, text, err='', valid=True, active=()):
        return recovery.classify(self.row, text, err, set(active),
                                 lambda _: (valid, 'valid ROOT output' if valid else 'invalid ROOT output'))[0]

    def test_failure_categories(self):
        success = event('005', 'Job terminated.\nNormal termination (return value 0)')
        self.assertEqual(self.classify(success), 'complete')
        self.assertEqual(self.classify(success, valid=False), 'invalid_output')
        self.assertEqual(self.classify(success, 'TLS error: connection lost'), 'transient')
        wall = event('009', 'Job removed by SYSTEM_PERIODIC_REMOVE due to wall time exceeded allowed max.')
        self.assertEqual(self.classify(wall), 'walltime')
        self.assertEqual(self.classify(wall, active=['10.0']), 'active')
        self.assertEqual(self.classify(success + event('000', 'submitted')), 'active')
        self.assertEqual(self.classify(event('012', 'held')), 'held')
        self.assertEqual(self.classify(event('005', 'Normal termination (return value 1)')), 'manual')
        self.assertEqual(self.classify(''), 'manual')
        self.assertEqual(self.classify(event('005', 'Normal termination (return value 1)'), 'failed to transfer output'), 'transfer_failure')
        self.assertEqual(self.classify(event('012', 'failed to transfer output')), 'held_transfer_failure')

    def campaign(self, termination):
        (self.base/'log').mkdir()
        (self.base/'log/sample_job0001.log').write_text(termination)
        (self.base/'log/sample_job0001.err').write_text('TLS error: test\n')
        sub = self.base/'sample.sub'
        sub.write_text('executable = /repo/processing/run_process_forest_condor.sh\n'
                       '+JobFlavour = "microcentury"\n' + recovery.QUEUE_HEADER + '\n'
                       f'{self.row.job_name} {self.row.input_list} {self.row.output_file}\n)\n')
        return sub

    def run_main(self, *options, queue=None, validator=None):
        def run(command, **kwargs):
            if command[0] == 'condor_q':
                return subprocess.CompletedProcess(command, 0, json.dumps(queue or []), '')
            return subprocess.CompletedProcess(command, 0, '20.0 - 20.0\n', '')
        with patch.object(sys, 'argv', ['retry', str(self.base), *options]), \
             patch.object(recovery.shutil, 'which', return_value='/bin/mock'), \
             patch.object(recovery.subprocess, 'run', side_effect=run), \
             patch.object(output_check, 'validate', side_effect=validator or (lambda _: (False, 'missing output'))), \
             contextlib.redirect_stdout(io.StringIO()):
            return recovery.main()

    def test_recovery_end_to_end_and_idempotence(self):
        sub = self.campaign(event('005', 'Normal termination (return value 74)'))
        self.run_main('--submit')
        journal = json.loads(sub.with_suffix('.recovery.json').read_text())
        self.assertEqual(len(journal[self.row.job_name]), 1)
        self.run_main('--submit')  # scheduler log not yet delivered: no duplicate
        self.assertEqual(len(list(self.base.glob('*_retry_*.sub'))), 1)
        (self.base/'log/sample_job0001.log').write_text(event('005', 'Normal termination (return value 0)', '20.0'))
        (self.base/'log/sample_job0001.err').write_text('JETANALYSIS_COMPLETED\n')
        self.assertEqual(self.run_main('--submit', validator=lambda _: (True, 'valid ROOT output')), 0)
        self.assertEqual(len(list(self.base.glob('*_retry_*.sub'))), 1)
        _, _, _, rows = recovery.parse_queue(next(self.base.glob('*_retry_*.sub')))
        self.assertEqual(rows[0].input_list.read_text(), 'a.root\nb.root\n')
        self.assertEqual(rows[0].output_file, self.row.output_file)

    def test_walltime_requires_strictly_longer_flavour(self):
        self.campaign(event('009', 'wall time exceeded allowed max'))
        self.run_main('--submit')
        self.run_main('--submit', '--walltime-flavour', 'microcentury')
        self.assertFalse(list(self.base.glob('*_retry_*.sub')))
        self.run_main('--submit', '--walltime-flavour', 'longlunch')
        self.assertIn('+JobFlavour = "longlunch"', next(self.base.glob('*_retry_*.sub')).read_text())

    def test_longer_flavour_survives_transient_retry(self):
        self.campaign(event('009', 'wall time exceeded allowed max'))
        self.run_main('--submit', '--walltime-flavour', 'longlunch')
        (self.base/'log/sample_job0001.log').write_text(event('005', 'Normal termination (return value 74)', '20.0'))
        (self.base/'log/sample_job0001.err').write_text('TLS error: failed')
        self.run_main('--submit')
        retries = list(self.base.glob('*_retry_*.sub'))
        self.assertEqual(len(retries), 2)
        self.assertTrue(all('+JobFlavour = "longlunch"' in path.read_text() for path in retries))

    def test_limit_and_active_queue_prevent_submission(self):
        self.campaign(event('005', 'Normal termination (return value 74)'))
        self.run_main('--submit', '--max-retries', '0')
        self.run_main('--submit', queue=[dict(ClusterId=10, ProcId=0, JobStatus=2)])
        self.assertFalse(list(self.base.glob('*_retry_*.sub')))

    def test_uncertain_submission_blocks_retry(self):
        sub = self.campaign(event('005', 'Normal termination (return value 74)'))
        sub.with_suffix('.recovery.json').write_text(json.dumps({self.row.job_name: [dict(state='uncertain')]}))
        self.run_main('--submit')
        self.assertFalse(list(self.base.glob('*_retry_*.sub')))

    def test_duplicate_output_campaigns_are_blocked(self):
        sub = self.campaign(event('005', 'Normal termination (return value 74)'))
        (self.base/'second.sub').write_text(sub.read_text())
        self.run_main('--submit')
        self.assertFalse(list(self.base.glob('*_retry_*.sub')))

    def test_changed_manifest_is_blocked(self):
        sub = self.campaign(event('005', 'Normal termination (return value 74)'))
        sub.with_suffix('.manifest.json').write_text(json.dumps(dict(submit_sha256='changed', jobs={})))
        self.run_main('--submit')
        self.assertFalse(list(self.base.glob('*_retry_*.sub')))

    def test_repeated_failure_exhausts_limit(self):
        sub = self.campaign(event('005', 'Normal termination (return value 74)'))
        self.run_main('--submit', '--max-retries', '1')
        (self.base/'log/sample_job0001.log').write_text(event('005', 'Normal termination (return value 74)', '20.0'))
        (self.base/'log/sample_job0001.err').write_text('TLS error: persistent')
        self.run_main('--submit', '--max-retries', '1')
        self.assertEqual(len(json.loads(sub.with_suffix('.recovery.json').read_text())[self.row.job_name]), 1)

    def test_missing_and_corrupt_root(self):
        self.assertFalse(output_check.validate(self.row.output_file)[0])
        self.row.output_file.write_text('not a ROOT file')
        self.assertFalse(output_check.validate(self.row.output_file)[0])


class RootOutputIntegrationTest(unittest.TestCase):
    @unittest.skipUnless(__import__('shutil').which('root'), 'ROOT CLI unavailable')
    def test_real_root_and_worker_success(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            macro = directory/'make.C'
            macro.write_text('#include <TFile.h>\n#include <TH1D.h>\n'
                             'void make(){ TFile f("' + str(directory/'out.root') + '","RECREATE");'
                             'TH1D h("hVzRaw","",10,0,10); h.Fill(1); h.Write(); f.Close(); }\n')
            subprocess.run(['root', '-l', '-b', '-q', str(macro)], check=True, capture_output=True)
            self.assertTrue(output_check.validate(directory/'out.root')[0])
            executable = directory/'analysis'
            executable.write_text('#!/bin/bash\ncp "' + str(directory/'out.root') + '" "$2"\n')
            executable.chmod(0o755)
            inputs = directory/'input.list'
            inputs.write_text('a.root\n')
            command = [str(PROCESSING/'run_process_forest_condor.sh'), str(PROCESSING.parent),
                       str(executable), str(inputs), str(directory/'worker.root'),
                       '0', '1', '0', '1', '2', '0', 'pPb', '0', '0']
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('JETANALYSIS_COMPLETED', result.stderr)
            self.assertEqual(subprocess.run(command, capture_output=True).returncode, 73)


if __name__ == '__main__':
    unittest.main()
