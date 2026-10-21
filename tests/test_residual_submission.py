import contextlib
import io
import fnmatch
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
PROCESSING = ROOT / 'processing'


def module(name):
    spec = importlib.util.spec_from_file_location(name, PROCESSING / f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


class ResidualSubmissionTest(unittest.TestCase):
    def test_submit_and_worker_modes(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            inputs = directory / 'inputs.list'
            inputs.write_text('/eos/cms/store/input.root\n')
            proxy = directory / 'proxy'
            proxy.write_text('test proxy\n')
            capture = directory / 'argv'
            executable = directory / 'analysis'
            executable.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$@" > ' + str(capture) + '\n')
            executable.chmod(0o755)
            for mode, suffix in [(2, ''), (1, '_oldResiduals'), (0, '_noResiduals')]:
                work = directory / f'work{mode}'
                command = [sys.executable, str(PROCESSING / 'submit_process_forest_condor.py'),
                           str(inputs), str(directory / 'out'), '--trigger-id', '1',
                           '--proxy', str(proxy), '--executable', str(executable),
                           '--work-dir', str(work), '--dry-run']
                if mode != 2:
                    command += ['--residual-mode', str(mode)]
                result = subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                text = next(work.rglob('*.sub')).read_text()
                self.assertEqual(next(work.rglob('*.list')).read_text(),
                                 'root://eoscms.cern.ch//eos/cms/store/input.root\n')
                self.assertNotIn('request_disk', text)
                self.assertNotIn('runtime.tar.gz', text)
                name = f'Jet60_Pbgoing_data_jetId{suffix}_job0001.root'
                self.assertIn(name, text)
                arguments = next(line for line in text.splitlines() if line.startswith('arguments = '))
                self.assertTrue(arguments.endswith(f'pPb 0 0 {mode}'), arguments)
                # Execute the generated worker arguments, replacing Condor macros.
                tokens = arguments.removeprefix('arguments = ').split()
                tokens = [str(inputs) if t == '$(input_list)' else
                          str(directory / name) if t == '$(output_file)' else t for t in tokens]
                result = subprocess.run([str(PROCESSING / 'run_process_forest_condor.sh'), *tokens],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(capture.read_text().splitlines()[-3:], ['0', '0', str(mode)])
            # Previously generated submit files retain their implicit mode 2.
            result = subprocess.run([str(PROCESSING / 'run_process_forest_condor.sh'), *tokens[:-1]],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(capture.read_text().splitlines()[-1], '2')

    def test_xrootd_input_mapping(self):
        submit = module('submit_process_forest_condor')
        cases = {
            '/eos/cms/store/a.root': 'root://eoscms.cern.ch//eos/cms/store/a.root',
            '/eos/user/g/user/a.root': 'root://eosuser.cern.ch//eos/user/g/user/a.root',
            '/eos/home-g/user/a.root': 'root://eosuser.cern.ch//eos/user/g/user/a.root',
            'root://host//path/a.root?x=1': 'root://host//path/a.root?x=1',
            '/tmp/a.root': '/tmp/a.root',
            'relative.root': 'relative.root',
        }
        for source, expected in cases.items():
            self.assertEqual(submit.xrootd_input(source), expected)
        with self.assertRaisesRegex(SystemExit, 'Unsupported EOS input namespace'):
            submit.xrootd_input('/eos/unknown/a.root')
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'input.list'
            content = '# comment\n\n/eos/cms/store/a.root\n/tmp/b.root\n/eos/cms/store/a.root\n'
            path.write_text(content)
            self.assertEqual(submit.read_inputs(path), [cases['/eos/cms/store/a.root'],
                                                       '/tmp/b.root', cases['/eos/cms/store/a.root']])
            self.assertEqual(path.read_text(), content)

    def test_mc_names_do_not_change(self):
        submit = module('submit_process_forest_condor')
        with patch.object(sys, 'argv', ['submit', 'input.list', 'out', '--mc-type', '1', '--pt-hat', '540']):
            args = submit.parse_args()
        for mode in [0, 1, 2]:
            args.residual_mode = mode
            self.assertEqual(submit.make_name(args, None, 540),
                             'NoTrigger_Pbgoing_embedding_ptHat540_jetId')

    def test_all_data_wrapper_forwards_mode(self):
        submit = module('submit_all_data_condor')
        for mode in [0, 1, 2]:
            with patch.object(sys, 'argv', ['submit', '--residual-mode', str(mode)]):
                args = submit.parse_args()
            for beam in ['pgoing', 'Pbgoing']:
                for trigger in range(4):
                    command = submit.submission_command(args, beam, trigger)
                    self.assertEqual(command[command.index('--residual-mode') + 1], str(mode))

    def test_merge_patterns_separate_variants(self):
        merge = module('merge_data_outputs')
        for trigger in ['MB', 'Jet60']:
            for mode, suffix in [(2, ''), (1, '_oldResiduals'), (0, '_noResiduals')]:
                pattern = merge.input_pattern(trigger, 'Pbgoing', 'jetId', 'dz1p0', mode)
                prefix = 'MB_PD1' if trigger == 'MB' else trigger
                for other in ['', '_oldResiduals', '_noResiduals']:
                    name = f'{prefix}_Pbgoing_data_jetId{other}_job0001.root'
                    self.assertEqual(fnmatch.fnmatchcase(name, pattern), other == suffix)
                with tempfile.TemporaryDirectory() as temporary:
                    directory = Path(temporary)
                    for other in ['', '_oldResiduals', '_noResiduals']:
                        (directory / f'{prefix}_Pbgoing_data_jetId{other}_job0001.root').touch()
                    result = subprocess.run([sys.executable, str(PROCESSING / 'merge_data_outputs.py'),
                                             str(directory), '--trigger', trigger, '--beam', 'Pbgoing',
                                             '--residual-mode', str(mode), '--dry-run'],
                                            capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn('Selected 1 source file(s)', result.stdout)
                    self.assertIn(f'{trigger.lower()}_Pbgoing_ak4_jetId{suffix}.root', result.stdout)

    def test_invalid_submission_modes(self):
        submit = module('submit_process_forest_condor')
        for value in ['-1', '3', 'bad']:
            with patch.object(sys, 'argv', ['submit', 'input.list', 'out', '--residual-mode', value]), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    submit.parse_args()
                self.assertEqual(error.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
