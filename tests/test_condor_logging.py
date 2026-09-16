import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSING = PROJECT_ROOT / "processing"


class CondorLoggingTest(unittest.TestCase):
    def _common_files(self, directory: Path):
        proxy = directory / "proxy"
        proxy.write_text("test proxy\n")
        inputs = directory / "inputs.txt"
        inputs.write_text("/store/test/input.root\n")
        return proxy, inputs

    def test_process_submit_uses_bounded_log_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            proxy, inputs = self._common_files(directory)
            executable = directory / "processForestSimple"
            executable.write_text("#!/usr/bin/env bash\nexit 0\n")
            executable.chmod(0o755)
            work = directory / "condor"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(PROCESSING / "submit_process_forest_condor.py"),
                    str(inputs), str(directory / "output"),
                    "--executable", str(executable),
                    "--proxy", str(proxy), "--work-dir", str(work),
                    "--mc-type", "0", "--beam", "Pbgoing",
                    "--trigger-id", "1", "--dry-run",
                ],
                check=False, capture_output=True, text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            submit = next(work.rglob("*.sub")).read_text()
            self.assertIn("output = /dev/null", submit)
            self.assertNotIn("$(job_name).out", submit)
            self.assertIn("error = ", submit)
            self.assertIn("$(job_name).err", submit)
            self.assertIn("$(job_name).log", submit)
            self.assertIn("stream_output = False", submit)
            self.assertIn("stream_error = False", submit)
            self.assertIn("notification = Never", submit)

    def test_clone_submit_uses_bounded_log_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            proxy, inputs = self._common_files(directory)
            macro = directory / "cloneForest.C"
            macro.write_text("void cloneForest() {}\n")
            work = directory / "condor"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(PROCESSING / "submit_clone_forest_condor.py"),
                    str(inputs), str(directory / "output"),
                    "--macro", str(macro), "--proxy", str(proxy),
                    "--work-dir", str(work), "--dry-run",
                ],
                check=False, capture_output=True, text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            submit = next(work.rglob("*.sub")).read_text()
            self.assertIn("output = /dev/null", submit)
            self.assertNotIn("$(job_name).out", submit)
            self.assertIn("$(job_name).err", submit)
            self.assertIn("$(job_name).log", submit)

    def test_worker_caps_failure_output_and_preserves_read_signature(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            executable = directory / "fake-analysis"
            executable.write_text(
                "#!/usr/bin/env bash\n"
                "printf 'x%.0s' {1..100000}\n"
                "echo\n"
                "echo 'SysError in <TFile::ReadBuffer>: Input/output error'\n"
                "exit 0\n"
            )
            executable.chmod(0o755)
            input_list = directory / "input.list"
            input_list.write_text("input.root\n")
            completed = subprocess.run(
                [
                    str(PROCESSING / "run_process_forest_condor.sh"),
                    str(PROJECT_ROOT), str(executable), str(input_list),
                    str(directory / "output.root"), "0", "1", "0", "1",
                    "2", "0",
                ],
                check=False, capture_output=True,
            )
            self.assertEqual(completed.returncode, 74)
            self.assertEqual(completed.stdout, b"")
            self.assertLess(len(completed.stderr), 70_000)
            self.assertIn(b"SysError in <TFile::ReadBuffer>", completed.stderr)


if __name__ == "__main__":
    unittest.main()
