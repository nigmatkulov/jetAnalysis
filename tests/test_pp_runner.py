"""Local runner checks without ROOT or production outputs."""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pp_analyzeMc_cp5tune as runner


class PpRunnerTests(unittest.TestCase):
    def test_dry_run_discovers_only_root_files_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            inputs = base / "forest"
            inputs.mkdir()
            for name in ("HiForest_2.root", "HiForest_1.root", "notes.txt"):
                (inputs / name).touch()
            output = base / "output"
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                status = runner.main(["--input-dir", str(inputs), "--output-dir", str(output),
                                      "--workers", "2", "--dry-run"])
            self.assertEqual(status, 0)
            text = stream.getvalue()
            self.assertIn("Found 2 files", text)
            self.assertIn(" 2 1 0 0 2 0 pp", text)
            self.assertNotIn("notes.txt", text)
            self.assertFalse(output.exists())

    def test_failure_is_reported_and_remaining_files_run_in_parallel(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            inputs = base / "forest"
            inputs.mkdir()
            for name in ("HiForest_1.root", "HiForest_2.root"):
                (inputs / name).touch()
            def task(command):
                return Path(command[2]), 1 if "_1.root" in command[1] else 0
            with patch.object(runner, "run_file", side_effect=task) as run, \
                    contextlib.redirect_stdout(io.StringIO()):
                status = runner.main(["--input-dir", str(inputs), "--output-dir", str(base / "out"),
                                      "--executable", "/usr/bin/true", "--workers", "2"])
            self.assertEqual(status, 1)
            self.assertEqual(run.call_count, 2)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            inputs = base / "forest"
            inputs.mkdir()
            (inputs / "HiForest_1.root").touch()
            existing = base / "HiForest_1_hist.root"
            existing.write_text("preserve")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()), \
                    self.assertRaises(SystemExit) as caught:
                runner.main(["--input-dir", str(inputs), "--output-dir", str(base),
                             "--executable", "/usr/bin/true"])
            self.assertEqual(caught.exception.code, 2)
            self.assertEqual(existing.read_text(), "preserve")


if __name__ == "__main__":
    unittest.main()
