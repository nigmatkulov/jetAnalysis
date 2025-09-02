#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import shutil
import subprocess

def main():
    # --- Parse command line arguments ---
    parser = argparse.ArgumentParser(description="Resubmit failed condor jobs.")
    parser.add_argument("-log", required=True, help="Path to log directory")
    parser.add_argument("-sub", required=True, help="Path to submission directory")
    args = parser.parse_args()

    # --- Paths ---
    processing_dir_path = os.path.expanduser("~/soft/jetAnalysis/processing")
    log_dir = os.path.abspath(args.log)
    sub_dir = os.path.abspath(args.sub)

    # --- Check that directories exist ---
    for d in [processing_dir_path, log_dir, sub_dir]:
        if not os.path.isdir(d):
            print(f"ERROR: Directory does not exist: {d}")
            sys.exit(1)

    # --- Search for error patterns ---
    error_patterns = ["Input/output error", "scramv1: command not found"]
    unique_filenames = set()

    for err_file in glob.glob(os.path.join(log_dir, "*.err")):
        try:
            with open(err_file, "r", errors="ignore") as f:
                content = f.read()
                if any(pat in content for pat in error_patterns):
                    # Strip path and extension
                    base = os.path.basename(err_file)
                    name_no_ext = os.path.splitext(base)[0]
                    unique_filenames.add(name_no_ext)
        except Exception as e:
            print(f"WARNING: Could not read {err_file}: {e}")

    if not unique_filenames:
        print("No matching error jobs found.")
        return

    # --- Process each unique filename ---
    for fname in sorted(unique_filenames):
        sub_filename = f"pPb8160_{fname}.sub"
        src_sub_path = os.path.join(sub_dir, sub_filename)
        dst_sub_path = os.path.join(processing_dir_path, sub_filename)

        if not os.path.isfile(src_sub_path):
            print(f"WARNING: Submission file not found: {src_sub_path}")
            continue

        # Copy to processing dir
        shutil.copy2(src_sub_path, dst_sub_path)
        print(f"Copied: {src_sub_path} -> {dst_sub_path}")

        # Submit the job
        try:
            subprocess.run(["condor_submit", dst_sub_path], check=True)
            print(f"Submitted: {sub_filename}")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: Failed to submit {sub_filename}: {e}")

if __name__ == "__main__":
    main()

