#!/usr/bin/env python3
"""
Collect screening results from all chunk TSV files into a single results file.
Then run analysis to generate top 500 reports.

Usage:
  python3 collect_results.py [--db npatlas|coconut|both]

Looks for chunk files in ~/final_project/screening_results/
Outputs combined TSV and runs analyze_screening_results.py automatically.
"""
import argparse
import glob
import os
import subprocess
import sys


def collect_chunks(db_name, results_dir):
    """Merge all chunk TSV files for a database into one file."""
    pattern = os.path.join(results_dir, f"{db_name}_chunk_*.tsv")
    chunk_files = sorted(glob.glob(pattern), key=lambda x: int(x.split("_chunk_")[1].split(".")[0]))

    if not chunk_files:
        print(f"  No chunk files found for {db_name}")
        return None

    output_path = os.path.join(results_dir, f"{db_name}_all_results.tsv")
    total = 0
    ok_count = 0

    with open(output_path, "w") as fout:
        fout.write("name\tsmiles\tscore\tstatus\n")

        for chunk_file in chunk_files:
            with open(chunk_file) as fin:
                header = fin.readline()
                for line in fin:
                    fout.write(line)
                    total += 1
                    if "\tok\n" in line or "\tok\t" in line:
                        ok_count += 1

    print(f"  {db_name}: {len(chunk_files)} chunks, {total} compounds, {ok_count} successfully docked")
    print(f"  Combined results: {output_path}")
    return output_path


def check_progress(db_name, results_dir, total_compounds):
    """Check how many compounds have been processed so far."""
    pattern = os.path.join(results_dir, f"{db_name}_chunk_*.tsv")
    chunk_files = glob.glob(pattern)
    processed = 0
    for f in chunk_files:
        with open(f) as fin:
            processed += sum(1 for line in fin) - 1
    pct = (processed / total_compounds * 100) if total_compounds > 0 else 0
    return processed, pct


def main():
    parser = argparse.ArgumentParser(description="Collect and analyze screening results")
    parser.add_argument("--db", default="both", choices=["npatlas", "coconut", "both"])
    parser.add_argument("--progress", action="store_true", help="Just show progress, don't merge")
    args = parser.parse_args()

    project_dir = os.path.expanduser("~/final_project")
    results_dir = os.path.join(project_dir, "screening_results")
    db_dir = os.path.join(project_dir, "np_databases")
    analyze_script = os.path.join(project_dir, "analyze_screening_results.py")

    databases = []
    if args.db in ("npatlas", "both"):
        databases.append("npatlas")
    if args.db in ("coconut", "both"):
        databases.append("coconut")

    if args.progress:
        print("=" * 50)
        print("Screening Progress")
        print("=" * 50)
        for db_name in databases:
            smi_path = os.path.join(db_dir, f"{db_name}.smi")
            total = sum(1 for _ in open(smi_path)) if os.path.exists(smi_path) else 0
            processed, pct = check_progress(db_name, results_dir, total)
            print(f"  {db_name}: {processed}/{total} ({pct:.1f}%)")

            log_pattern = os.path.join(results_dir, "logs", f"{db_name}_*.err")
            err_files = glob.glob(log_pattern)
            errors = 0
            for ef in err_files:
                if os.path.getsize(ef) > 0:
                    errors += 1
            if errors:
                print(f"    {errors} tasks had errors (check logs/)")
        return

    print("=" * 50)
    print("Collecting Screening Results")
    print("=" * 50)

    result_files = []
    for db_name in databases:
        print(f"\n{db_name}:")
        result_path = collect_chunks(db_name, results_dir)
        if result_path:
            result_files.append(result_path)

    if not result_files:
        print("\nNo results found. Make sure screening jobs have completed.")
        print("Check progress: python3 collect_results.py --progress")
        return

    combined_path = os.path.join(results_dir, "all_databases_results.tsv")
    if len(result_files) > 1:
        print(f"\nMerging all databases into {combined_path}")
        with open(combined_path, "w") as fout:
            fout.write("name\tsmiles\tscore\tstatus\n")
            for rf in result_files:
                with open(rf) as fin:
                    fin.readline()
                    for line in fin:
                        fout.write(line)
        result_files = [combined_path]

    analysis_input = result_files[0]
    print(f"\nRunning analysis on: {analysis_input}")

    if os.path.exists(analyze_script):
        subprocess.run([
            sys.executable, analyze_script,
            analysis_input, "--top", "500"
        ])
    else:
        print(f"Analysis script not found at {analyze_script}")
        print(f"Download it: curl -sL https://raw.githubusercontent.com/anirudhatalmale6-alt/autogrow-filter-fix/main/analyze_screening_results.py -o {analyze_script}")


if __name__ == "__main__":
    main()
