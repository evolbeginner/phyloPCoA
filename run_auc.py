#!/usr/bin/env python3
"""Run phyloPCoA simulations for exponent 0 and a comparison exponent."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import math
import random
import shlex
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from tqdm import tqdm


METHODS = ("pcoa_2", "pcoa_3", "pcoa_1", "pcoa_4")
METHOD_LABELS = {
    "pcoa_2": "CLR distance (raw)",
    "pcoa_3": "CLR distance (phyloPCoA-adjusted)",
    "pcoa_1": "BC distance (raw)",
    "pcoa_4": "BC distance (phyloPCoA-adjusted)",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare exponent 0 with a selected exponent and calculate AUC from adonis p-values."
    )
    parser.add_argument("--exponent", type=float, default=1,
                        help="Comparison exponent (default: 1; exponent 0 is the baseline)")
    parser.add_argument("--script", type=Path,
                        default=Path.home() / "lab-tools/phyloPCoA/phyloPCoA.R",
                        help="Path to phyloPCoA.R")
    parser.add_argument("--outdir", type=Path, default=Path("auc_runs"),
                        help="Directory for individual simulation outputs and summary files")
    parser.add_argument("--runs", type=int, default=100,
                        help="Number of replicates per exponent (default: 100)")
    parser.add_argument("--seed", type=int,
                        help="Base seed for consecutive simulation seeds (default: random seeds)")
    parser.add_argument("--cpu", type=int, default=1,
                        help="Maximum number of simulations to run in parallel (default: 1)")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    parser.add_argument("--bd", default="1,1,0.01,10")
    parser.add_argument("-T", type=int, default=50, dest="tnum")
    parser.add_argument("-B", type=int, default=30, dest="bnum")
    parser.add_argument("--binary-rate", default=".010")
    parser.add_argument("--pagel_lam_sim",
                        help="Simulated Pagel's lambda (number or beta:alpha,beta)")
    parser.add_argument("--pagel_lam_mode",
                        choices=("global", "per_feature", "hierarchical", "beta", "none", "auto"),
                        help="Pagel's lambda estimation mode")
    return parser.parse_args()


def read_pvalues(path):
    found = {}
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"pcoa", "P_value"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Expected columns {sorted(required)} in {path}")
        for row in reader:
            method = row["pcoa"]
            if method in METHODS:
                try:
                    pvalue = float(row["P_value"])
                except (TypeError, ValueError):
                    continue
                if math.isfinite(pvalue):
                    found[method] = pvalue
    missing = set(METHODS) - found.keys()
    if missing:
        raise ValueError(f"Missing valid {', '.join(sorted(missing))} rows in {path}")
    return found


def auc(scores):
    """Mann-Whitney ROC AUC with average ranks for tied scores."""
    positives = [score for label, score in scores if label == 1]
    negatives = [score for label, score in scores if label == 0]
    if not positives or not negatives:
        raise ValueError("AUC requires observations from both exponent groups")
    ordered = sorted((score, label) for label, score in scores)
    rank_sum = 0.0
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][0] == ordered[i][0]:
            j += 1
        average_rank = ((i + 1) + j) / 2
        rank_sum += average_rank * sum(label for _, label in ordered[i:j])
        i = j
    return (rank_sum - len(positives) * (len(positives) + 1) / 2) / (len(positives) * len(negatives))


def roc_points(scores):
    """Return false-positive and true-positive rates, grouping tied scores."""
    positives = sum(label == 1 for label, _ in scores)
    negatives = sum(label == 0 for label, _ in scores)
    if not positives or not negatives:
        raise ValueError("ROC requires observations from both exponent groups")

    points = [(0.0, 0.0)]
    true_positives = false_positives = 0
    ordered = sorted(scores, key=lambda item: item[1], reverse=True)
    i = 0
    while i < len(ordered):
        j = i
        while j < len(ordered) and ordered[j][1] == ordered[i][1]:
            if ordered[j][0] == 1:
                true_positives += 1
            else:
                false_positives += 1
            j += 1
        points.append((false_positives / negatives, true_positives / positives))
        i = j
    return points


def write_roc_pdf(path, observations, methods):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 7))
    for method in methods:
        points = roc_points(observations[method])
        false_positive_rates, true_positive_rates = zip(*points)
        area = auc(observations[method])
        ax.plot(false_positive_rates, true_positive_rates,
                label=f"{METHOD_LABELS[method]} (AUC = {area:.3f})")
    ax.plot((0, 1), (0, 1), linestyle="--", color="gray", label="Chance")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="False positive rate",
           ylabel="True positive rate", title="ROC curves")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def archive_replicates(outdir, exponent):
    """Compress replicate directories after retaining the first run."""
    label = format(exponent, "g")
    run_dirs = sorted(
        path for path in outdir.glob(f"exp{label}_run*") if path.is_dir()
    )
    if len(run_dirs) < 2:
        return None

    folders_to_archive = run_dirs[10:]
    archive_path = outdir / f"exp{label}_runs.tar.xz"
    with tarfile.open(archive_path, "w:xz") as archive:
        for run_dir in folders_to_archive:
            archive.add(run_dir, arcname=run_dir.name)

    for run_dir in folders_to_archive:
        shutil.rmtree(run_dir)
    return archive_path


def main():
    args = parse_args()
    if args.exponent == 0:
        raise SystemExit("--exponent must differ from the baseline exponent 0")
    if args.runs < 1:
        raise SystemExit("--runs must be at least 1")
    if args.cpu < 1:
        raise SystemExit("--cpu must be at least 1")
    if args.seed is not None and not (0 <= args.seed <= 2**31 - 2 * args.runs):
        raise SystemExit("--seed must allow one valid R integer seed per simulation")
    if not args.script.is_file():
        raise SystemExit(f"R script not found: {args.script}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    observations = {method: [] for method in METHODS}
    score_path = args.outdir / "auc_scores.tsv"
    exponent_label = format(args.exponent, "g")

    def simulation_command(exponent, run, seed):
        exponent_text = format(exponent, "g")
        run_dir = args.outdir.resolve() / f"exp{exponent_text}_run{run:04d}"
        command = [
            args.rscript, str(args.script.resolve()), "--bd", args.bd,
            "-T", str(args.tnum), "-B", str(args.bnum),
            "--outdir", str(run_dir), "--force",
            "--exponent", exponent_text, "--sim",
            "--binary_rate", args.binary_rate,
        ]
        if args.pagel_lam_sim is not None:
            command.extend(("--pagel_lam_sim", args.pagel_lam_sim))
        if args.pagel_lam_mode is not None:
            command.extend(("--pagel_lam_mode", args.pagel_lam_mode))
        if seed is not None:
            command.extend(("--seed", str(seed)))
        return command, run_dir

    tasks = [(exponent, run) for exponent in (0, args.exponent)
             for run in range(1, args.runs + 1)]
    seeds = (random.SystemRandom().sample(range(1, 2**31), len(tasks))
             if args.seed is None else range(args.seed, args.seed + len(tasks)))
    commands = {
        task: simulation_command(*task, seed)
        for task, seed in zip(tasks, seeds)
    }
    command_path = args.outdir / "phyloPCoA_commands.txt"
    command_path.write_text(
        "\n".join(shlex.join(commands[task][0]) for task in tasks) + "\n"
    )
    print(f"Commands: {command_path}", flush=True)

    def run_simulation(exponent, run):
        command, run_dir = commands[(exponent, run)]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        table = run_dir / "adonis" / "by_trait_all_axes.tbl"
        return exponent, run, read_pvalues(table)

    results = {}
    with ThreadPoolExecutor(max_workers=args.cpu) as executor:
        futures = {executor.submit(run_simulation, exponent, run): (exponent, run)
                   for exponent, run in tasks}
        with tqdm(total=len(futures), desc="Simulations", unit="run") as progress:
            for future in as_completed(futures):
                exponent, run = futures[future]
                results[(exponent, run)] = future.result()
                progress.update(1)

    with score_path.open("w", newline="") as score_file:
        writer = csv.writer(score_file, delimiter="\t")
        writer.writerow(("exponent", "run", "method", "p_value",
                         f"score_for_exponent_{exponent_label}"))
        for exponent, run in tasks:
            _, _, pvalues = results[(exponent, run)]
            for method, pvalue in pvalues.items():
                score = 1.0 - pvalue
                observations[method].append((int(exponent == args.exponent), score))
                writer.writerow((exponent, run, method, pvalue, score))

    summary_path = args.outdir / "auc_summary.tsv"
    with summary_path.open("w", newline="") as summary_file:
        writer = csv.writer(summary_file, delimiter="\t")
        writer.writerow(("method", "n_per_exponent", "score", "auc"))
        for method in METHODS:
            auc_value = auc(observations[method])
            writer.writerow((method, args.runs, "1 - P_value", auc_value))
            print(f"{method}: AUC = {auc_value:.6f}")
    roc_outputs = (
        ("ROC_CLR.pdf", ("pcoa_3", "pcoa_2")),
        ("ROC_BC.pdf", ("pcoa_1", "pcoa_4")),
        ("ROC_all.pdf", METHODS),
    )
    for filename, methods in roc_outputs:
        roc_path = args.outdir / filename
        write_roc_pdf(roc_path, observations, methods)
        print(f"ROC:     {roc_path}")
    for exponent in (0, args.exponent):
        archive_path = archive_replicates(args.outdir, exponent)
        if archive_path:
            print(f"Archived: {archive_path}")
    print(f"Scores:  {score_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        sys.exit(f"Simulation failed (exit code {error.returncode}); see the run output above.")
