"""
Match paths across CSV files where paths are reverse-traversals of each other
with y1 and y2 swapped.

CSV format per row: x (complex), log_y1 (complex), log_y2 (complex), <ignored>

Two paths A and B match if (after reversing A):
  - x_A[i] matches x_B[j] within X_TOL (hard distance gate)
  - exp(log_y1_A[i]) matches exp(log_y2_B[j]) within Y_TOL  (y swapped)
  - exp(log_y2_A[i]) matches exp(log_y1_B[j]) within Y_TOL
  - The pairing i->j over all matched points is order-preserving
    (corr(i, j) > 0.5), confirming true reversal not alignment
  - Sufficient arc-length of the overlap region is covered (>= ALIGN_FRACTION)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import KDTree
import multiprocessing as mp
from itertools import combinations
import sys
import os

# ---------------------------------------------------------------------------
# Hardcoded paths
# ---------------------------------------------------------------------------
DATA_DIR     = "data/path_data"
FILE_PATTERN = "path_data_*"

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
X_TOL          = 1e-1  # absolute distance gate for x matching
Y_TOL          = 1e-2  # absolute tolerance for exp(log_y) values
Y_LOG_TOL      = 1e-1  # absolute tolerance for |(log_y2 - log_y1) - (log_y1' - log_y2')|
ALIGN_FRACTION = 0.05   # min fraction of overlap arc-length that must be covered

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def parse_complex(s: str) -> complex:
    return complex(s.replace(" ", ""))

def load_path(filepath: str) -> dict | None:
    """Load a path file. Returns dict with arrays, or None on error."""
    try:
        rows = []
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = [p.strip() for p in line.split(",")]
                if len(parts) < 3:
                    continue
                x      = parse_complex(parts[0])
                log_y1 = parse_complex(parts[1])
                log_y2 = parse_complex(parts[2])
                # parts[3] (last column) is ignored
                rows.append((x, log_y1, log_y2))
        if not rows:
            return None
        x      = np.array([r[0] for r in rows], dtype=complex)
        log_y1 = np.array([r[1] for r in rows], dtype=complex)
        log_y2 = np.array([r[2] for r in rows], dtype=complex)
        y1 = np.exp(log_y1)
        y2 = np.exp(log_y2)
        seg_lengths = np.abs(np.diff(x))
        arc = np.concatenate([[0.0], np.cumsum(seg_lengths)])
        total_arc = arc[-1]
        return dict(
            path=filepath,
            x=x,
            y1=y1,
            y2=y2,
            log_y1=log_y1,
            log_y2=log_y2,
            arc=arc,
            total_arc=total_arc,
            n=len(x),
        )
    except Exception as e:
        print(f"Error loading {filepath}: {e}", file=sys.stderr)
        return None

# ---------------------------------------------------------------------------
# Arc-length coverage
# ---------------------------------------------------------------------------

def _arc_coverage(matched_i: np.ndarray, seg: np.ndarray, n: int) -> tuple[float, float]:
    """
    Compute arc covered by matched points and the overlap-region denominator.
    Each matched point gets credit for its Voronoi cell (half-segment each side).
    Denominator = arc from first to last matched point.
    Returns (covered, fraction).
    """
    if len(matched_i) == 0:
        return 0.0, 0.0
    covered = 0.0
    for i in matched_i:
        left  = seg[i - 1] / 2 if i > 0     else 0.0
        right = seg[i]     / 2 if i < n - 1 else 0.0
        covered += left + right
    i_first, i_last = matched_i[0], matched_i[-1]
    overlap_arc = float(np.sum(seg[i_first:i_last])) if i_last > i_first else 0.0
    overlap_arc = max(overlap_arc, covered)
    fraction = covered / overlap_arc if overlap_arc > 0 else 0.0
    return covered, fraction

# ---------------------------------------------------------------------------
# Direction check
# ---------------------------------------------------------------------------

def _is_reversed(matched_i: np.ndarray, matched_j: np.ndarray) -> tuple[bool, float]:
    """
    Check that the pairing runs in the correct direction for a true reversal.
    matched_i indexes into reversed-A (always ascending from np.where).
    For a true reversal, matched_j must also ascend as matched_i ascends.
    An aligned pair gives matched_j descending -> negative correlation.
    Returns (direction_ok, correlation).
    Division by zero (constant i or j) is handled explicitly.
    """
    if len(matched_i) < 2:
        return True, 1.0
    std_i = matched_i.std()
    std_j = matched_j.std()
    if std_i == 0 or std_j == 0:
        # All matches map to a single point — cannot determine direction
        return False, float('nan')
    corr = np.cov(matched_i.astype(float), matched_j.astype(float))[0, 1] / (std_i * std_j)
    return corr > 0.5, float(corr)

# ---------------------------------------------------------------------------
# Core matching logic
# ---------------------------------------------------------------------------

def align_and_check(A: dict, B: dict) -> tuple[bool, float]:
    """
    Check if path A reversed matches path B with y1<->y2 swap.
    Returns (matched, arc_fraction).

    Match conditions (all must hold simultaneously):
      |exp(log_y1_A) - exp(log_y2_B)| < Y_TOL
      |exp(log_y2_A) - exp(log_y1_B)| < Y_TOL
      |(log_y2_A - log_y1_A) - (log_y1_B - log_y2_B)| < Y_LOG_TOL
    """
    xA    = A["x"][::-1]
    y1A   = A["y1"][::-1]
    y2A   = A["y2"][::-1]
    ly1A  = A["log_y1"][::-1]
    ly2A  = A["log_y2"][::-1]
    xB    = B["x"]
    y1B   = B["y1"]
    y2B   = B["y2"]
    ly1B  = B["log_y1"]
    ly2B  = B["log_y2"]
    nA, nB = len(xA), len(xB)

    pts_B = np.column_stack([xB.real, xB.imag])
    tree  = KDTree(pts_B)
    pts_A = np.column_stack([xA.real, xA.imag])
    dists, idxs = tree.query(pts_A, k=1, distance_upper_bound=X_TOL)

    # Collect matched pairs: x within X_TOL, both exp(y) within Y_TOL,
    # and log_y2 - log_y1 matches (log_y1_B - log_y2_B) within Y_LOG_TOL
    matched_i = []
    matched_j = []
    for i in np.where(dists < X_TOL)[0]:
        j = idxs[i]
        if j < nB and abs(y1A[i] - y2B[j]) < Y_TOL and abs(y2A[i] - y1B[j]) < Y_TOL                 and abs((ly2A[i] - ly1A[i]) - (ly1B[j] - ly2B[j])) < Y_LOG_TOL:
            matched_i.append(i)
            matched_j.append(j)

    if not matched_i:
        return False, 0.0

    matched_i = np.array(matched_i)
    matched_j = np.array(matched_j)

    # Direction check over all matched points together
    direction_ok, _ = _is_reversed(matched_i, matched_j)
    if not direction_ok:
        return False, 0.0

    seg_A = np.abs(np.diff(xA))
    covered, fraction = _arc_coverage(matched_i, seg_A, nA)
    return fraction >= ALIGN_FRACTION, fraction

# ---------------------------------------------------------------------------
# Pairwise check worker
# ---------------------------------------------------------------------------

def check_pair(args):
    i, j, paths_data = args
    A = paths_data[i]
    B = paths_data[j]
    matched, frac = align_and_check(A, B)
    if matched:
        return (A["path"], B["path"], frac)
    return None

# ---------------------------------------------------------------------------
# Full scan
# ---------------------------------------------------------------------------

def find_matching_paths(n_workers: int | None = None) -> list[tuple]:
    files = sorted(Path(DATA_DIR).glob(FILE_PATTERN))
    print(f"Found {len(files)} files in {DATA_DIR}/.")

    print("Loading files...")
    paths_data = []
    for f in files:
        d = load_path(str(f))
        if d is not None:
            paths_data.append(d)
    print(f"Loaded {len(paths_data)} valid paths.")

    pairs = list(combinations(range(len(paths_data)), 2))
    print(f"Checking {len(pairs)} pairs...")

    if n_workers is None:
        n_workers = max(1, os.cpu_count() - 1)

    worker_args = [(i, j, paths_data) for i, j in pairs]

    matches = []
    with mp.Pool(n_workers) as pool:
        for result in pool.imap_unordered(check_pair, worker_args, chunksize=64):
            if result is not None:
                matches.append(result)
                print(f"  MATCH: {Path(result[0]).name} <-> {Path(result[1]).name}  arc={result[2]:.3f}")

    return matches

# ---------------------------------------------------------------------------
# Inspect a specific pair
# ---------------------------------------------------------------------------

def inspect_pair(k1: int, k2: int) -> None:
    """
    Load path_data_{k1}.csv and path_data_{k2}.csv and print a detailed match report.
    """
    file_a = str(Path(DATA_DIR) / f"path_data_{k1}.csv")
    file_b = str(Path(DATA_DIR) / f"path_data_{k2}.csv")

    A = load_path(file_a)
    B = load_path(file_b)
    if A is None:
        print(f"Could not load {file_a}", file=sys.stderr); return
    if B is None:
        print(f"Could not load {file_b}", file=sys.stderr); return

    xA   = A["x"][::-1]
    y1A  = A["y1"][::-1]
    y2A  = A["y2"][::-1]
    ly1A = A["log_y1"][::-1]
    ly2A = A["log_y2"][::-1]
    xB   = B["x"]
    y1B  = B["y1"]
    y2B  = B["y2"]
    ly1B = B["log_y1"]
    ly2B = B["log_y2"]
    nA, nB = len(xA), len(xB)

    pts_B = np.column_stack([xB.real, xB.imag])
    tree  = KDTree(pts_B)
    pts_A = np.column_stack([xA.real, xA.imag])
    # Wide net for inspection so near-misses are visible too
    dists, idxs = tree.query(pts_A, k=1, distance_upper_bound=X_TOL * 1000)

    seg = np.abs(np.diff(xA))

    print(f"\n{'='*75}")
    print(f"  Inspecting: path_data_{k1}.csv (reversed)  <->  path_data_{k2}.csv")
    print(f"  Points: A={nA}, B={nB}  |  Arc: A={A['total_arc']:.6f}, B={B['total_arc']:.6f}")
    print(f"  X_TOL={X_TOL:.2e}  Y_TOL={Y_TOL:.2e}  Y_LOG_TOL={Y_LOG_TOL:.2e}  ALIGN_FRACTION={ALIGN_FRACTION:.2f}")
    print(f"{'='*75}")
    print(f"  {'i_A':>5}  {'j_B':>5}  {'x_dist':>10}  {'|y1A-y2B|':>12}  {'|y2A-y1B|':>12}  {'|dly_A-dly_B|':>14}  {'ok':>5}  {'arc_contrib':>11}")
    print(f"  {'-'*90}")

    matched_i = []
    matched_j = []
    for i in range(nA):
        if dists[i] >= X_TOL * 1000:
            continue
        j = idxs[i]
        if j >= nB:
            continue
        x_dist  = dists[i]
        dy1     = abs(y1A[i] - y2B[j])
        dy2     = abs(y2A[i] - y1B[j])
        dlog    = abs((ly2A[i] - ly1A[i]) - (ly1B[j] - ly2B[j]))
        x_ok    = x_dist < X_TOL
        y_ok    = dy1 < Y_TOL and dy2 < Y_TOL and dlog < Y_LOG_TOL
        left    = seg[i - 1] / 2 if i > 0      else 0.0
        right   = seg[i]     / 2 if i < nA - 1 else 0.0
        contrib = left + right if (x_ok and y_ok) else 0.0
        if x_ok and y_ok:
            matched_i.append(i)
            matched_j.append(j)
        marker = "✓" if (x_ok and y_ok) else ("~" if x_ok else " ")
        print(f"  {marker} {i:>5}  {j:>5}  {x_dist:>10.2e}  {dy1:>12.2e}  {dy2:>12.2e}  {dlog:>14.2e}  {str(x_ok and y_ok):>5}  {contrib:>11.6f}")

    print(f"  {'-'*80}")

    matched_i = np.array(matched_i) if matched_i else np.array([], dtype=int)
    matched_j = np.array(matched_j) if matched_j else np.array([], dtype=int)

    covered, frac_arc = _arc_coverage(matched_i, seg, nA)
    direction_ok, corr = _is_reversed(matched_i, matched_j)

    matched = direction_ok and frac_arc >= ALIGN_FRACTION

    print(f"  Matched points   : {len(matched_i)} / {nA}")
    print(f"  Direction corr   : {corr:.3f}  [must be > +0.5]  {'ok' if direction_ok else 'FAIL — paths not reversed!'}")
    print(f"  Arc coverage     : {covered:.6f}  ({frac_arc*100:.1f}% of overlap)  [threshold: {ALIGN_FRACTION*100:.0f}%]  {'ok' if frac_arc >= ALIGN_FRACTION else 'FAIL'}")
    print(f"  Result           : {'MATCH ✓' if matched else 'NO MATCH ✗'}")
    print(f"{'='*75}\n")

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Find reverse-matching paths in data/path_data/.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_scan = subparsers.add_parser("scan", help="Check all pairs of paths.")
    p_scan.add_argument("--workers", type=int, default=None)
    p_scan.add_argument("--output",  default=None, help="Save matches to this CSV file")

    p_inspect = subparsers.add_parser("inspect", help="Inspect a specific pair.")
    p_inspect.add_argument("k1", type=int, help="Index of first path  (path_data_{k1}.csv)")
    p_inspect.add_argument("k2", type=int, help="Index of second path (path_data_{k2}.csv)")

    for p in (p_scan, p_inspect):
        p.add_argument("--x-tol",          type=float, default=X_TOL)
        p.add_argument("--y-tol",          type=float, default=Y_TOL)
        p.add_argument("--y-log-tol",      type=float, default=Y_LOG_TOL)
        p.add_argument("--match-fraction", type=float, default=ALIGN_FRACTION)

    args = parser.parse_args()
    X_TOL          = args.x_tol
    Y_TOL          = args.y_tol
    Y_LOG_TOL      = args.y_log_tol
    ALIGN_FRACTION = args.match_fraction

    if args.command == "inspect":
        inspect_pair(args.k1, args.k2)
    elif args.command == "scan":
        matches = find_matching_paths(args.workers)
        print(f"\nFound {len(matches)} matching pairs.")
        if args.output:
            df = pd.DataFrame(matches, columns=["file_a", "file_b", "arc_fraction"])
            df.to_csv(args.output, index=False)
            print(f"Saved to {args.output}")
        else:
            for m in matches:
                print(f"{m[0]}  <->  {m[1]}  (arc={m[2]:.3f})")
