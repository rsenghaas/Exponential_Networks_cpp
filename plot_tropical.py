"""
plot_tropical.py
Reads tropical_segments.csv produced by tropical_gen and plots the diagram.

Usage:
    python plot_tropical.py                          # default input/output
    python plot_tropical.py --input path/to/csv      # custom CSV
    python plot_tropical.py --output path/to/out.pdf # custom output
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import numpy as np


def load_segments(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ("x0", "y0", "x1", "y1"):
        if col not in df.columns:
            raise ValueError(f"Expected column '{col}' not found in {path}")
    return df


def plot_tropical(df: pd.DataFrame, output: str) -> None:
    # Build a LineCollection for efficient rendering
    segments = [
        [(row.x0, row.y0), (row.x1, row.y1)]
        for _, row in df.iterrows()
    ]
    lc = mc.LineCollection(segments, colors="steelblue", linewidths=1.5)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.add_collection(lc)

    # Auto-fit axes with a small margin
    all_x = np.concatenate([df.x0.values, df.x1.values])
    all_y = np.concatenate([df.y0.values, df.y1.values])
    margin_x = 0.1 * (all_x.max() - all_x.min() or 1)
    margin_y = 0.1 * (all_y.max() - all_y.min() or 1)
    ax.set_xlim(all_x.min() - margin_x, all_x.max() + margin_x)
    ax.set_ylim(all_y.min() - margin_y, all_y.max() + margin_y)

    
    plt.axis([-10, 10, -10,10])
    # ax.set_xlabel("x")
    # ax.set_ylabel("y")
    # ax.set_title("Tropical diagram")
    # ax.set_aspect("equal")
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(output, dpi=1000)
    print(f"Saved to {output}")
    # plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot tropical diagram from CSV.")
    parser.add_argument("--input",  default="tropical/tropical_segments.csv",
                        help="Path to input CSV (default: tropical/tropical_segments.csv)")
    parser.add_argument("--output", default="graphics/toric/toric.pdf",
                        help="Path to output file (default: graphics/toric/toric.pdf)")
    args = parser.parse_args()

    df = load_segments(args.input)
    print(f"Loaded {len(df)} segments from {args.input}")
    plot_tropical(df, args.output)


if __name__ == "__main__":
    main()


"""
This code is great. The next issue is that there will be lots and lots and lots of line segments which in the end will give just a few effective straigt lines. Can you try to improve the code such that the effective straight lines are extracted, i.e. before writing to the file could you check if the segment connects/overlaps (up to error) withan existing segment with the same slope. Instead of just appending merge them together to a single segment. Maybe it's better to first collect the line segments and than write the few line segments in the end when everything is reduced down to the essential lines.
"""
