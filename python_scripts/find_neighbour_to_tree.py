"""Given a distance matrix computed with no SNP cutoff, produce a neighbour joining tree."""
import argparse
import time
from pathlib import Path
from collections import defaultdict

# These are optional dependencies to reduce the footprint of the library where not required.
try:
    import numpy as np
    from biotite.sequence.phylo import neighbor_joining
except ImportError:
    print("Please install FN6 with `pip install fn6[tree]` to use this script.")
    exit(1)


class Distance(object):
    def __init__(self, run1: str, run2: str, dist: int):
        self.run_id1 = run1
        self.run_id2 = run2
        self.dist = int(dist)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True, help="Path to the input file containing distances between runs")
    parser.add_argument("--output", type=Path, required=True, help="Path to the output file where the Newick tree will be saved")
    args = parser.parse_args()

    with open(args.input) as f:
        distances = [Distance(*[x for x in line.split(" ")]) for line in f]

    runs = sorted(
        list(
            set([d.run_id1 for d in distances]).union(
                set([d.run_id2 for d in distances])
            )
        )
    )
    # Track distances between samples and related samples symmetrically
    dists = {}
    run_clusters = defaultdict(set)
    for dist in distances:
        if dist.run_id1 != dist.run_id2:
            dists[(dist.run_id1, dist.run_id2)] = dist.dist
            dists[(dist.run_id2, dist.run_id1)] = dist.dist
            run_clusters[dist.run_id1].add(dist.run_id2)
            run_clusters[dist.run_id2].add(dist.run_id1)
    start = time.time()
    if len(runs) > 4:
        distances = []
        for run1 in runs:
            row = []
            for run2 in runs:
                if run1 == run2:
                    row.append(0)
                else:
                    d = dists.get((run1, run2))
                    if d is not None:
                        row.append(d)
                    else:
                        # Arbitrarily high distance for disjoint items
                        # Not that there should be any...
                        row.append(999999)
            distances.append(row)

    distances = np.array(distances)
    tree = neighbor_joining(distances)
    print(f"Built tree in {time.time() - start:.2f} seconds")
    with open(args.output, "w") as f:
        f.write(tree.to_newick(include_distance=True, labels=runs))


if __name__ == "__main__":
    main()
