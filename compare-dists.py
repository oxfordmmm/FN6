"""Helper for comparing 2 sets of distances for functional equality.
Saves can be produced in any order, both in terms of rows and samples
"""
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two distance files for functional equality")
    parser.add_argument("file1", help="First distance file")
    parser.add_argument("file2", help="Second distance file")
    args = parser.parse_args()
    with open(args.file1, "r") as f:
        dists1 = [line.strip().split(" ") for line in f]
        dists1 = {(sorted([s1, s2])[0], sorted([s1, s2])[1], dist) for s1, s2, dist in dists1}

    with open(args.file2, "r") as f:
        dists2 = [line.strip().split(" ") for line in f]
        dists2 = {(sorted([s1, s2])[0], sorted([s1, s2])[1], dist) for s1, s2, dist in dists2}
    
    if dists1 == dists2:
        print("The distance files are functionally equivalent.")
    else:
        print("The distance files differ.")
        only_in_file1 = dists1 - dists2
        only_in_file2 = dists2 - dists1
        in_both = dists1 & dists2

        if only_in_file1:
            print("\nDistances only in file 1:")
            for dist in only_in_file1:
                print(dist)

        if only_in_file2:
            print("\nDistances only in file 2:")
            for dist in only_in_file2:
                print(dist)
        
        print(f"\nDistances in both files: {len(in_both)}")
