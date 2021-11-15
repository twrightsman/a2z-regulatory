#!/usr/bin/env python

import argparse
import sys
import warnings

import h5py
import numpy as np


def main(config: str, species: str) -> None:
    in_path = f"tmp/tf-modisco/{config}/{species.replace(' ', '_')}.hdf5"
    hdf5_results = h5py.File(in_path, "r")

    with open(f"tmp/tf-modisco/{config}/{species.replace(' ', '_')}.meme", "w") as meme_out:
        meme_out.write("MEME version 5\n\n")
        meme_out.write("ALPHABET= ACGT\n\n")
        meme_out.write("strands: + -\n\n")
        meme_out.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")

        metacluster_names = [name.decode("utf-8") for name in hdf5_results["metaclustering_results"]["all_metacluster_names"]]

        pattern_to_agg_score = {}
        for metacluster_name in metacluster_names:
            metacluster_grp = hdf5_results["metacluster_idx_to_submetacluster_results"][metacluster_name]
            all_pattern_names = [name.decode("utf-8") for name in metacluster_grp["seqlets_to_patterns_result"]["patterns"]["all_pattern_names"]]
            for pattern_name in all_pattern_names:
                pattern = metacluster_grp["seqlets_to_patterns_result"]["patterns"][pattern_name]
                pattern_scores = pattern["task0_contrib_scores"]["fwd"][:]
                pattern_to_agg_score[(metacluster_name, pattern_name)] = np.abs(pattern_scores).sum()

        # sort patterns by the sum of the absolute value of their contribution scores
        for metacluster_name, pattern_name in sorted(pattern_to_agg_score, key = pattern_to_agg_score.__getitem__):
            metacluster_grp = hdf5_results["metacluster_idx_to_submetacluster_results"][metacluster_name]
            pattern = metacluster_grp["seqlets_to_patterns_result"]["patterns"][pattern_name]
            pattern_scores = pattern["task0_contrib_scores"]["fwd"][:]

            meme_out.write(f"MOTIF {metacluster_name}:{pattern_name}\n")
            meme_out.write(f"letter-probability matrix: alength= 4 w= {pattern_scores.shape[0]} nsites= {len(pattern['seqlets_and_alnmts']['seqlets'])}\n")

            # remove negative weights
            positive = np.where(pattern_scores > 0, pattern_scores, 0)
            # convert to PFM
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', r'invalid value encountered in true_divide')
                pfm_plus = positive / np.expand_dims(positive.sum(axis = 1), axis = 1)
            # drop NaNs to background frequency
            pfm_plus[np.isnan(pfm_plus)] = 0.25
            # get rid of negative zeros
            pfm_plus[pfm_plus==0.] = 0.
            for position in pfm_plus:
                meme_out.write(" " + "  ".join([f"{e:.6f}" for e in position]) + "\n")

            meme_out.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = sys.argv[0],
        description = "Extract seqlets from TF-MoDISco run to minimal MEME format"
    )

    parser.add_argument(
        "config",
        help = "configuration name",
    )

    parser.add_argument(
        "species",
        help = "species to test"
    )

    args = parser.parse_args()

    main(args.config, args.species)

