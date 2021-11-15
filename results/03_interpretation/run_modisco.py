#!/usr/bin/env python

import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys

from Bio import SeqIO
import h5py
import matplotlib.pyplot as plt
import modisco
import numpy as np
import pandas as pd
import tensorflow as tf

from a2z.data import IntervalDataset, seq_one_hot_encode
from a2z.models import allow_growth, get_batch_saliency_map_and_prediction


def main(config_name: str, species: str) -> None:
    with open(f"../01_models/configs/{config_name}.json") as config_file:
        config = json.load(config_file)

    # grab a held-out model
    model = tf.keras.models.load_model(Path("../01_models/tmp/results") / config_name / species.replace(" ", "_") / "0" / "model")

    reference = {rec.id: rec for rec in SeqIO.parse(Path(config['paths']['genomes']) / f"{species.replace(' ', '_')}.fa", "fasta")}

    split_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n" + json.dumps(config['splitting'], indent = 2) + "\n").encode('utf-8')).hexdigest()
    split_data_path = f"../01_models/tmp/data/split/{split_hash}.tsv"

    data = pd.read_table(
        split_data_path,
        header = 0,
        dtype = {
            'seqid': 'category',
            'start': int,
            'end': int,
            'target': np.byte,
            'distance_class': 'category',
            'species': 'category',
            'is_test': bool,
            'is_val': bool,
            'is_train': bool
        }
    )

    test_data = data.loc[data['is_test'] & (data['species'] == species)].copy().reset_index()
    test_data['seq'] = test_data.apply(lambda row: str(reference[row['seqid']][row['start']:row['end']].seq), axis = 1)

    # generate contribution scores and hypothetical contribution scores (~45 sec)
    contrib_scores = []
    hyp_contrib_scores = []
    onehot_data = []
    test_dataset = IntervalDataset(test_data, no_target = True)
    for batch_idx in range(len(test_dataset)):
        inp = test_dataset[batch_idx][0]
        saliency_maps, predictions = get_batch_saliency_map_and_prediction(inp, model)
        contrib_scores.append(saliency_maps * inp)
        # mean-normalize gradients
        hyp_contrib_scores.append(saliency_maps - saliency_maps.mean(axis=2)[:, :, np.newaxis])
        onehot_data.append(inp)
    contrib_scores = np.concatenate(contrib_scores)
    hyp_contrib_scores = np.concatenate(hyp_contrib_scores)

    onehot_data = np.concatenate(onehot_data)
    task_to_scores = {'task0': contrib_scores}
    task_to_hyp_scores = {'task0': hyp_contrib_scores}

    null_per_pos_scores = modisco.coordproducers.LaplaceNullDist(num_to_samp = 5000)
    tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
        sliding_window_size = 15,
        flank_size = 5,
        target_seqlet_fdr = 0.15,
        seqlets_to_patterns_factory =
            modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
                trim_to_window_size = 15,
                initial_flank_to_add = 5,
                final_min_cluster_size = 60,
                n_cores = 10
            )
    )(
        task_names = ['task0'],
        contrib_scores = task_to_scores,
        hypothetical_contribs = task_to_hyp_scores,
        one_hot = onehot_data,
        null_per_pos_scores = null_per_pos_scores
    )

    out_path = Path("tmp/tf-modisco") / config_name / f"{species.replace(' ', '_')}.hdf5"
    if not out_path.parent.exists():
        out_path.parent.mkdir(parents = True, exist_ok = True)

    results_h5file = h5py.File(out_path, "w")
    tfmodisco_results.save_hdf5(results_h5file)
    results_h5file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = sys.argv[0],
        description = "Run TF-MoDISco on a given model configuration and species"
    )

    parser.add_argument(
        "config",
        help = "configuration name",
    )

    parser.add_argument(
        "species",
        help = "species to test"
    )

    logging.basicConfig(
        level = logging.INFO,
        format = '[{asctime} {levelname}] {message}',
        style = '{'
    )

    args = parser.parse_args()

    allow_growth()

    main(args.config, args.species)

