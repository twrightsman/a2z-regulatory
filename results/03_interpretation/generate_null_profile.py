#!/usr/bin/env python

import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from pyfaidx import Fasta
import tensorflow as tf

from a2z.data import IntervalDataset, seq_one_hot_encode
from a2z.models import allow_growth
from a2z.models.attribution import generate_null_profile


NUCLEOTIDES = np.array(["A", "C", "G", "T"])
NUCLEOTIDES.flags.writeable = False


def reverse_1hot(seq_1hot: np.ndarray, alphabet: np.ndarray = NUCLEOTIDES) -> str:
     return "".join(alphabet[np.argmax(seq_1hot, axis = 1)])


def main(config_name: str, species: str):
    with open(f"../01_models/configs/{config_name}.json") as config_file:
            config = json.load(config_file)

    # grab a held-out model
    model_path = Path("../01_models/tmp/results") / config_name / species.replace(" ", "_") / "0/model"
    model = tf.keras.models.load_model(model_path)

    reference = Fasta(str(Path(config['paths']['genomes']) / f"{species.replace(' ', '_')}.fa"), as_raw = True, rebuild = False)

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

    test_data = data.loc[data['is_test'] & (data['target'] == 1) & (data['species'] == species)].copy().reset_index()
    test_data['seq'] = test_data.apply(lambda row: str(reference[row['seqid']][row['start']:row['end']]), axis = 1)
    test_dataset = IntervalDataset(test_data, no_target = True)

    # generate null sequences from profile model
    onehot_seqs = np.concatenate([batch[0] for batch in test_dataset])
    null_size = 1000

    null_set, null_set_preds = generate_null_profile(onehot_seqs, n = null_size, model = model)

    pd.DataFrame({'seq': map(reverse_1hot, null_set), 'pred': null_set_preds}).to_csv(
        sys.stdout,
        sep = "\t",
        index = False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = sys.argv[0],
        description = "Run Global Importance Analysis on given MEME file by sampling kmers"
    )

    parser.add_argument(
        "config",
        help = "configuration name"
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

    allow_growth()

    args = parser.parse_args()

    main(config_name = args.config, species = args.species)
