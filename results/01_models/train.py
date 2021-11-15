#!/usr/bin/env python3

import argparse
import datetime
import hashlib
import json
import logging
from pathlib import Path
import sys
from typing import Optional

import pandas as pd
from pyfaidx import Fasta
import numpy as np
import tensorflow as tf

from a2z.data import IntervalReferenceDataset
from a2z.data.motifs import MinimalMEME
from a2z.models import allow_growth, train_model, validate_model, build_model


def train_and_validate(config, out_dir, train_species, test_species, intervals, genome_dir):
    logging.info("Starting trainable")

    with open(out_dir / 'config.json', 'w') as config_file:
        json.dump(config, config_file, indent = 2)

    references = {}
    for species in (test_species | train_species):
        references[species] = Fasta(str(genome_dir / f"{species.replace(' ', '_')}.fa"), as_raw = True, rebuild = False)

    # split training and test data
    model_name = config['training']['architecture']

    # set up training data
    training_intervals = intervals.loc[(intervals['species'].isin(train_species)) & (intervals['is_train'] == 1), :].sample(frac = 1)
    validation_intervals = intervals.loc[(intervals['species'].isin(train_species)) & (intervals['is_val'] == 1), :]
    test_intervals = intervals.loc[(intervals['species'].isin(test_species)) & (intervals['is_test'] == 1), :]
    logging.info("Set up model-specific data")

    # set up species weights
    if ('negative_upweighting' in config) and config['negative_upweighting']:
        test_class_counts = intervals.groupby(['species', 'target'])['is_test'].sum().astype(np.int64)
        species_weights = {}
        for species in intervals['species'].unique():
            species_weights[species] = {}
            species_weights[species][1] = 1
            species_weights[species][0] = test_class_counts[species][0] / test_class_counts[species][1]
        logging.debug("Species weights %s" % species_weights)
    else:
        species_weights = None

    training_data = IntervalReferenceDataset(training_intervals, references, species_weights = species_weights)
    validation_data = IntervalReferenceDataset(validation_intervals, references)
    test_data = IntervalReferenceDataset(test_intervals, references)

    window_size = config['preprocessing']['window_size']

    logging.debug(f"Using window size of {window_size}")
    model = build_model(model_name, window_size = window_size, **config['training']['hyperparameters'])
    model.compile(optimizer = 'Adam', loss = 'BCE')
    logging.info("Built datasets and compiled the model")

    training_start = datetime.datetime.now()
    train_model(
        model,
        training_data,
        validation_data,
        out_dir,
        overwrite = True,
        epochs = config['training']['epochs'],
        validation_freq = 5,
        xtra_callbacks = [
            tf.keras.callbacks.TensorBoard(
                log_dir = out_dir / 'tensorboard',
                write_graph = False
            )
        ]
    )
    logging.info("Trained the model in %s" % (datetime.datetime.now() - training_start))

    validation_start = datetime.datetime.now()
    validate_model(model, test_data, out_dir)
    logging.info("Validated the model in %s" % (datetime.datetime.now() - validation_start))


def main(config, outdir, test_species, train_species):
    split_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n" + json.dumps(config['splitting'], indent = 2) + "\n").encode('utf-8')).hexdigest()
    logging.info("Loading data")
    intervals = pd.read_table(
        f"tmp/data/split/{split_hash}.tsv",
        header = 0,
        dtype = {
            'seqid': str,
            'start': int,
            'end': int,
            'distance_class': 'category',
            'species': str,
            'target': np.byte,
            'is_train': np.byte,
            'is_test': np.byte,
            'is_val': np.byte
        }
    )
    logging.info("Data loaded")

    genome_dir = Path(config['paths']['genomes'])
    train_and_validate(config, outdir, train_species, test_species, intervals, genome_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Train a single model configuration"
    )

    parser.add_argument(
        "config_path",
        help = "path to config",
        type = Path,
    )

    parser.add_argument(
        "outdir",
        help = "path to output directory",
        type = Path,
    )

    parser.add_argument(
        "test_species",
        help = "list of test species"
    )

    parser.add_argument(
        "train_species",
        help = "list of train species"
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="output progress and other informative messages",
        action="count",
        dest="verbosity",
        default=0,
    )

    args = parser.parse_args()

    logging.basicConfig(
        format="{asctime} [{module}:{levelname}] {message}",
        style="{",
        level=max(logging.DEBUG, logging.WARNING - (args.verbosity * 10)),
        filename = Path(args.outdir) / "train.log"
    )

    # configure GPU for shared training
    allow_growth()

    with open(args.config_path) as config_file:
      config = json.load(config_file)

    main(
        config = config,
        outdir = args.outdir,
        test_species = set(args.test_species.split(',')),
        train_species = set(args.train_species.split(','))
    )

