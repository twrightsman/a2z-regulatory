#!/usr/bin/env python

import argparse
import functools
import hashlib
import json
import logging
from multiprocessing import Pool
from pathlib import Path
import random
from statistics import mean
import sys

from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import pandas as pd
import tensorflow as tf

from a2z.data import seq_one_hot_encode
from a2z.data.motifs import MinimalMEME, Motif
from a2z.models import allow_growth
from a2z.models.attribution import global_importance_analysis


def motif_pGIA(motif: Motif, model: tf.keras.Model, null_set: np.ndarray, null_set_preds: np.ndarray) -> (list, list):
    motif_consensus_1hot = np.eye(4)[motif.pfm.argmax(axis = 1)]
    motif_consensus_1hot_rev = np.flip(motif_consensus_1hot, axis = [0, 1])

    gia = functools.partial(
        global_importance_analysis,
        motif = motif_consensus_1hot,
        model = model,
        null_set = null_set,
        null_set_preds = null_set_preds
    )
    gia_rev = functools.partial(
        global_importance_analysis,
        motif = motif_consensus_1hot_rev,
        model = model,
        null_set = null_set,
        null_set_preds = null_set_preds
    )

    positional_gia = []
    positional_gia_rev = []
    for pos in range(0, null_set.shape[1] - len(motif)):
        positional_gia.append(gia(embed_location = pos))
        positional_gia_rev.append(gia_rev(embed_location = pos))

    return positional_gia, positional_gia_rev


def main(config_name: str, species: str, meme_file_path: Path, null_file_path: Path):
    # grab a held-out model
    model_path = Path("../01_models/tmp/results") / config_name / species.replace(" ", "_") / "0/model"
    model = tf.keras.models.load_model(model_path)

    # load in null profile
    nulls = pd.read_table(
        null_file_path,
        header = 0,
        dtype = {
            'seq': str,
            'pred': float
        }
    )
    null_set = np.stack(nulls['seq'].map(seq_one_hot_encode))
    null_set_preds = nulls['pred'].to_numpy()

    # load in JASPAR motif
    motif = list(MinimalMEME(meme_file_path).motifs.values())[0]

    positional_gia, positional_gia_rev = motif_pGIA(motif, model, null_set, null_set_preds)
    positional_gia_both = positional_gia + positional_gia_rev

    print("\t".join([motif.name if motif.name is not None else motif.identifier] + [str(f) for f in positional_gia_both]))


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

    parser.add_argument(
        "null_file_path",
        help = "path to TSV file of null sequences and their predictions",
        type = Path
    )

    parser.add_argument(
        "meme_path",
        help = "path to MEME file to test global importance",
        type = Path
    )

    logging.basicConfig(
        level = logging.INFO,
        format = '[{asctime} {levelname}] {message}',
        style = '{'
    )

    allow_growth()

    args = parser.parse_args()

    main(
        config_name = args.config,
        species = args.species,
        meme_file_path = args.meme_path,
        null_file_path = args.null_file_path
    )
