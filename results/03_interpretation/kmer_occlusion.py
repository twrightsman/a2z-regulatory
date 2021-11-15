#!/usr/bin/env python

import argparse
import hashlib
import io
import json
import logging
from pathlib import Path
import subprocess
import sys
import tempfile
import warnings

from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate
import scipy.stats
import tensorflow as tf

from a2z.data import IntervalDataset, seq_one_hot_encode
from a2z.models import allow_growth
from a2z.utils import sliding_window_starts, q_values


def motif_scan(effects: pd.DataFrame, motifs_path: Path = Path("../../data/JASPAR/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt")) -> pd.DataFrame:
    with tempfile.TemporaryFile() as seq_file:
        for i, seq in enumerate(effects.index):
            seq_file.write(f">k{i}\n".encode('utf-8'))
            seq_file.write(f"{seq}\n".encode('utf-8'))

        fimo = subprocess.Popen(
            ['fimo', '--max-strand', '--skip-matched-sequence', str(motifs_path), '/dev/stdin'],
            stdin = seq_file,
            encoding = 'utf-8',
            stdout = subprocess.PIPE
        )

    motif_scan_results = pd.read_table(
        io.StringIO(fimo.stdout.read()),
        comment = '#',
        header = 0,
        usecols = ['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p-value'],
        dtype = {
            'motif_id': 'category',
            'motif_alt_id': 'category',
            'sequence_name': 'category',
            'start': int,
            'stop': int,
            'strand': 'category',
            'score': float,
            'p-value': float
        }
    )

    return motif_scan_results


def do_fisher_exact(row, n_high, n_null):
    contingency_table = np.array([
        [row.kmer_matches_high, n_high - row.kmer_matches_high],
        [row.kmer_matches_null, n_null - row.kmer_matches_null]
    ])

    odds_ratio, p_value = scipy.stats.fisher_exact(contingency_table)
    return p_value
    if p_value < (alpha / n_tests) and ((row.kmer_matches_high / n_high) > (row.kmer_matches_null / n_null)):
        return True
    return False


def main(config_name: str, species: str, k: int, skip_fimo: bool):
    with open(f"../01_models/configs/{config_name}.json") as config_file:
        config = json.load(config_file)

    # grab a held-out model
    model = tf.keras.models.load_model(Path("../01_models/tmp/results") / config_name / species.replace(" ", "_") / "0/model")

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

    test_data = data.loc[data['is_test'] & (data['target'] == 1) & (data['species'] == species)].copy().reset_index()
    test_data['seq'] = test_data.apply(lambda row: str(reference[row['seqid']][row['start']:row['end']].seq), axis = 1)
    test_dataset = IntervalDataset(test_data, no_target = True)
    test_data['prediction'] = model.predict(test_dataset).flatten()

    mutagenesis_results = {
        'seqid': [],
        'start': [],
        'masked_seq': [],
        'delta_pred': []
    }

    for row in test_data.itertuples():
        seq = str(reference[row.seqid][row.start : row.end].seq)
        mutants = {'seq': []}
        for i in sliding_window_starts(0, len(seq), k, 1):
            mutagenesis_results['seqid'].append(row.seqid)
            mutagenesis_results['start'].append(row.start)
            mutagenesis_results['masked_seq'].append(seq[i:i+k])
            mutants['seq'].append(seq[:i] + ''.join(['N'] * k) + seq[i+k:])
        mutants = pd.DataFrame(mutants)
        mutant_dataset = IntervalDataset(mutants, no_target = True)
        mutant_predictions = model.predict(mutant_dataset).flatten()
        delta_pred = mutant_predictions - row.prediction
        mutagenesis_results['delta_pred'] += list(delta_pred)

    mutagenesis_results = pd.DataFrame(mutagenesis_results)

    seq_effects = mutagenesis_results.groupby('masked_seq').mean().sort_values('delta_pred')['delta_pred']

    out_dir = Path(f"tmp/kmer-occlusion") / config_name / species.replace(' ', '_') / str(k)
    out_dir.mkdir(parents = True)

    cutoff = np.quantile(seq_effects, q = 0.05)

    fig, ax = plt.subplots()
    ax.hist(seq_effects, bins = np.arange(-1, 1, 0.01))
    ax.axvline(x = 0, color = 'red')
    ax.axvline(x = cutoff, linestyle = '--')
    ax.set_xlabel("$\Delta_{Prediction}$")
    ax.set_ylabel("Count")
    fig.savefig(out_dir / 'delta_prediction_hist.png', dpi = 150)

    high_effects = seq_effects.loc[seq_effects <= cutoff]
    n_high = len(high_effects)
    zero_index = (seq_effects - 0).abs().reset_index().sort_values('delta_pred').index[0]
    radius = int(len(seq_effects) * 0.25)
    null_effects = seq_effects.iloc[zero_index - radius:zero_index + radius]
    n_null = len(null_effects)

    fig, ax = plt.subplots()
    ax.hist(high_effects, bins = np.arange(-1, 1, 0.01), color = "red", label = "High")
    ax.hist(null_effects, bins = np.arange(-1, 1, 0.01), color = "gray", label = "Null")
    ax.set_xlabel("$\Delta_{Prediction}$")
    ax.set_ylabel("Count")
    ax.legend()
    fig.savefig(out_dir / 'delta_prediction_cutoff.png', dpi = 150)

    high_effects.to_csv(
        out_dir / "high_effect_kmers.tsv",
        sep = "\t",
        index_label = "kmer"
    )

    if not skip_fimo:
        high_effect_scan_results = motif_scan(high_effects)

        unknown_high_effect_kmers = high_effects.iloc[high_effects.reset_index().index.difference(high_effect_scan_results['sequence_name'].dtype.categories.str.lstrip("k").astype(int))].index.to_list()
        with open(out_dir / 'unknown_high_effect_kmers.txt', 'w') as out_file:
            out_file.write("\n".join(unknown_high_effect_kmers))

        null_effect_scan_results = motif_scan(null_effects)

        high_effect_kmer_counts = (high_effect_scan_results.groupby(['motif_alt_id', 'sequence_name']).size() >= 1).reset_index().groupby('motif_alt_id').sum().rename(columns = {0: 'kmer_matches'})
        null_effect_kmer_counts = (null_effect_scan_results.groupby(['motif_alt_id', 'sequence_name']).size() >= 1).reset_index().groupby('motif_alt_id').sum().rename(columns = {0: 'kmer_matches'})
        joined_kmer_counts = high_effect_kmer_counts.join(null_effect_kmer_counts, lsuffix = '_high', rsuffix = '_null', how = 'outer').fillna(0, downcast = 'infer')
        joined_kmer_counts['p-value'] = joined_kmer_counts.apply(do_fisher_exact, axis = 1, n_high = n_high, n_null = n_null)
        joined_kmer_counts['q-value'] = q_values(joined_kmer_counts['p-value'])

        joined_kmer_counts.sort_values("q-value").to_csv(
            out_dir / 'motifs.tsv',
            sep = "\t"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = sys.argv[0],
        description = "Run kmer occlusion on given model configuration and test species"
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
        "-k", "--kmer-size",
        help = "kmer size to occlude (default: 10)",
        type = int,
        default = 10,
    )

    parser.add_argument(
        "--skip-fimo",
        help = "skip FIMO scanning/filtering for JASPAR-hitting k-mers",
        action = "store_true"
    )

    logging.basicConfig(
        level = logging.INFO,
        format = '[{asctime} {levelname}] {message}',
        style = '{'
    )

    allow_growth()

    args = parser.parse_args()

    main(config_name = args.config, species = args.species, k = args.kmer_size, skip_fimo = args.skip_fimo)

