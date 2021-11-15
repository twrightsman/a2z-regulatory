#!/usr/bin/env python
import hashlib
import json
import logging
from pathlib import Path
import random
import sys

import numpy as np
from numpy.random import SeedSequence, PCG64
import pandas as pd
from pyfaidx import Fasta

# set seed for reproducibility
entropy = 45572032134201216385363710516788190149
seed = SeedSequence(entropy)
rng = PCG64(seed)
pyrng = random.Random(x = 4242)

logging.basicConfig(
    format="{asctime} [{module}:{levelname}] {message}",
    style="{",
    level=logging.INFO
)


# get command-line arguments
if len(sys.argv) != 2:
    logging.error("Usage: %s config.json" % sys.argv[0])
    sys.exit(1)

# read config
with open(sys.argv[1]) as config_file:
    config = json.load(config_file)

window_size = config['preprocessing']['window_size']
data_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n").encode('utf-8')).hexdigest()
all_data_path = Path('tmp/data') / f"{data_hash}.tsv"

split_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n" + json.dumps(config['splitting'], indent = 2) + "\n").encode('utf-8')).hexdigest()
out_dir = all_data_path.parent / "split"
out_path = out_dir / f"{split_hash}.tsv"
if out_path.exists():
    logging.error("%s exists, exiting" % out_path)
    sys.exit(1)
if not out_dir.exists():
    out_dir.mkdir()

# read data
names = ['seqid', 'start', 'end', 'target', 'species']
if 'annotations' in config['paths']:
    names.insert(4, 'distance_class')
all_data = pd.read_table(
    all_data_path,
    names = names,
    dtype = {
        'seqid': 'category',
        'start': int,
        'end': int,
        'target': np.byte,
        'distance_class': 'category',
        'species': 'category'
    }
)

test_frac = config['splitting']['test_frac']
all_data['is_test'] = 0
all_data['is_val'] = 0
if ('chromosome' in config['splitting']) and config['splitting']['chromosome']:
    logging.info("Using hold-out chromosome as test and validation")
    # use hold-out chromosomes as test and validation sets
    indices = None
    for species in all_data['species'].unique():
        reference = Fasta(str(Path(config['paths']['genomes']) / f"{species.replace(' ', '_')}.fa"), rebuild = False)
        # grab random chromosome bigger than 1,000,000bp and with at least 5 open chromatin regions (minimal sanity checks)
        chromosomes = list(reference.keys())
        pyrng.shuffle(chromosomes)
        i = -1
        while (len(reference[chromosomes[i]]) < 1_000_000) or (len(all_data[(all_data['species'] == species) & (all_data['seqid'] == chromosomes[i]) & (all_data['target'] == 1)]) < 5):
            i += 1
            if i >= len(chromosomes):
                logging.error("Couldn't find a suitable hold-out chromosome for %s" % species)
                sys.exit(1)
        test_chr = chromosomes[i]
        logging.info("Selected '%s' as hold-out chromosome for %s" % (test_chr, species))
        if indices is None:
            indices = all_data[(all_data['species'] == species) & (all_data['seqid'] == test_chr)].sample(frac = 1, random_state = rng).index
        else:
            indices = indices.append(all_data[(all_data['species'] == species) & (all_data['seqid'] == test_chr)].sample(frac = 1, random_state = rng).index)
    # shuffle the indices after concatenting the sampled species indices
    indices = indices.to_series().sample(frac = 1, random_state = rng).index
else:
    # sample equal-sized test set and validation set
    indices = all_data.sample(frac = test_frac * 2, random_state = rng).index
midpoint = len(indices) // 2
all_data.loc[indices.to_series().iloc[:midpoint].index, 'is_test'] = 1
all_data.loc[indices.to_series().iloc[midpoint:].index, 'is_val'] = 1

# downsample training negatives for each species
balancing = config['splitting']['balancing']
all_data['is_train'] = 0
if balancing:
    training_data = all_data.loc[(all_data['is_test'] == 0) & (all_data['is_val'] == 0)]
    smallest_class = training_data.groupby(balancing).size().min()
    for balance_class, subdata in training_data.groupby(balancing):
        downsample = subdata.sample(n = smallest_class, random_state = rng).index
        all_data.loc[downsample, 'is_train'] = 1
else:
    all_data['is_train'] = (all_data['is_test'] == 0) & (all_data['is_val'] == 0)

# filter anything not in any of the three sets and write to file
all_data.loc[(all_data['is_test'] == 1) | (all_data['is_val'] == 1) | (all_data['is_train'] == 1)].to_csv(
    out_path,
    sep = "\t",
    header = True,
    index = False
)

