#!/usr/bin/env python

import hashlib
import json
import logging
from pathlib import Path
import sys

import pandas as pd

CLADEMAP = {
    'Arabidopsis thaliana': {'Angiosperm', 'Dicot'},
    'Asparagus officinalis': {'Angiosperm', 'Monocot'},
    'Brachypodium distachyon': {'Angiosperm', 'Monocot', 'Grass'},
    'Eutrema salsugineum': {'Angiosperm', 'Dicot'},
    'Glycine max': {'Angiosperm', 'Dicot'},
    'Hordeum vulgare': {'Angiosperm', 'Monocot', 'Grass'},
    'Oryza sativa': {'Angiosperm', 'Monocot', 'Grass'},
    'Phaseolus vulgaris': {'Angiosperm', 'Dicot'},
    'Populus trichocarpa': {'Angiosperm', 'Dicot'},
    'Setaria viridis': {'Angiosperm', 'Monocot', 'Grass'},
    'Sorghum bicolor': {'Angiosperm', 'Monocot', 'Grass'},
    'Spirodela polyrhiza': {'Angiosperm', 'Monocot'},
    'Vitis vinifera': {'Angiosperm', 'Dicot'},
    'Zea mays': {'Angiosperm', 'Monocot', 'Grass'}
}

logging.basicConfig(level = logging.DEBUG)

if len(sys.argv) != 2:
    logging.error("Usage: %s config.json" % sys.argv[0])
    sys.exit(1)

config_path = Path(sys.argv[1])
with open(config_path) as config_file:
    config = json.load(config_file)

cmdline = 'mkdir -p "{out_dir}" && TF_CPP_MIN_LOG_LEVEL=1 python train.py -v -v "{config_path}" "{out_dir}" "{test_species}" "{train_species}" 2> "{out_dir}/train.err" > "{out_dir}/train.out"'

split_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n" + json.dumps(config['splitting'], indent = 2) + "\n").encode('utf-8')).hexdigest()
config_hash = hashlib.sha256((json.dumps(config, indent = 2) + "\n").encode('utf-8')).hexdigest()
samples = config['training']['samples']

# subset species to within clade
clade = config['training']['clade']
clade_set = set(filter(lambda s: clade in CLADEMAP[s], CLADEMAP.keys()))

split_data_path = f"tmp/data/split/{split_hash}.tsv"
available_species = set(pd.read_table(split_data_path, header = 0, usecols = ['species'])['species'].unique())

available_clade_species = available_species & clade_set
if not available_clade_species:
    logging.error("No species in the data (%s) are in the clade '%s'" % (split_data_path, clade))
    sys.exit(1)

if config['training']['strategy'] == 'drop-one':
    # train on all other species in given clade
    for test_species in available_clade_species:
        for sample in range(samples):
            out_dir = Path(f"tmp/results/{config_path.stem}/{test_species.replace(' ', '_')}/{sample}")
            if out_dir.exists():
                logging.warning("%s exists, skipping that job" % out_dir)
            else:
                sys.stdout.write(cmdline.format(
                    config_path = config_path,
                    test_species = test_species,
                    train_species = ','.join(available_clade_species - {test_species}),
                    out_dir = out_dir
                ) + "\n")
elif config['training']['strategy'] == 'within':
    # train within species in given clade
    for test_species in available_clade_species:
        for sample in range(samples):
            out_dir = Path(f"tmp/results/{config_path.stem}/{test_species.replace(' ', '_')}/{sample}")
            if (out_dir / 'metrics.json').exists():
                logging.warning("%s/metrics.json exists, skipping that job" % out_dir)
            else:
                sys.stdout.write(cmdline.format(
                    config_path = config_path,
                    test_species = test_species,
                    train_species = test_species,
                    out_dir = out_dir
                ) + "\n")
else:
    logging.error("Unknown training strategy '%s'" % config['training']['strategy'])
    sys.exit(1)

