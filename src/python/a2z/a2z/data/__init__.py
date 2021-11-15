"""
data
"""

import logging
from math import ceil
from typing import Optional, Dict

from Bio import Seq
import numpy as np
import pandas as pd
from tensorflow import keras as k


BASE_1HOT = {
    "A": np.array([1, 0, 0, 0]),
    "C": np.array([0, 1, 0, 0]),
    "G": np.array([0, 0, 1, 0]),
    "T": np.array([0, 0, 0, 1]),
    "W": np.array([0.5, 0, 0, 0.5]),
    "S": np.array([0, 0.5, 0.5, 0]),
    "M": np.array([0.5, 0.5, 0, 0]),
    "K": np.array([0, 0, 0.5, 0.5]),
    "R": np.array([0.5, 0, 0.5, 0]),
    "Y": np.array([0, 0.5, 0, 0.5]),
    "B": np.array([0, 1.0 / 3, 1.0 / 3, 1.0 / 3]),
    "D": np.array([1.0 / 3, 0, 1.0 / 3, 1.0 / 3]),
    "H": np.array([1.0 / 3, 1.0 / 3, 0, 1.0 / 3]),
    "V": np.array([1.0 / 3, 1.0 / 3, 1.0 / 3, 0]),
    "N": np.array([0.25, 0.25, 0.25, 0.25]),
}


def seq_one_hot_encode(seq: str, pad_to: Optional[int] = None):
    """ one-hot encodes a DNA sequence """
    if pad_to:
        encoded = np.zeros(shape=(pad_to, 4))
        if len(seq) > pad_to:
            seq = seq[:pad_to]
    else:
        encoded = np.zeros(shape=(len(seq), 4))

    for i, base in enumerate(seq.upper()):
        try:
            encoded[i, :] = BASE_1HOT[base]
        except KeyError:
            logging.error(
                f"Unrecognized base encountered during one-hot " f"encoding: '{base}'"
            )
    return encoded


class IntervalDataset(k.utils.Sequence):
    def __init__(self, intervals_df: pd.DataFrame, batch_size = 128, pad_to = None, no_target = False):
        self._df = intervals_df
        self.batch_size = batch_size
        self.pad_to = pad_to
        self.no_target = no_target

    def __len__(self):
        return ceil(len(self._df) / self.batch_size)

    def __getitem__(self, idx: int):
        batch_start_idx = idx * self.batch_size
        batch_intervals = self._df.iloc[batch_start_idx : batch_start_idx + self.batch_size]
        try:
            x = np.stack([seq_one_hot_encode(seq, self.pad_to) for seq in batch_intervals['seq']])
        except ValueError as e:
            logging.error("ValueError encountered while trying to one-hot encode sequence batches, are all of your sequences the same length?")
            raise e
        if self.no_target:
            y = None
        else:
            y = np.expand_dims(batch_intervals['target'].to_numpy(), axis = 1)

        return (x, y)


class IntervalReferenceDataset(k.utils.Sequence):
    def __init__(self, intervals_df: pd.DataFrame, references: Dict[str, Seq.Seq], batch_size = 128, pad_to = None, no_target = False, random_shift: bool = False, species_weights: Optional[Dict[str, float]] = None):
        self._df = intervals_df
        self._references = references
        self.batch_size = batch_size
        self.pad_to = pad_to
        self.no_target = no_target
        self.random_shift = random_shift
        self.species_weights = species_weights
        if self.random_shift:
            self._rng = np.random.default_rng()

    def __len__(self):
        return ceil(len(self._df) / self.batch_size)

    def __getitem__(self, idx: int):
        batch_start_idx = idx * self.batch_size
        batch_intervals = self._df.iloc[batch_start_idx : batch_start_idx + self.batch_size]

        if self.random_shift:
            # draw shift values from a normal distribution, clip between [-5, 5]
            shifts = self._rng.normal(size = (self.batch_size,)).round().clip(-5, 5).astype(np.int8)
        else:
            shifts = np.zeros(shape = (self.batch_size,), dtype = np.int8)

        try:
            seqs = [self._references[row.species][row.seqid][row.start+shifts[i]:row.end+shifts[i]] for i, row in enumerate(batch_intervals.itertuples())]
            x = np.stack([seq_one_hot_encode(seq, self.pad_to) for seq in seqs])
        except ValueError as e:
            logging.error("ValueError encountered while trying to one-hot encode sequence batches, are all of your sequences the same length?")
            raise e
        if self.no_target:
            y = None
        else:
            y = np.expand_dims(batch_intervals['target'].to_numpy(), axis = 1)

        if self.species_weights is not None:
            weights = np.array([self.species_weights[row.species][row.target] for row in batch_intervals.itertuples()])
            return (x, y, weights)
        return (x, y)
