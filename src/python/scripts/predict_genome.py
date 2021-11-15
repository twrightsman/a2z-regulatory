#!/usr/bin/env python3
from Bio import SeqIO
import numpy as np
import tensorflow as tf
import tensorflow.keras as k

import argparse
import logging
import math
from pathlib import Path
import sys
from typing import Optional

from a2z.data import seq_one_hot_encode
from a2z.models import allow_growth


def sliding_window_starts(start, end, length, stride):
    """ Generator for start of sliding intervals that are fully-enclosed within given interval """
    # skip short intervals
    if end - start >= length:
        yield from range(start, end - length + 1, stride)


class BatchedOneHotSlidingWindows(k.utils.Sequence):
    def __init__(self, sequence: str, length: int, stride: int = 1, batch_size: int = 128):
        self.sequence = sequence
        self.length = length
        self.stride = stride
        self.batch_size = batch_size

    def __len__(self):
        return math.ceil((math.floor((len(self.sequence) - self.length) / self.stride) + 1) / self.batch_size)

    def __getitem__(self, idx):
        start = (idx * self.batch_size) * self.stride
        end = (((idx + 1) * self.batch_size) * self.stride) + self.length - 1
        end = min(end, len(self.sequence))
        x = np.stack(list(map(seq_one_hot_encode, (self.sequence[s : s + self.length] for s in sliding_window_starts(start, end, self.length, self.stride)))))
        return (x, None)


def predict_genome(sequence_path: Path, model_path: Path, stride: Optional[int], batch_size: int) -> None:
    allow_growth()

    model = k.models.load_model(model_path)
    logging.info("Loaded model")
    if len(model.inputs) != 1:
        logging.error("Model has incorrect number of input layers (!=1)")
        sys.exit(1)
    elif len(model.inputs[0].shape) != 3:
        logging.error("Model input layer has incorrect number of dimensions (!=3)")

    window_size = model.inputs[0].shape[1]
    logging.debug(f"Model has window size of {window_size}bp")
    if stride is None:
        stride = window_size

    records = {rec.id: rec for rec in SeqIO.parse(sequence_path, 'fasta')}
    logging.info("Loaded sequence")

    for seqid, record in records.items():
        if len(record) < window_size:
            logging.warning("Skipping record '{}' because it was too short (<{} bp)".format(seqid, window_size))
        else:
            inputs = BatchedOneHotSlidingWindows(str(record.seq), length = window_size, stride = stride, batch_size = batch_size)
            logging.info("Starting predictions on '{}'".format(seqid))
            outputs = model.predict(inputs)
            logging.info("Predictions done")
            for i, start in enumerate(sliding_window_starts(0, len(record), length = window_size, stride = stride)):
                sys.stdout.write("{}\t{}\t{}\t{}\n".format(seqid, start, start + window_size, "\t".join([str(output) for output in outputs[i]])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Predict model outputs in sliding windows along a sequence and output BED3+N format to stdout",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="output progress and other informative messages",
        action="count",
        dest="verbosity",
        default=0
    )

    parser.add_argument(
        "model",
        type = Path,
        help = "path to tf.keras model directory"
    )

    parser.add_argument(
        "sequence",
        type = Path,
        help = "path to FASTA"
    )

    parser.add_argument(
        "--stride",
        type = int,
        default = None,
        help = "stride for the sliding windows (default: the window length, non-overlapping)"
    )

    parser.add_argument(
        "--batch-size",
        type = int,
        default = 128,
        help = "batch size to use for predictions (default: 128)"
    )

    args = parser.parse_args()

    logging.basicConfig(
        format="{asctime} [{module}:{levelname}] {message}",
        style="{",
        level=max(logging.DEBUG, logging.WARNING - (args.verbosity * 10))
    )

    predict_genome(
        sequence_path = args.sequence,
        model_path = args.model,
        stride = args.stride,
        batch_size = args.batch_size
    )

