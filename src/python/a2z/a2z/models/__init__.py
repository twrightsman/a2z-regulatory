"""
models package

Code related to model architecture, training, and evaluation
"""

import json
import logging
from pathlib import Path
from typing import Tuple, Optional, List

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import precision_recall_curve, classification_report, roc_auc_score, average_precision_score, f1_score
from sklearn.metrics import roc_curve as sk_roc_curve
import tensorflow as tf
import tensorflow.keras as keras
import tensorflow.keras.layers as kl

from ..plots import pr_curve, tpr_curve, roc_curve


def allow_growth():
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        try:
            # Currently, memory growth needs to be the same across GPUs
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
            logical_gpus = tf.config.list_logical_devices('GPU')
            logging.debug("{} Physical GPUs, {} Logical GPUs".format(len(gpus), len(logical_gpus)))
        except RuntimeError as e:
            # Memory growth must be set before GPUs have been initialized
            print(e)


def get_batch_saliency_map_and_prediction(batch_X: np.ndarray, model: keras.Model) -> Tuple[np.ndarray, np.ndarray]:
    """Computes the saliency map and predictions for a given batch and model.
    The shape of the first output NumPy array is the same as the input array.
    The shape of the second depends on the model output layer.

    Credit
    ------
    https://blog.ai.ovgu.de/posts/jens/2019/001_tf20_pitfalls/index.html
    https://fairyonice.github.io/Saliency-Map-with-keras-vis.html
    """
    inp = tf.convert_to_tensor(batch_X, dtype = batch_X.dtype)
    with tf.GradientTape() as tape:
        # GradientTape doesn't record operations on constants without explicitly being told to watch them
        tape.watch(inp)
        # Run a batch through the model to record operations and the predictions
        pred = model(inp)
    # compute the derivative of the prediction w.r.t. the input one-hot batch
    saliency_map = tape.gradient(pred, inp).numpy()

    return saliency_map, pred


class TrainingHistoryEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstnace(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def train_model(model: tf.keras.Model, training_data: tf.data.Dataset, validation_data: tf.data.Dataset, out_dir: Path, overwrite = False, epochs = 10, validation_freq = 5, xtra_callbacks: Optional[List] = None):
    if not out_dir.exists():
        out_dir.mkdir(parents = True)

    callbacks = [
        tf.keras.callbacks.LambdaCallback(
            on_epoch_end = lambda epoch, logs: logging.info(f"Finished epoch {epoch}")
        )
    ]
    if xtra_callbacks is not None:
        callbacks += xtra_callbacks

    model_dir = out_dir / 'model'
    if model_dir.exists() and not overwrite:
        model = tf.keras.models.load_model(model_dir)
    else:
        logging.info("Starting training")
        train_history = model.fit(
            training_data,
            validation_data = validation_data,
            validation_freq = validation_freq,
            epochs = epochs,
            callbacks = callbacks,
            workers = 1,
            verbose = 0
        )
        logging.info("Training done")

        model.save(model_dir)
        # also save in HDF5 format
        model.save(out_dir / 'model.h5')

        # save training history
        with open(out_dir / 'train_history.json', 'w') as history_out_file:
            json.dump(train_history.history, history_out_file, cls = TrainingHistoryEncoder, indent = 2)

        # plot training curves
        fig, ax = plt.subplots()
        ax.plot(train_history.epoch, train_history.history['loss'], color = 'orange', label = 'Training')
        ax.plot(np.arange(validation_freq - 1, epochs, validation_freq), train_history.history['val_loss'], color = 'lightblue', label = 'Validation')
        ax.set_xlabel('Epoch')
        ax.set_ylabel(f"Loss ({model.loss})")
        ax.set_title(f"{model.name} Training Curve")
        ax.legend()
        fig.savefig(out_dir / 'training.svg')


def validate_model(model: tf.keras.Model, validation_data: keras.utils.Sequence, out_dir: Path, target_names = None):
    if not out_dir.exists():
        out_dir.mkdir(parents = True)

    # compute predictions
    logging.info("Starting validation, computing predictions")
    y_hat = model.predict(validation_data).flatten()
    y = np.concatenate([batch[1].flatten() for batch in validation_data])
    logging.info("Finished predictions")

    # classification report
    cls_report = classification_report(y, y_hat.round(), target_names = target_names, output_dict = True)
    auROC = roc_auc_score(y, y_hat)
    auPR = average_precision_score(y, y_hat)
    # use 0.5 as cutoff for F1 score
    f1 = f1_score(y, y_hat.round())
    logging.info("Classification report done")

    # Overall PR curve
    pr_curve(y, y_hat).savefig(out_dir / 'pr.svg')
    logging.info("PR curve done")

    # Threshold versus PR curve
    tpr_curve(y, y_hat).savefig(out_dir / 'tpr.svg')
    logging.info("TPR curve done")

    # Overall ROC curve
    roc_curve(y, y_hat).savefig(out_dir / 'roc.svg')
    logging.info("ROC curve done")

    metrics = {
        'classification_report': cls_report,
        'auROC': auROC,
        'auPR': auPR,
        'f1_score': f1
    }

    with open(out_dir / 'metrics.json', 'w') as metrics_file:
        json.dump(metrics, metrics_file, indent = 2)


def build_DanQ(window_size: int = 300, kernel: np.ndarray = None, conv_activation = 'relu') -> keras.Model:
    DanQ = keras.Sequential([
        kl.Conv1D(320, 26, activation = conv_activation, input_shape = (window_size, 4)),
        kl.Dropout(0.2),
        kl.MaxPool1D(13, 13),
        kl.Bidirectional(kl.LSTM(320)),
        kl.Dropout(0.5),
        kl.Dense(925),
        kl.Dense(1, activation = 'sigmoid')
    ], name = "DanQ")

    if kernel is not None:
        # initialize first layer kernel with motifs
        DanQ.layers[0].kernel.assign(kernel)

    return DanQ


def build_Basset(window_size: int = 300) -> keras.Model:
    Basset = keras.Sequential([
        kl.Convolution1D(
            filters = 300,
            kernel_size = 19,
            input_shape = (window_size, 4)
        ),
        kl.BatchNormalization(),
        kl.Activation('relu'),
        kl.MaxPooling1D(
            pool_size = 3
        ),
        kl.Convolution1D(
            filters = 200,
            kernel_size = 11
        ),
        kl.BatchNormalization(),
        kl.Activation('relu'),
        kl.MaxPooling1D(
            pool_size = 4
        ),
        kl.Convolution1D(
            filters = 200,
            kernel_size = 7
        ),
        kl.BatchNormalization(),
        kl.Activation(
            activation = 'relu'
        ),
        kl.MaxPooling1D(
            pool_size = 4
        ),
        kl.Flatten(),
        kl.Dense(
            units = 1000
        ),
        kl.Activation('relu'),
        kl.Dropout(
            rate = 0.3
        ),
        kl.Dense(
            units = 1000
        ),
        kl.Activation('relu'),
        kl.Dropout(
            rate = 0.3
        ),
        kl.Dense(
            units = 1
        ),
        kl.Activation('sigmoid')
    ], name = "Basset")

    return Basset


def build_CharPlant(window_size: int = 300) -> keras.Model:
    CharPlant = keras.Sequential([
        kl.Convolution1D(
            filters = 200,
            kernel_size = 19,
            padding = 'same',
            input_shape = (window_size, 4)
        ),
        kl.Activation('relu'),
        kl.Dropout(
            rate = 0.6
        ),
        kl.Convolution1D(
            filters = 100,
            kernel_size = 11,
            padding = 'same'
        ),
        kl.Activation('relu'),
        kl.Dropout(
            rate = 0.6
        ),
        kl.Flatten(),
        kl.Dense(
            units = 200,
            kernel_initializer = 'normal',
            activation = 'relu'
        ),
        kl.Dropout(
            rate = 0.6
        ),
        kl.Dense(
            units = 1,
            kernel_initializer = 'normal',
            activation = 'linear'
        ),
        kl.Activation('sigmoid')
    ], name = "CharPlant")

    return CharPlant


def build_DeeperDeepSEA(window_size: int = 600, n_targets: int = 1) -> keras.Model:
    DeeperDeepSEA = keras.Sequential([
        kl.Conv1D(
            filters = 320,
            kernel_size = 8,
            input_shape = (window_size, 4)
        ),
        kl.Activation('relu'),
        kl.Conv1D(
            filters = 320,
            kernel_size = 8
        ),
        kl.Activation('relu'),
        kl.MaxPool1D(
            pool_size = 4,
            strides = 4
        ),
        kl.BatchNormalization(),
        kl.Conv1D(
            filters = 480,
            kernel_size = 8
        ),
        kl.Activation('relu'),
        kl.Conv1D(
            filters = 480,
            kernel_size = 8
        ),
        kl.Activation('relu'),
        kl.MaxPool1D(
            pool_size = 4,
            strides = 4
        ),
        kl.BatchNormalization(),
        kl.Dropout(rate = 0.2),
        kl.Conv1D(
            filters = 960,
            kernel_size = 8
        ),
        kl.Activation('relu'),
        kl.Conv1D(
            filters = 960,
            kernel_size = 8
        ),
        kl.Activation('relu'),
        kl.BatchNormalization(),
        kl.Dropout(rate = 0.2),
        kl.Flatten(),
        kl.Dense(n_targets),
        kl.Activation('relu'),
        kl.BatchNormalization(),
        kl.Dense(n_targets, activation = 'sigmoid')
    ], name = "DeeperDeepSEA")

    return DeeperDeepSEA


def build_model(model_name: str, **kwargs) -> keras.Model:
    if model_name == 'DanQ':
        window_size = kwargs.get('window_size', 300)
        kernel = kwargs.get('kernel', None)
        conv_activation = kwargs.get('conv_activation', 'relu') 
        return build_DanQ(window_size = window_size, kernel = kernel, conv_activation = conv_activation)
    elif model_name == 'Basset':
        window_size = kwargs.get('window_size', 300)
        return build_Basset(window_size = window_size)
    elif model_name == 'CharPlant':
        window_size = kwargs.get('window_size', 300)
        return build_CharPlant(window_size = window_size)
    elif model_name == 'DeeperDeepSEA':
        # cannot be built with window size 300
        window_size = kwargs.get('window_size', 600)
        return build_DeeperDeepSEA(window_size = window_size)
    else:
        raise ValueError("Unknown model '{}'".format(model_name))
