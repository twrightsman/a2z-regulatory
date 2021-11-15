"""
a2z.plots

Functions to return various plots for evaluating models
"""

from collections import Counter
from typing import Iterable

import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, f1_score
from sklearn.metrics import roc_curve as sk_roc_curve


def pr_curve(y_true: Iterable, scores: Iterable) -> Figure:
    fig, ax = plt.subplots()

    # determine baseline
    true_cnt = Counter(y_true)
    n = sum(true_cnt.values())
    n_true = true_cnt[1]
    base_y = n_true / n

    precision, recall, thresholds = precision_recall_curve(y_true, scores)
    ax.plot(recall, precision)
    ax.plot([0, 1], [base_y] * 2, color="gray", linestyle="--", label="Baseline")
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend()

    return fig


def tpr_curve(y_true: Iterable, scores: Iterable) -> Figure:
    """ F1 baseline credit: https://stats.stackexchange.com/questions/390200/what-is-the-baseline-of-the-f1-score-for-a-binary-classifier  """
    fig, ax = plt.subplots()

    # determine baseline
    true_cnt = Counter(y_true)
    n = sum(true_cnt.values())
    n_true = true_cnt[1]
    base_y = n_true / n

    precision, recall, thresholds = precision_recall_curve(y_true, scores)
    f1 = 2 * ((precision * recall) / (precision + recall))
    f1_baseline = (2 * base_y) / (base_y + 1)
    ax.plot(thresholds, precision[:-1], label = "Precision")
    ax.plot(thresholds, recall[:-1], label = "Recall")
    ax.plot(thresholds, f1[:-1], label = "F1 Score")
    ax.plot([0, 1], [base_y] * 2, color="gray", linestyle="--", label="Precision Baseline")
    ax.plot([0, 1], [f1_baseline] * 2, color="lightgray", linestyle = "--", label="F1 Baseline")
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_xlabel("Threshold")
    ax.set_ylabel("Metric")
    ax.legend()

    return fig


def roc_curve(y_true: Iterable, scores: Iterable) -> Figure:
    fig, ax = plt.subplots()

    fpr, tpr, thresholds = sk_roc_curve(y_true, scores)
    ax.plot(fpr, tpr)
    ax.plot([0, 1], [0, 1], color="gray", linestyle="--", label="Baseline")
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.legend()

    return fig


def make_placeholder(fig: Figure) -> None:
    """ Adds a stamp over the plot signalling it is a placeholder plot """
    fig.add_artist(FancyBboxPatch(
        xy = (0.35, 0.45),
        width = 0.3,
        height = 0.1,
        boxstyle = 'Round, pad=0.015',
        linewidth = 3,
        edgecolor = 'red',
        facecolor = 'lightpink',
        alpha = 0.5
    ))
    fig.text(
        x = 0.5,
        y = 0.5,
        s = "Placeholder",
        ha = "center",
        va = "center",
        fontsize = 'xx-large',
        fontweight = 'bold',
        alpha = 0.5
    )

