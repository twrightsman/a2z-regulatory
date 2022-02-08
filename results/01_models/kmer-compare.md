---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
import hashlib
import itertools
import json
import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyfaidx import Fasta
import sklearn.feature_extraction.text
import sklearn.linear_model
import tensorflow as tf

from a2z.data import IntervalReferenceDataset
from a2z.models import allow_growth
from a2z.plots import roc_curve, pr_curve

allow_growth()
```

```{code-cell} ipython3
def createKmerSet(kmersize):
    '''
    write all possible kmers
    :param kmersize: integer, 8
    :return uniq_kmers: list of sorted unique kmers
    '''
    kmerSet = set()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = itertools.product(nucleotides, repeat=kmersize)
    for i in kmerall:
        kmer = ''.join(i)
        kmerSet.add(kmer)
    uniq_kmers = sorted(list(kmerSet))  
    return uniq_kmers


def compute_kmer_entropy(kmer):
    '''
    compute shannon entropy for each kmer
    :param kmer: string
    :return entropy: float
    '''
    prob = [float(kmer.count(c)) / len(kmer) for c in dict.fromkeys(list(kmer))]
    entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob ])
    return round(entropy, 2)


def make_stopwords(kmersize):
    '''
    write filtered out kmers
    :param kmersize: integer, 8
    :return stopwords: list of sorted low-complexity kmers
    '''
    kmersize_filter = {5:1.3, 6:1.3, 7:1.3, 8:1.3, 9:1.3, 10:1.3}
    limit_entropy = kmersize_filter.get(kmersize)
    kmerSet = set()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = itertools.product(nucleotides, repeat=kmersize)
    for n in kmerall:
        kmer = ''.join(n)
        
        if compute_kmer_entropy(kmer) < limit_entropy:
            kmerSet.add(make_newtoken(kmer))
        else:
            continue
    stopwords = sorted(list(kmerSet))
    return stopwords

  
def createNewtokenSet(kmersize):
    '''
    write all possible newtokens
    :param kmersize: integer, 8
    :return uniq_newtokens: list of sorted unique newtokens
    ''' 
    newtokenSet = set()
    uniq_kmers = createKmerSet(kmersize)
    for kmer in uniq_kmers:
        newtoken = make_newtoken(kmer)
        newtokenSet.add(newtoken)  
    uniq_newtokens = sorted(list(newtokenSet))
    return uniq_newtokens      


def make_newtoken(kmer):
    '''
    write a collapsed kmer and kmer reverse complementary as a newtoken
    :param kmer: string e.g., "AT"
    :return newtoken: string e.g., "atnta"
    :param kmer: string e.g., "TA"
    :return newtoken: string e.g., "atnta"
    '''
    kmer = str(kmer).lower()
    newtoken = "n".join(sorted([kmer,kmer.translate(str.maketrans('tagc', 'atcg'))[::-1]]))
    return newtoken

def write_ngrams(sequence, kmerlength = 8):
    '''
    write a bag of newtokens of size n
    :param sequence: string e.g., "ATCG"
    :param (intern) kmerlength e.g., 2
    :return newtoken_string: string e.g., "atnta" "gatc" "cgcg" 
    '''
    seq = str(sequence).lower()
    finalstart = (len(seq)-kmerlength)+1
    allkmers = [seq[start:(start+kmerlength)] for start in range(0,finalstart)]
    tokens = [make_newtoken(kmer) for kmer in allkmers if len(kmer) == kmerlength and "n" not in kmer]
    newtoken_string = " ".join(tokens)
    return newtoken_string
```

```{code-cell} ipython3
def plot_bok_vs_a2z(configuration: str, ax: matplotlib.axes.Axes):
    with open(f"configs/{configuration}.json") as config_file:
        config = json.load(config_file)

    split_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n" + json.dumps(config['splitting'], indent = 2) + "\n").encode('utf-8')).hexdigest()
    config_hash = hashlib.sha256((json.dumps(config, indent = 2) + "\n").encode('utf-8')).hexdigest()
    split_data_path = f"tmp/data/split/{split_hash}.tsv"

    intervals = pd.read_table(
        split_data_path,
        header = 0,
        dtype = {
            'seqid': str,
            'start': int,
            'end': int,
            'species': str,
            'target': np.byte,
            'is_train': np.byte,
            'is_test': np.byte,
            'is_val': np.byte
        }
    )
    intervals = intervals.loc[intervals['species'] == 'Zea mays']

    reference = Fasta(f"{config['paths']['genomes']}/Zea_mays.fa", as_raw = True, rebuild = False)

    intervals['sequence'] = intervals.apply(lambda row: reference[row.seqid][row.start:row.end], axis = 1)

    kmer_length = 8
    filtered = True

    all_tokens = createNewtokenSet(kmer_length)

    if kmer_length > 4:
        stpwrds = make_stopwords(kmer_length)
    else:
        stpwrds = []

    intervals['tokens'] = intervals['sequence'].map(write_ngrams)

    filtered_kmers = list(filter(lambda k: k in stpwrds, all_tokens)) if filtered else all_tokens

    vectorizer = sklearn.feature_extraction.text.TfidfVectorizer(sublinear_tf = True, vocabulary = filtered_kmers)
    X_TFIDF_DEV = vectorizer.fit_transform(intervals.loc[intervals['is_train'] == 1, 'tokens'])

    Y_DEV = np.asarray(intervals.loc[intervals['is_train'] == 1, 'target'])
    Y_holdout = np.asarray(intervals.loc[intervals['is_test'] == 1, 'target'])
    X_TFIDF_test = vectorizer.fit_transform(intervals.loc[intervals['is_test'] == 1, 'tokens'])

    TFIDF_LR = sklearn.linear_model.LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
               intercept_scaling=1, max_iter=100, multi_class='ovr', n_jobs=1,
               penalty='l2', random_state=None, solver='liblinear', tol=0.0001,
               verbose=0, warm_start=False)
    TFIDF_LR.fit(X_TFIDF_DEV, Y_DEV)

    LR_hold_TFIDF_pred = TFIDF_LR.predict(X_TFIDF_test) # y_pred
    LR_hold_TFIDF_prob = TFIDF_LR.predict_proba(X_TFIDF_test)[:,1] # y_score

    test_species = "Zea mays"
    model_path = f"tmp/results/{configuration}/{test_species.replace(' ', '_')}/0/model"
    model = tf.keras.models.load_model(model_path)

    data = IntervalReferenceDataset(intervals, references = {'Zea mays': reference})

    intervals['a2z_prediction'] = model.predict(data).flatten()

    bok_precision, bok_recall, bok_thresholds = sklearn.metrics.precision_recall_curve(Y_holdout, LR_hold_TFIDF_prob)
    bok_auPR = sklearn.metrics.average_precision_score(Y_holdout, LR_hold_TFIDF_prob)
    a2z_precision, a2z_recall, a2z_thresholds = sklearn.metrics.precision_recall_curve(intervals.loc[intervals['is_test'] == 1, 'target'], intervals.loc[intervals['is_test'] == 1, 'a2z_prediction'])
    a2z_auPR = sklearn.metrics.average_precision_score(intervals.loc[intervals['is_test'] == 1, 'target'], intervals.loc[intervals['is_test'] == 1, 'a2z_prediction'])
    ax.plot(bok_recall, bok_precision, label = f"bag-of-kmers (auPR = {bok_auPR:.2f})")
    ax.plot(a2z_recall, a2z_precision, label = f"a2z (auPR = {a2z_auPR:.2f})")
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend()
```

```{code-cell} ipython3
fig, axs = plt.subplots(constrained_layout = True, nrows = 1, ncols = 2, figsize = (12, 5))

%time plot_bok_vs_a2z("accessibility_base_within", axs[0])
axs[0].set_title("Accessibility")
%time plot_bok_vs_a2z("hypomethylation_base_within", axs[1])
axs[1].set_title("Methylation")

fig.savefig("figs/bok-vs-a2z.png", facecolor = "white", dpi = 300)
```

```{code-cell} ipython3

```
