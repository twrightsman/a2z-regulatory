import numpy as np
import tensorflow as tf


def generate_null_profile(seqs: np.ndarray, n: int, model: tf.keras.models.Model) -> (np.ndarray, np.ndarray):
    """ generate null sequences using a profile model and run predictions """
    seq_profile = np.mean(seqs, axis = 0)
    seq_profile /= seq_profile.sum(axis = 1, keepdims = True)

    null_set_profile = np.zeros(shape = (n, seqs.shape[1], 4))
    for i, seq in enumerate(null_set_profile):
        for pos in range(seqs.shape[1]):
            seq[pos][np.random.choice(4, p = seq_profile[pos])] = 1

    return null_set_profile, model.predict(null_set_profile).flatten()


def global_importance_analysis(motif: np.ndarray, embed_location: int, model: tf.keras.models.Model, null_set: np.ndarray, null_set_preds: np.ndarray) -> float:
    embedded_null_set = np.copy(null_set)
    for seq in embedded_null_set:
        seq[embed_location:embed_location+len(motif)] = motif
    embedded_null_set_preds = model.predict(embedded_null_set).flatten()
    global_importance = (embedded_null_set_preds - null_set_preds).mean()
    return global_importance

