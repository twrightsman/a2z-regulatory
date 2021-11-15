# Interpretation experiments on a2z

## TF-MoDISco

```
mkdir -p tmp/tf-modisco
./run_modisco.py accessibility_base "Zea mays"
./modisco_to_meme.py accessibility_base "Zea mays"
tomtom -oc tmp/tf-modisco/accessibility_base/Zea_mays.tomtom tmp/tf-modisco_accessibility_base_Zea_mays.meme ../../data/JASPAR/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt
```

### Importance scores

Example TF-MoDISco HDF5 output:

```
root
├─ contrib_scores
│   ├─ task0 (shape: (n_seqs, max_seq_length, 4))
│   └─ ...
└─ hyp_contrib_scores
    ├─ task0 (shape: (n_seqs, max_seq_length, 4))
    └─ ...
```

`contrib_scores` would be $gradient \times input$ while `hyp_contrib_scores` would be (mean-normalized) $gradient$ where each position sums to zero.

Sequences do *not* need to be all the same length for TF-MoDISco.

For discussion on why mean-normalized gradients, see [tfmodisco#5](https://github.com/kundajelab/tfmodisco/issues/5).

## kmer Occlusion

```
simple_gpu_scheduler --gpus 0 < jobs_occlusion > tmp/jobs.out 2> tmp/jobs.err
./generate_null_profile.py accessibility_base "Arabidopsis thaliana" > tmp/kmer-occlusion/accessibility_base/
Arabidopsis_thaliana/nulls.tsv 2> tmp/kmer-occlusion/accessibility_base/Arabidopsis_thaliana/nulls.err
./generate_null_profile.py accessibility_base "Zea mays" > tmp/kmer-occlusion/accessibility_base/Zea_mays/nulls.tsv 2> tmp/kmer-occlusion/accessibility_base/Zea_mays/nulls.err
./generate_null_profile.py hypomethylation_base "Arabidopsis thaliana" > tmp/kmer-occlusion/hypomethylation_base/Arabidopsis_thaliana/nulls.tsv 2> tmp/kmer-occlusion/hypomethylation_base/Arabidopsis_thaliana/nulls.err
./generate_null_profile.py hypomethylation_base "Zea mays" > tmp/kmer-occlusion/hypomethylation_base/Zea_mays/nulls.tsv 2> tmp/kmer-occlusion/hypomethylation_base/Zea_mays/nulls.err
simple_gpu_scheduler --gpus 0,1 < jobs_pGIA > tmp/jobs_pGIA.out 2> tmp/jobs_pGIA.err
parallel -j6 ./fimo_kmers.sh {} ../../data/JASPAR/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt ::: tmp/kmer-occlusion/*/*/10
```

You will then need to run the following notebooks in order:

1. `jaspar_matched_kmers.ipynb`
2. `JASPAR_pGIA.ipynb`
3. `mds.ipynb`
4. `fig3.ipynb`

