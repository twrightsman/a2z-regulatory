# Conda Build Instructions

## Install conda-build

```
$ conda install conda-build
```

## Build the conda package (required for peak-pipeline)

```
git clone https://bitbucket.org/bucklerlab/a2z-regulatory
conda-build -c conda-forge a2z-regulatory/src/python/a2z/dist/conda/portion
conda-build -c local -c bioconda -c conda-forge a2z-regulatory/src/python/a2z/dist/conda/a2z-regulatory
```

## Use the package (not required for peak-pipeline)

```
conda create -c local -c bioconda -c conda-forge -n a2z-regulatory a2z-regulatory
conda activate a2z-regulatory
```

## Troubleshooting

- Clear your `$CONDA_PREFIX/conda-bld` directory (remove it)
- Remove snakemake conda environments in `.snakemake` in working directory and re-run

