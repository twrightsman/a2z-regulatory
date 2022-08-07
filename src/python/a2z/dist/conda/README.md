# Conda Build Instructions

## Install conda-build

```
$ conda install conda-build
```

## Build the conda package

```
git clone https://github.com/twrightsman/a2z-regulatory
conda-build -c bioconda -c conda-forge a2z-regulatory/src/python/a2z/dist/conda/a2z-regulatory
```

## Use the package

```
conda create -c local -c bioconda -c conda-forge -n a2z-regulatory a2z-regulatory
conda activate a2z-regulatory
```

## Troubleshooting

- Remove the `$CONDA_PREFIX/conda-bld` directory
- Remove snakemake conda environments in `.snakemake` in working directory and re-run

