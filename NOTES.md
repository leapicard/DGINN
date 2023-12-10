To export conda/mamba environment:

```shell
mamba env export > environment.yml
```

To run locally, first create and activate a conda/mamba environment from environment.yml, then:

```shell
snakemake --cores 1 --configfile config.yaml
```

To run with docker:

```shell
docker build . -t dginn
docker run --rm -u $(id -u $USER) -v $(pwd):/opt/local dginn --cores 1 --configfile config.yaml
```

To run with Singularity/Apptainer:

```shell
apptainer build dginn.sif Apptainer
apptainer run dginn.sif --cores 1 --configfile config.yaml
```

# TODO

- allow to start the pipeline from intermediary step by providing input file. What about `Data` ?
- the possel/aln example has two input files
- the `all` rule should only have the last step files as output ?
