FROM mambaorg/micromamba:latest

ENV HOME=/home/mambauser

WORKDIR $HOME
COPY environment.yml .
RUN micromamba install --yes --name base -c conda-forge python=3.10 \
    && micromamba install --name base --yes -f environment.yml \
    && micromamba clean --all --yes
COPY ./etc ./etc/
COPY ./lib ./lib/
COPY Snakefile .

WORKDIR /local/
ENTRYPOINT ["micromamba", "run", "-n", "base", "snakemake", "-s", "/home/mambauser/Snakefile"]


