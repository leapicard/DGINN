FROM mambaorg/micromamba:latest


WORKDIR /opt
COPY environment.yml .
RUN micromamba install --yes --name base -c conda-forge python=3.10 \
    && micromamba install --name base --yes -f environment.yml \
    && micromamba clean --all --yes
COPY etc/ .
COPY lib/ .
COPY DGINN.py .

ENV HOME=/home/mambauser
RUN chmod 777 $HOME

WORKDIR /opt/local/
ENTRYPOINT ["micromamba", "run", "-n", "base", "python", "/opt/DGINN.py"]
CMD ["-h"]

