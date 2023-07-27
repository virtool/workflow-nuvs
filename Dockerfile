FROM debian:buster as prep
WORKDIR /build
RUN apt-get update && apt-get install -y make gcc zlib1g-dev wget unzip
RUN wget https://zlib.net/pigz/pigz-2.7.tar.gz && \
    tar -xzvf pigz-2.7.tar.gz && \
    cd pigz-2.7 && \
    make
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.2/bowtie2-2.3.2-legacy-linux-x86_64.zip && \
    unzip bowtie2-2.3.2-legacy-linux-x86_64.zip && \
    mkdir bowtie2 && \
    cp bowtie2-2.3.2-legacy/bowtie2* bowtie2

FROM python:3.10-buster as spades
WORKDIR /build
RUN wget https://github.com/ablab/spades/releases/download/v3.11.0/SPAdes-3.11.0-Linux.tar.gz
RUN tar -xvf SPAdes-3.11.0-Linux.tar.gz
RUN mv SPAdes-3.11.0-Linux spades
RUN sed -i 's/import collections/import collections\nimport collections.abc/g' spades/share/spades/pyyaml3/constructor.py
RUN sed -i 's/key, collections.Hashable/key, collections.abc.Hashable/g' spades/share/spades/pyyaml3/constructor.py

FROM python:3.10-buster as rust
WORKDIR /build
RUN apt-get update && apt-get install -y curl build-essential
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN pip install maturin==0.14.12
COPY src src
COPY Cargo.toml Cargo.lock pyproject.toml ./
RUN maturin build --release

FROM python:3.10-buster as pip
RUN pip install --upgrade pip
RUN pip install --user biopython
COPY --from=rust /build/target/wheels/nuvs_rust*.whl ./
RUN pip install --user nuvs_rust*.whl

FROM python:3.10-buster as base
WORKDIR /workflow
COPY --from=prep /build/bowtie2/* /usr/local/bin/
COPY --from=prep /build/FastQC /opt/fastqc
COPY --from=prep /build/pigz-2.7/pigz /usr/local/bin/pigz
COPY --from=spades /build/spades /opt/spades
COPY --from=pip /root/.local /root/.local
COPY poetry.lock pyproject.toml workflow.py ./
RUN curl -sSL https://install.python-poetry.org | python -
ENV PATH="/opt/spades/bin:/root/.local/bin:${PATH}"
RUN pip install --upgrade pip
RUN poetry export > requirements.txt
RUN pip install -r requirements.txt
RUN poetry export  --with dev > requirements.txt
RUN pip install -r requirements.txt

FROM base as test
COPY tests ./tests
RUN pytest