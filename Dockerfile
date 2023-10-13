FROM debian:buster as prep
WORKDIR /build
RUN apt-get update && apt-get install -y cmake gcc g++ make unzip wget zlib1g-dev
RUN wget https://zlib.net/pigz/pigz-2.8.tar.gz
RUN tar -xvf pigz-2.8.tar.gz
WORKDIR /build/pigz-2.8
RUN make
WORKDIR /build
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.2/bowtie2-2.3.2-legacy-linux-x86_64.zip
RUN unzip bowtie2-2.3.2-legacy-linux-x86_64.zip
RUN mkdir bowtie2
RUN cp bowtie2-2.3.2-legacy/bowtie2* bowtie2
RUN wget https://github.com/ablab/spades/releases/download/v3.11.0/SPAdes-3.11.0-Linux.tar.gz
RUN tar -xvf SPAdes-3.11.0-Linux.tar.gz
RUN mv SPAdes-3.11.0-Linux spades
RUN sed -i 's/import collections/import collections\nimport collections.abc/g' spades/share/spades/pyyaml3/constructor.py
RUN sed -i 's/key, collections.Hashable/key, collections.abc.Hashable/g' spades/share/spades/pyyaml3/constructor.py
RUN wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz
RUN tar -xf hmmer-3.2.1.tar.gz
WORKDIR /build/hmmer-3.2.1
RUN ./configure --prefix /build/hmmer
RUN make
RUN make install
WORKDIR /build
RUN wget https://github.com/relipmoc/skewer/archive/0.2.2.tar.gz
RUN tar -xf 0.2.2.tar.gz
WORKDIR /build/skewer-0.2.2
RUN make
RUN mv skewer /build

FROM python:3.10-buster as base
WORKDIR /app
COPY --from=prep /build/bowtie2/* /usr/local/bin/
COPY --from=prep /build/FastQC /opt/fastqc
COPY --from=prep /build/hmmer /opt/hmmer
COPY --from=prep /build/pigz-2.8/pigz /usr/local/bin/pigz
COPY --from=prep /build/skewer /usr/local/bin/
COPY --from=prep /build/spades /opt/spades
RUN chmod ugo+x /opt/fastqc/fastqc && \
    ln -fs /opt/fastqc/fastqc /usr/local/bin/fastqc && \
    for file in `ls /opt/hmmer/bin`; do ln -fs /opt/hmmer/bin/${file} /usr/local/bin/${file};  done
RUN apt-get update && \
    apt-get install -y --no-install-recommends curl build-essential default-jre && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean
RUN apt-get update && apt-get install -y curl build-essential
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -sSL https://install.python-poetry.org | python -
ENV PATH="/root/.cargo/bin:/root/.local/bin:/opt/spades/bin:${PATH}"
RUN pip install --upgrade pip
RUN pip install maturin==0.14.12
COPY src src
COPY Cargo.toml Cargo.lock poetry.lock pyproject.toml workflow.py ./
RUN maturin build --release
RUN poetry export > requirements.txt
RUN pip install -r requirements.txt
RUN pip install /app/target/wheels/nuvs_rust*.whl

FROM base as test
WORKDIR /app
ENV PATH="/opt/spades/bin:${PATH}"
RUN poetry export --with dev > requirements.txt
RUN pip install -r requirements.txt
COPY workflow.py ./
COPY tests ./tests
RUN pytest
