FROM debian:buster as prep
WORKDIR /build
RUN apt-get update && apt-get install -y cmake gcc g++ make unzip wget zlib1g-dev
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.2/bowtie2-2.3.2-legacy-linux-x86_64.zip
RUN unzip bowtie2-2.3.2-legacy-linux-x86_64.zip
RUN mkdir bowtie2
RUN cp bowtie2-2.3.2-legacy/bowtie2* bowtie2
RUN wget https://github.com/ablab/spades/releases/download/v3.11.0/SPAdes-3.11.0-Linux.tar.gz
RUN tar -xvf SPAdes-3.11.0-Linux.tar.gz
RUN mv SPAdes-3.11.0-Linux spades
RUN sed -i 's/import collections/import collections\nimport collections.abc/g' spades/share/spades/pyyaml3/constructor.py
RUN sed -i 's/key, collections.Hashable/key, collections.abc.Hashable/g' spades/share/spades/pyyaml3/constructor.py

FROM python:3.10-bullseye as build
WORKDIR /app
COPY --from=prep /build/bowtie2/* /usr/local/bin/
COPY --from=prep /build/spades /opt/spades
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/hmmer /opt/hmmer
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/skewer /usr/local/bin/
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/pigz /usr/local/bin/
RUN apt-get update && apt-get install -y --no-install-recommends curl build-essential default-jre
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -sSL https://install.python-poetry.org | python -
ENV PATH="/root/.cargo/bin:/root/.local/bin:/opt/spades/bin:/opt/hmmer/bin/:${PATH}" \
    POETRY_CACHE_DIR='/tmp/poetry_cache' \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1
COPY src src
COPY Cargo.toml Cargo.lock poetry.lock pyproject.toml ./
RUN poetry install --only rust
RUN poetry run maturin build --release
RUN poetry remove nuvs-rust && poetry add target/wheels/*.whl && poetry install
RUN poetry install --without dev --no-root && rm -rf $POETRY_CACHE_DIR

FROM build as test
ENV PATH="/root/.cargo/bin:/root/.local/bin:/opt/spades/bin:${PATH}" \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1
RUN poetry install --with dev
COPY example ./example
COPY tests ./tests
COPY workflow.py ./
ENTRYPOINT ["poetry", "run"]

FROM python:3.10-bullseye as base
WORKDIR /app
COPY --from=prep /build/bowtie2/* /usr/local/bin/
COPY --from=prep /build/spades /opt/spades
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/hmmer /opt/hmmer
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/skewer /usr/local/bin/
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/pigz /usr/local/bin/
ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:/opt/spades/bin:/opt/hmmer/bin:${PATH}"
COPY --from=build /app/.venv /app/.venv
COPY workflow.py VERSION* ./
