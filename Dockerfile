FROM debian:bookworm as bowtie2
WORKDIR /build
RUN apt-get update && apt-get install -y unzip wget
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip
RUN unzip bowtie2-2.5.4-linux-x86_64.zip
RUN mkdir bowtie2
RUN cp bowtie2-2.5.4-linux-x86_64/bowtie2* bowtie2

FROM debian:bookworm as spades
WORKDIR /build
RUN apt-get update && apt-get install -y build-essential cmake libbz2-dev wget zlib1g-dev
RUN wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5.tar.gz
RUN tar -xvf SPAdes-3.15.5.tar.gz
WORKDIR SPAdes-3.15.5
ENV PREFIX=/build/spades
RUN ./spades_compile.sh

FROM python:3.12-bookworm as build
WORKDIR /app
COPY --from=bowtie2 /build/bowtie2/* /usr/local/bin/
COPY --from=spades /build/spades /opt/spades
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/hmmer /opt/hmmer
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/skewer /usr/local/bin/
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/pigz /usr/local/bin/
RUN apt-get update && apt-get install -y --no-install-recommends curl build-essential default-jre
RUN curl -sSL https://install.python-poetry.org | python -
ENV PATH="/root/.local/bin:/opt/spades/bin:/opt/hmmer/bin/:${PATH}" \
    POETRY_CACHE_DIR='/tmp/poetry_cache' \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1
COPY poetry.lock pyproject.toml ./
RUN poetry install --without dev --no-root && rm -rf $POETRY_CACHE_DIR

FROM build as test
ENV PATH="/root/.local/bin:/opt/spades/bin:${PATH}" \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1
RUN poetry install --with dev
COPY example ./example
COPY tests ./tests
COPY utils.py workflow.py ./
ENTRYPOINT ["poetry", "run"]

FROM python:3.12-bookworm as base
WORKDIR /app
COPY --from=bowtie2 /build/bowtie2/* /usr/local/bin/
COPY --from=spades /build/spades /opt/spades
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/hmmer /opt/hmmer
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/skewer /usr/local/bin/
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/pigz /usr/local/bin/
ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:/opt/spades/bin:/opt/hmmer/bin:${PATH}"
COPY --from=build /app/.venv /app/.venv
COPY utils.py workflow.py VERSION* ./
