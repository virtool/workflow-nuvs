FROM debian:bookworm AS spades
WORKDIR /build
RUN apt-get update && apt-get install -y build-essential cmake libbz2-dev wget zlib1g-dev
RUN wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5.tar.gz
RUN tar -xvf SPAdes-3.15.5.tar.gz
WORKDIR SPAdes-3.15.5
ENV PREFIX=/build/spades
RUN ./spades_compile.sh

FROM python:3.13-bookworm AS deps
WORKDIR /app
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/bowtie2/2.5.4/bowtie* /usr/local/bin/
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/pigz/2.8/pigz /usr/local/bin/
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/hmmer/3.3.2/ /opt/hmmer
COPY --from=spades /build/spades /opt/spades

FROM python:3.13-bookworm AS uv
WORKDIR /app
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.cargo/bin:/root/.local/bin:${PATH}" \
    UV_CACHE_DIR='/tmp/uv_cache'
COPY uv.lock pyproject.toml README.md ./
RUN uv sync

FROM deps AS dev
WORKDIR /app
ENV PATH="/root/.cargo/bin:/root/.local/bin:${PATH}" \
    UV_CACHE_DIR='/tmp/uv_cache'
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -LsSf https://astral.sh/uv/install.sh | sh

FROM dev AS test
COPY uv.lock pyproject.toml ./
COPY README.md ./
RUN uv sync
COPY example ./example
COPY tests ./tests
COPY utils.py workflow.py VERSION* ./

FROM deps AS base
WORKDIR /app
ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:/opt/spades/bin:/opt/hmmer/bin:${PATH}"
COPY --from=uv /app/.venv /app/.venv
COPY utils.py workflow.py VERSION* ./
