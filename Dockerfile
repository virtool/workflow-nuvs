FROM debian:bookworm as bowtie2
WORKDIR /build
RUN apt-get update && apt-get install -y build-essential cmake wget zlib1g-dev
RUN wget https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.5.4.tar.gz
RUN tar -xvf v2.5.4.tar.gz
WORKDIR bowtie2-2.5.4
RUN make
RUN mkdir /build/bowtie2
RUN cp bowtie2* /build/bowtie2/

FROM debian:bookworm as spades
WORKDIR /build
RUN apt-get update && apt-get install -y build-essential cmake libbz2-dev wget zlib1g-dev
RUN wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5.tar.gz
RUN tar -xvf SPAdes-3.15.5.tar.gz
WORKDIR SPAdes-3.15.5
ENV PREFIX=/build/spades
RUN ./spades_compile.sh

FROM python:3.12.3-bookworm as deps
WORKDIR /app
COPY --from=bowtie2 /build/bowtie2/* /usr/local/bin/
COPY --from=spades /build/spades /opt/spades
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/hmmer /opt/hmmer
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/skewer /usr/local/bin/
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/pigz /usr/local/bin/
RUN apt-get update && apt-get install -y --no-install-recommends curl build-essential default-jre

FROM deps as build
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:/opt/spades/bin:/opt/hmmer/bin/:${PATH}"
COPY uv.lock pyproject.toml ./
RUN uv sync --no-dev

FROM deps as test
ARG USER_ID
ARG GROUP_ID
RUN addgroup --gid $GROUP_ID appgroup
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID appuser
RUN chown -R appuser:appgroup /app
USER appuser
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/home/appuser/.local/bin:/opt/spades/bin:/opt/hmmer/bin/:${PATH}"
COPY uv.lock pyproject.toml ./
RUN uv sync

FROM python:3.12.3-bookworm as base
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
