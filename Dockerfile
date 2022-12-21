FROM virtool/workflow:5.2.1 as spades
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
RUN pip install maturin==0.14.5
COPY src src
COPY Cargo.toml Cargo.lock ./
RUN maturin build --release
RUN mv target/wheels/nuvs_rust*.whl ./

FROM virtool/workflow:5.2.1 as base
WORKDIR /workflow
COPY workflow.py /workflow/workflow.py
RUN pip install --upgrade pip
COPY --from=rust /build/nuvs_rust*.whl .
RUN pip install biopython maturin nuvs_rust*.whl
COPY --from=spades /build/spades /opt/spades
ENV PATH="/opt/spades/bin:${PATH}"

FROM base as test
ADD tests ./tests
RUN pip install pytest pytest-asyncio==0.15.1 pytest-aiohttp==0.3.0 pytest-xdist pytest-regressions pydantic-factories
RUN pytest
