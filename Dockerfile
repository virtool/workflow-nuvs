# SPAdes
FROM alpine:3.14 as spades
WORKDIR /build
RUN wget https://github.com/ablab/spades/releases/download/v3.11.0/SPAdes-3.11.0-Linux.tar.gz && \
    tar -xvf SPAdes-3.11.0-Linux.tar.gz && \
    mv SPAdes-3.11.0-Linux spades

FROM virtool/workflow:4.2.2 as build
WORKDIR /workflow
COPY --from=spades /build/spades /opt/spades
RUN ln -fs /opt/spades/bin/spades.py /usr/local/bin/spades.py
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
COPY src src
COPY Cargo.toml Cargo.toml
COPY pyproject.toml pyproject.toml
COPY poetry.lock poetry.lock
RUN poetry install
RUN poetry run maturin build
RUN poetry add target/wheels/nuvs_rust*.whl
COPY workflow.py /workflow/workflow.py
