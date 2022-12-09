# SPAdes
FROM virtool/workflow:5.2.1 as spades
WORKDIR /build
RUN wget https://github.com/ablab/spades/releases/download/v3.11.0/SPAdes-3.11.0-Linux.tar.gz
RUN tar -xvf SPAdes-3.11.0-Linux.tar.gz
RUN mv SPAdes-3.11.0-Linux spades
RUN sed -i 's/import collections/import collections\nimport collections.abc/g' spades/share/spades/pyyaml3/constructor.py
RUN sed -i 's/key, collections.Hashable/key, collections.abc.Hashable/g' spades/share/spades/pyyaml3/constructor.py

FROM virtool/workflow:5.2.1 as base
WORKDIR /workflow
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
COPY src src
COPY Cargo.lock Cargo.toml pyproject.toml poetry.lock ./
RUN poetry install
RUN poetry run maturin build
RUN poetry add target/wheels/nuvs_rust*.whl
COPY workflow.py /workflow/workflow.py
COPY --from=spades /build/spades /opt/spades
RUN ln -fs /opt/spades/bin/spades.py /usr/local/bin/spades.py

FROM base as test
COPY tests ./tests/
RUN poetry run pytest --cache-clear
