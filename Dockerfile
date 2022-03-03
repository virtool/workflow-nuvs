# SPAdes
FROM alpine:3.14 as spades
WORKDIR /build
RUN wget https://github.com/ablab/spades/releases/download/v3.11.0/SPAdes-3.11.0-Linux.tar.gz && \
    tar -xvf SPAdes-3.11.0-Linux.tar.gz && \
    mv SPAdes-3.11.0-Linux spades

FROM virtool/workflow:2.1.2 as build
WORKDIR /workflow
COPY --from=spades /build/spades /opt/spades
RUN ln -fs /opt/spades/bin/spades.py /usr/local/bin/spades.py
COPY workflow.py /workflow/workflow.py
