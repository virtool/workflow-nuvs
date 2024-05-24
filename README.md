workflow-nuvs
-------------

A workflow for identifying novel viruses in Virtool.

## Steps

1. Eliminate sample reads that map to any default (representative) isolate of any OTU.
2. Eliminate sample reads that map to the configured subtraction.
3. Repair paired reads if some pair members were lost in elimination.
4. Assemble the remaining reads using [SPAdes](https://github.com/ablab/spades).
5. Calculate ORFs from the assembled contigs.
6. Use ORFs as input for [HMMER](http://hmmer.org/) to detect viral motifs using profile hidden Markov models derived from the vFAM project.

## Contributing

### Running Tests

To run the tests for this package, use the following command:

```sh
docker run $(docker build -q .) pytest
```

If you want to run specific tests, you can specify the test file or directory as an
argument to `pytest`:

```sh
docker run $(docker build -q .) pytest tests/test_workflow.py
```

If you get an error like:
```
Unable to find image 'pytest:latest' locally
```

Your build is likely failing. Run the build separately to make sure it works:
```sh
docker build .
```

### Commits

Read [our guide](https://dev.virtool.ca/en/latest/commits_releases.html#commits) on
writing commits for Virtool repositories.
