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

#### Setup (one-time)

Set environment variables for proper file permissions:
```sh
export USER_ID=$(id -u) GROUP_ID=$(id -g)
```

Build the test container:
```sh
docker compose build
```

#### Running Tests

Run all tests:
```sh
docker compose run --rm app poetry run pytest
```

Run specific tests:
```sh
docker compose run --rm app poetry run pytest tests/test_workflow.py
```

#### Why This Approach?

- **No rebuilds**: Uses bind mounts to get latest code without rebuilding
- **No permission issues**: Container runs as your user via USER_ID/GROUP_ID
- **Fast iterations**: Only Poetry and pytest execution, no Docker build time

### Commits

Read [our guide](https://dev.virtool.ca/en/latest/commits_releases.html#commits) on
writing commits for Virtool repositories.
