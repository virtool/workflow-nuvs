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

### Commits

All commits must follow the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0) specification.

These standardized commit messages are used to automatically publish releases using [`semantic-release`](https://semantic-release.gitbook.io/semantic-release)
after commits are merged to `main` from successful PRs.

**Example**

```text
feat: add API support for assigning labels to existing samples
```

Descriptive bodies and footers are required where necessary to describe the impact of the commit. Use bullets where appropriate.

Additional Requirements
1. **Write in the imperative**. For example, _"fix bug"_, not _"fixed bug"_ or _"fixes bug"_.
2. **Don't refer to issues or code reviews**. For example, don't write something like this: _"make style changes requested in review"_.
Instead, _"update styles to improve accessibility"_.
3. **Commits are not your personal journal**. For example, don't write something like this: _"got server running again"_
or _"oops. fixed my code smell"_.

From Tim Pope: [A Note About Git Commit Messages](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
