## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
- [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/wslh-bio/spriggan/tree/main/.github/CONTRIBUTING.md)
- [ ] If necessary, also make a PR on the wslh-bio/spriggan _branch_ on the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository.
- [ ] Ensure the test suite passes (`nextflow run . -profile test,docker --outdir <OUTDIR>`).
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
