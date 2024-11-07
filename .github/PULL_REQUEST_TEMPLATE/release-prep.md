---
name: Release preparation
about: Prepare a new release for this package
title: 'Release : '
labels: 'release-prep'
assignees: ''

---

To create a new release, please do the following:
- [ ] make any changes to the repo needed to accommodate the release
    - [ ] check `environment.yml` for any updates needed; create any issues needed to update versions in next release cycle
- [ ] run [`rever`](https://regro.github.io/rever-docs/index.html#rever-releaser-of-versions) to generate the `CHANGELOG` entry from `news` items: `rever <version_number>`
- [ ] merge this PR into `main`


After you have merged this PR, please also:
- [ ]  create a release here in GitHub, including generating detailed release notes
    - [ ] in the summary release notes, also link to the `CHANGELOG` entry for this release in the `gufe` docs
- [ ]  await automated PR on `conda-forge` [`gufe-feedstock` `meta.yaml`](https://github.com/conda-forge/gufe-feedstock/blob/main/recipe/meta.yaml)
    - [ ] review and merge
