---
name: Release preparation
about: Prepare a new release for this package
title: 'Release : '
labels: 'release-prep'
assignees: ''

---

To create a new release, please do the following:

- [ ] make any changes to the repo needed to accommodate the release
- [ ] run `rever` to generate the `CHANGELOG` entry from `news` items
- [ ] merge this PR into `main`


After you have merged this PR, please also:

- [ ]  create a release here in GitHub, including generating detailed release notes
    - [ ] in the summary release notes, also link to the `CHANGELOG` entry for this release in the `gufe` docs
- [ ]  await automated PR on `conda-forge` `gufe-feedstock`
    - [ ] review and merge
