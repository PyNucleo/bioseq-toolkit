# Architecture documentation

This directory describes the current repository architecture. Source code and
tests remain authoritative; the model, views, and narratives are synchronized
descriptions rather than a separate specification.

## Contents

- `model.c4` defines shared LikeC4 elements, contracts, and relationships.
- `views/system-overview.c4` shows entry points, modules, externals, tests, and benchmarks.
- `views/search-workflow.c4` and the search detail views show regular/indexed candidate retrieval and optional refinement.
- `views/alignment-workflow.c4` and the alignment detail views show the linear-gap Needleman-Wunsch and Smith-Waterman workflows.
- `views/fasta-and-database-flow.c4` and `views/fasta-database-normalization-detail.c4` show the shared parser, UniProt fetch, writer, and normalization boundaries.
- `views/translation-and-benchmarks.c4` and `views/translation-detail.c4` show accountable translation outcomes and benchmark consumers.
- `audit-notes/*.md` documents current module contracts and evidence.
- `open_questions.md` contains only unresolved or partially resolved behavior.

## Preview and validation

From the repository root, the documented LikeC4 workflow is:

```bash
npx likec4 start architecture
```

The local validation command used by this repository is:

```bash
likec4 validate architecture
```

At this synchronization revision, that command reported `Valid (12 files)`,
covering the model and every view.

LikeC4 is documentation tooling, not a Bioseq runtime dependency. If the tool
is unavailable, static element-reference checks are the minimum fallback and
the missing render/parse validation must be reported.

## Verification baseline

At documentation synchronization base revision `3ce3f74`, `python -m pytest -q`
completed with `199 passed`. That count is revision-scoped. The architecture
documents reflect current source signatures and tests at that revision.
