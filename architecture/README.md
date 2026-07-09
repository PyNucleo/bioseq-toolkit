# Architecture documentation

This directory documents the current repository architecture for manual audit. It is intentionally private-first and public-ready: the goal is to make current behavior traceable, not to market the project or describe planned features.

Scope rules used while writing these files:

- Only repository files were used as evidence: source code, tests, configuration, benchmark scripts, data paths, and existing docs.
- Production code, tests, dependency configuration, and the main `README.md` were not intentionally modified.
- Unclear, inconsistent, or untested behavior is recorded in `open_questions.md`.
- Diagrams use LikeC4 source files. No renderer dependency was added.

## Contents

- `model.c4` defines shared LikeC4 elements and relationships.
- `views/system-overview.c4` shows entry points, main modules, external services, tests, and benchmarks.
- `views/search-workflow.c4` traces single-query and multi-query k-mer search paths.
- `views/alignment-workflow.c4` separates Needleman-Wunsch and Smith-Waterman workflows.
- `views/fasta-and-database-flow.c4` traces FASTA parsing, writing, UniProt fetching, and `SequenceDatabase` normalization.
- `views/translation-and-benchmarks.c4` shows translation and benchmark relationships.
- `views/search-pipeline-detail.c4` drills into `search`, `multi_search`, regular multi-search, indexed multi-search, ranking, and refinement.
- `views/kmer-candidate-retrieval-detail.c4` drills into scan-based k-mer retrieval and indexed candidate lookup.
- `views/smith-waterman-detail.c4` drills into local alignment, score-only Smith-Waterman, refinement, score matrices, and traceback structures.
- `views/needleman-wunsch-detail.c4` drills into global alignment, score matrix filling, movement matrix creation, traceback, and structured result assembly.
- `views/fasta-database-normalization-detail.c4` drills into FASTA parsing/writing/fetching and database normalization.
- `views/translation-detail.c4` drills into FASTA-driven DNA filtering, transcription, translation, and translation result dictionaries.
- `audit-notes/*.md` contains module-by-module audit notes using a fixed evidence structure.
- `open_questions.md` lists unresolved behavior and verification risks.

## Suggested preview

If LikeC4 is installed locally, preview from the repository root:

```bash
npx likec4 start architecture
```

If `npx likec4` is not available, the `.c4` files are still plain source files and can be reviewed directly.

The LikeC4 source files were validated in this workspace with:

```bash
cmd /c likec4 validate architecture
```

## Verification performed

Repository structure, source modules, tests, benchmark scripts, benchmark reports, data paths, and package metadata were inspected. The full test command was attempted with:

```bash
python -m pytest
python -m pytest --basetemp .pytest-tmp -p no:cacheprovider
```

Both runs were blocked by Windows permission errors while pytest created or cleaned temporary directories. Before the temp-directory failure, 78 tests passed and the remaining 22 errors were setup errors tied to temp-path access, not project assertion failures. The second run left `.pytest-tmp` inaccessible to normal workspace commands.
