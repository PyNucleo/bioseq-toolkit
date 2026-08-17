# UniProt fetch demonstration

This case study records the observed behavior of a dated fetch demonstration.
UniProt service state and accession histories can change, so the recorded counts
are not a promise about current remote results.

The command reads a plain-text file containing one accession per nonempty line,
writes successfully parsed records to FASTA, and prints a human-readable
summary:

```bash
bioseq fetch-uniprot -f accessions.txt -o sequences.fasta
```

Useful flags are:

- `--strict`: raise on the first expected HTTP, requests-layer/network, or empty
  HTTP-200 response failure instead of collecting it and continuing;
- `--full-header`: write each stored FASTA header; otherwise select a short
  header independently per record, preferring accession then ID; and
- `--show-failed`: print details for failures collected in non-strict mode.

Successful nonempty response bodies use the same strict structural FASTA parser
as local files. A malformed nonempty body can therefore raise parser
`ValueError` even in non-strict mode; non-strict mode does not convert every
possible malformed body into a `failed` entry.

## Historical observation

In one dated run over the included 39-line accession set, 33 records were
written, 6 accessions were reported as operational failures, and the output
contained 56,497 amino-acid residues. Repeated elapsed times were approximately
40-64 seconds. These values are historical observations, not a reproducible
algorithmic benchmark or a statement about current accession status.

The earlier version of this note asserted a specific demerge/secondary-
accession history for `Q9V3G5`. That provenance was not reliably established,
so the claim and the purported replacement accessions have been removed rather
than replaced with another unverified history. The fetcher deliberately does
not guess replacement records for unresolved accessions.
