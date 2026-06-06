````markdown
## UniProt Fetching Notes

The `fetch-uniprot` command can retrieve protein sequences from a text file of UniProt accessions and write the successfully fetched records to a FASTA file.

Example workflow:

```bash
bioseq fetch-uniprot -f accessions.txt -o sequences.fasta
````

Each non-empty line in the input file is treated as one UniProt accession. Successfully fetched records are parsed into FASTA format and written to the output file. Accessions that do not return valid FASTA records are reported separately instead of being silently written or replaced.

### Accession validity

Not every valid-looking UniProt accession maps directly to a current, unique FASTA sequence. Some accessions may be secondary, merged, demerged, deleted, or otherwise affected by UniProt entry-history changes.

For example, an accession may exist historically but no longer represent one unique current protein entry. In that case, UniProt may not return a normal FASTA record for the accession. The fetch command treats responses as valid only when they contain FASTA text beginning with a `>` header line. One such observed example is the accession Q9V3G5, which is an entry that has been demerged. Its accession has been set as a secondary accession in UniProt for P8262 and P8263, and such, was not included in the fetched results file and instead was added into failed accessions with a "reason" key.

This behavior avoids silently choosing a replacement sequence when an accession is ambiguous.

### Runtime note

In one test accession file containing 39 UniProt accessions:

* 33 accessions returned valid FASTA records
* 6 accessions failed or returned invalid/non-FASTA responses
* 56,497 total amino-acid residues were written
* repeated runs took approximately 40–64 seconds

This runtime should not be interpreted as a pure algorithmic benchmark. The command depends on remote UniProt requests, so elapsed time is affected by network latency, UniProt server response time, accession history behavior, and the number of requested accessions. For datasets of this size, local parsing and FASTA writing are expected to be minor contributors compared with remote fetching.

