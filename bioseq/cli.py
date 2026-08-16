import argparse
import json
import textwrap
from bioseq.pipelines.search_pipeline import search, multi_search
from bioseq.alignment.smith_waterman import local_alignment
from bioseq.alignment.needleman_wunsch import global_alignment
from bioseq.fasta_io import fetch_uniprot_sequences, write_fasta_records

def format_failed_accessions(failed_accessions):
    final_text = ""

    for num, failed_record in enumerate(failed_accessions, start=1):
        reason_string = textwrap.fill(
    failed_record["reason"],
    width=88,
    initial_indent="   Reason: ",
    subsequent_indent="           ",
)
        final_text += (f'{num}. {failed_record["accession"]}\n'
                       f'   Status code: {failed_record["status_code"]}\n'
                       f'{reason_string}\n\n'
        )
    
    return final_text

def main():
    
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest="command", required=True)

    search_parser = subparsers.add_parser("search", help = "Run kmer search")
    local_alignment_parser = subparsers.add_parser("align-local", help = "Run local alignment on two sequences")
    global_alignment_parser = subparsers.add_parser("align-global", help = "Run global alignment on two sequences")
    fetch_uniprot_parser = subparsers.add_parser("fetch-uniprot", help = "Fetch UniProt sequences through their accessions.", description="Takes the path for a FASTA file containing valid ids, returns a list of FASTA sequences with their header's details.")
    multi_query_search_parser = subparsers.add_parser(
        "multi-search",
        help="Perform a search of top hits for query sequences inside a FASTA file.",
        description=(
            "Search multiple queries in indexed mode by default or regular mode with "
            "--regular. --refine is supported in both modes; refinement scoring "
            "options are inactive when refinement is disabled."
        ),
    )

    search_parser.add_argument('-q', '--query', type=str, required=True)
    search_parser.add_argument('-d', '--database', type=str, required=True)
    search_parser.add_argument('-k', '--kmer-size', type=int, default=3)
    search_parser.add_argument('-t', '--threshold',type=int, default=3)
    search_parser.add_argument('-n', '--top-n-hits', type=int, default=10)
    search_parser.add_argument('-r', '--refine', action='store_true')
    search_parser.add_argument('-m', '--match', type=int, default=1)
    search_parser.add_argument('-u', '--mismatch', type=int, default=-1)
    search_parser.add_argument('-g', '--gap-penalty', type=int, default=-2)
    search_parser.add_argument('-x', '--matrix', type=str, default=None)

    local_alignment_parser.add_argument('-s1', '--sequence1', type=str, required=True)
    local_alignment_parser.add_argument('-s2', '--sequence2', type=str, required=True)
    local_alignment_parser.add_argument('-m', '--match', type=int, default=1)
    local_alignment_parser.add_argument('-u', '--mismatch', type=int, default=-1)
    local_alignment_parser.add_argument('-g', '--gap-penalty', type=int, default=-2)
    local_alignment_parser.add_argument('-x', '--matrix', type=str, default=None)
    local_alignment_parser.add_argument('-ra', '--return-all', action='store_true')

    global_alignment_parser.add_argument('-s1', '--sequence1', type=str, required=True)
    global_alignment_parser.add_argument('-s2', '--sequence2', type=str, required=True)
    global_alignment_parser.add_argument('-m', '--match', type=int, default=1)
    global_alignment_parser.add_argument('-u', '--mismatch', type=int, default=-1)
    global_alignment_parser.add_argument('-g', '--gap-penalty', type=int, default=-2)
    global_alignment_parser.add_argument('-x', '--matrix', type=str, default=None)
    global_alignment_parser.add_argument('-ra', '--return-all', action='store_true')

    fetch_uniprot_parser.add_argument("-f", "--file-path", type=str, required=True)
    fetch_uniprot_parser.add_argument("-o", "--result-path", type=str, required=True)
    fetch_uniprot_parser.add_argument("--full-header", action="store_true", default=False)
    fetch_uniprot_parser.add_argument("-s", "--strict", action = "store_true")
    fetch_uniprot_parser.add_argument("--show-failed", action="store_true", default=False)

    multi_query_search_parser.add_argument("-q", "--query-sequences", type=str, required=True)
    multi_query_search_parser.add_argument("-d", "--database", type=str, required=True)
    multi_query_search_parser.add_argument("-k", "--kmer-size", type=int, default=3)
    multi_query_search_parser.add_argument("-t", "--threshold", type=int, default=1)
    multi_query_search_parser.add_argument("-n", "--top-n-hits", type=int, default=10)
    multi_query_search_parser.add_argument("--regular", dest="indexed", action="store_false", default=True)
    multi_query_search_parser.add_argument(
        "-r",
        "--refine",
        dest="refine_result",
        action="store_true",
        default=False,
        help="Refine selected hits in indexed or regular mode.",
    )
    multi_query_search_parser.add_argument(
        '-m',
        '--match',
        type=int,
        default=1,
        help="Refinement match score when no substitution matrix is selected (default: 1).",
    )
    multi_query_search_parser.add_argument(
        '-u',
        '--mismatch',
        type=int,
        default=-1,
        help="Refinement mismatch score when no substitution matrix is selected (default: -1).",
    )
    multi_query_search_parser.add_argument(
        '-g',
        '--gap-penalty',
        type=int,
        default=-2,
        help="Linear refinement gap penalty in simple or matrix mode (default: -2).",
    )
    multi_query_search_parser.add_argument(
        '-x',
        '--matrix',
        type=str,
        default=None,
        help="Optional substitution-matrix name; replaces simple match/mismatch scoring.",
    )

    args = parser.parse_args()

    if args.command == "search":
        results = search(
            args.query,
            args.database,
            args.kmer_size,
            args.threshold,
            args.top_n_hits,
            args.refine,
            args.match,
            args.mismatch,
            args.gap_penalty,
            args.matrix
        )
        print(json.dumps(results, indent=2))

    elif args.command == "align-local":
        results = local_alignment(
            args.sequence1,
            args.sequence2,
            args.match,
            args.mismatch,
            args.gap_penalty,
            args.matrix,
            args.return_all
        )

        print(json.dumps(results, indent=2))

    elif args.command == "align-global":
        results = global_alignment(
            args.sequence1,
            args.sequence2,
            args.match,
            args.mismatch,
            args.gap_penalty,
            args.matrix,
            args.return_all
        )

        print(json.dumps(results, indent=2))

    elif args.command == "fetch-uniprot":

        result = fetch_uniprot_sequences(args.file_path, args.strict)
        records = result["records"]
        write_fasta_records(records, args.result_path ,args.full_header)
  
        print("Fetch Summary:")
        print("--------------\n")
        print(f"Successfully fetched accessions: {len(result['records'])}")
        print(f"Failed accessions: {len(result['failed'])}")
        print(f"Output written to: {args.result_path}\n")
        if args.show_failed:
            print("Detailed breakdown of failed accessions:")
            print("----------------------------------------\n")
            print(format_failed_accessions(result["failed"]))

    elif args.command == "multi-search":

        result = multi_search(
                query_fasta=args.query_sequences,
                database=args.database,
                k=args.kmer_size,
                threshold=args.threshold,
                top_n_hits=args.top_n_hits,
                 indexed=args.indexed,
                 refinement=args.refine_result,
                 match_score=args.match,
                 mismatch_score=args.mismatch,
                 gap_penalty=args.gap_penalty,
                 matrix=args.matrix,
             )

        print(json.dumps(result, indent=2)) 

if __name__ == "__main__":
    main()
