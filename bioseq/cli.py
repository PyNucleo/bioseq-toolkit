import argparse
import json

from bioseq.pipelines.search_pipeline import search
from bioseq.alignment.smith_waterman import local_alignment
from bioseq.alignment.needleman_wunsch import global_alignment


def main():
    
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest="command", required=True)

    search_parser = subparsers.add_parser("search", help = "Run kmer search")
    local_alignment_parser = subparsers.add_parser("align-local", help = "Run local alignment on two sequences")
    global_alignment_parser = subparsers.add_parser("align-global", help = "Run global alignment on two sequences")

    search_parser.add_argument('-q', '--query', type=str, required=True)
    search_parser.add_argument('-d', '--database', type=str, required=True)
    search_parser.add_argument('-k', '--kmer-size', type=int, default=3)
    search_parser.add_argument('-t', '--threshold',type=int, default=3)
    search_parser.add_argument('-n', '--top-n-hits', type=int, default=10)
    search_parser.add_argument('-r', '--refine', action='store_true')

    local_alignment_parser.add_argument('-s1', '--sequence1', type=str, required=True)
    local_alignment_parser.add_argument('-s2', '--sequence2', type=str, required=True)
    local_alignment_parser.add_argument('-m', '--match', type=int, default=2)
    local_alignment_parser.add_argument('-u', '--mismatch', type=int, default=-1)
    local_alignment_parser.add_argument('-g', '--gap-penalty', type=int, default=-2)
    local_alignment_parser.add_argument('-x', '--matrix', type=str, default=None)
    local_alignment_parser.add_argument('-ra', '--return-all', action='store_true')

    global_alignment_parser.add_argument('-s1', '--sequence1', type=str, required=True)
    global_alignment_parser.add_argument('-s2', '--sequence2', type=str, required=True)
    global_alignment_parser.add_argument('-m', '--match', type=int, default=2)
    global_alignment_parser.add_argument('-u', '--mismatch', type=int, default=-1)
    global_alignment_parser.add_argument('-g', '--gap-penalty', type=int, default=-2)
    global_alignment_parser.add_argument('-x', '--matrix', type=str, default=None)
    global_alignment_parser.add_argument('-ra', '--return-all', action='store_true')

    args = parser.parse_args()

    if args.command == "search":
        results = search(
            args.query,
            args.database,
            args.kmer_size,
            args.threshold,
            args.top_n_hits,
            args.refine
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

if __name__ == "__main__":
    main()