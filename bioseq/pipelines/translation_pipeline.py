from ..fasta_io import read_fasta_records
from ..sequence_utils import gc_content, mrna_template
from ..translation import translate_sequence
from ..validators import report_invalid_symbols_and_positions


def process_fasta_sequences(file_path):
    file_records = read_fasta_records(file_path)
    accepted = []
    rejected = []

    for record_position, record in enumerate(file_records, start=1):
        sequence = record["sequence"].upper()
        invalid_report = report_invalid_symbols_and_positions(sequence)

        if invalid_report["invalid_positions"]:
            rejected.append({
                "id": record["id"],
                "record_position": record_position,
                "sequence": sequence,
                "reason_code": "unsupported_dna_symbols",
                "reason": "Sequence contains unsupported DNA symbols.",
                **invalid_report,
            })
            continue

        mrna_sequence = mrna_template(sequence)
        accepted.append({
            "id": record["id"],
            "record_position": record_position,
            "sequence": sequence,
            "length": len(sequence),
            "gc_content": gc_content(sequence),
            "transcribed_strand": mrna_sequence,
            "amino_acid_chain": str(translate_sequence(mrna_sequence)),
        })

    return {
        "accepted": accepted,
        "rejected": rejected,
        "summary": {
            "total_records": len(accepted) + len(rejected),
            "accepted_records": len(accepted),
            "rejected_records": len(rejected),
        },
    }
