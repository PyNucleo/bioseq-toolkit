from bioseq.search.kmer_search import kmer_search
from bioseq.pipelines.search_pipeline import search

from database.database_utils import normalize_database
import time

def benchmark_kmer_search_run_time(file, k = 3, threshold = 3):

    start = time.perf_counter()

    query = "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels"

    db = normalize_database(file)

    kmer_search(query, db, k, threshold)

    end = time.perf_counter()

    interval = end - start

    return interval

def benchmark_kmer_refinement_run_time(file, k=3, threshold=1, top_n_hits=10, refinement=True):
    
    start = time.perf_counter()

    query = "ttltnpqkaairsswskfmdngvsngqgfymdlfkahpetltpfkslfggltlaqlqdnpkmkaqslvfcngmssfvdhlddndmlvvliqkmaklhnnrgirasdlrtaydilihymedhnhmvggakdawevfvgficktlgdymkels"

    search(query, file, k, threshold, top_n_hits, refinement)

    end = time.perf_counter()

    interval = end - start
    
    return interval

