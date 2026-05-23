from bioseq.pipelines.search_pipeline import search

results = search(
    query="ATGCG",
    database=[
        "ATGCGT",
        "ATGCGA",
        "GGGGGG"
    ],
    k=3
)

print(results)