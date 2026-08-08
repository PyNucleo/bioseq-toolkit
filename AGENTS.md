## graphify

This project has a local knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

When the user types `/graphify`, use the installed graphify skill or instructions before doing anything else.

Rules:
- Use Graphify first for broad, unfamiliar, cross-module, architectural, dependency, and blast-radius exploration. When graphify-out/graph.json exists, use `graphify query "<question>"` for orientation; use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts.
- For narrow work where the relevant files or functions are already known, direct repository inspection is appropriate and Graphify is optional.
- Graphify is an orientation aid only. Current repository source and tests remain authoritative, and material Graphify conclusions must be verified against source before making changes.
- Dirty graphify-out/ files are expected after hooks or incremental updates; they are not a reason to skip Graphify when it is useful.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
