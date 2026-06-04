# virusHanter2 documentation

Long-form documentation for `virusHanter2`. The package overview,
quick-start commands and full-feature run example live in the
top-level [`README.md`](../README.md); this directory holds material
that is too long to keep there.

## For operators

| Topic | File |
|---|---|
| **Database setup — all 9 databases, sources, build commands** | [DATABASE_SETUP.md](DATABASE_SETUP.md) |
| Pipeline stages, the `{assembler}` wildcard, output tree | [PIPELINE.md](PIPELINE.md) |
| Config schema and every opt-in flag | [CONFIGURATION.md](CONFIGURATION.md) |
| Reference databases — production paths, refresh cadence | [REFERENCE_DBS.md](REFERENCE_DBS.md) |
| **Rebuild the classification databases with one snapshot** | [REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md) |
| Per-(sample, virus) CSV schema + multi-run merge | [PER_VIRUS_OUTPUT.md](PER_VIRUS_OUTPUT.md) |

## Conventions and provenance

| Topic | File |
|---|---|
| Parity invariants and intentional breaks against the original `virusHanter` | [PARITY_NOTES.md](PARITY_NOTES.md) |
| Project conventions for AI assistants | [../CLAUDE.md](../CLAUDE.md) |

## Related repositories

- [`reportHanter`](../../reportHanter) — the HTML rendering package.
  See its [`docs/README.md`](../../reportHanter/docs/README.md) for
  the report-side documentation.

## Testing and smoke

- [Test harness](../test/README.md) — smoke fixtures and end-to-end
  smoke runner.
