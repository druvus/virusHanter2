"""
Pipeline-side helpers for virusHanter2.

These are the small utilities that are called directly from the Snakemake
rule files: sample discovery, FASTQ/FASTA wrangling, Kraken2 post-processing,
and the BLASTN driver. All HTML-report and plotting helpers that used to live
here have moved into the reportHanter package (`reporthanter` on PyPI), which
the report rules invoke either by CLI or by library import.
"""

from pathlib import Path
import getpass
import hashlib
import os
import pickle
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd


def _cache_user() -> str:
    """Best-effort current-user tag for namespacing the cache dir.

    Falls back to the numeric uid (or "user") if the username cannot be
    resolved, so the call never raises in an odd execution environment.
    """
    try:
        return getpass.getuser()
    except Exception:  # noqa: BLE001 - any resolution failure -> fallback
        return str(os.getuid()) if hasattr(os, "getuid") else "user"


def _taxdump_cache_path(kind: str, *sources: str) -> Path:
    """Per-user cache path for a parsed taxdump structure.

    The cache key folds in each source path together with its mtime and
    size, so a refreshed taxdump (new mtime/size) misses the stale cache
    and is re-parsed automatically. The cache lives under the OS
    cache/temp dir, never on the (often read-only, often external and
    slow) reference-database volume.
    """
    parts: list[str] = [kind]
    for source in sources:
        try:
            st = os.stat(source)
            # Full-resolution mtime (nanoseconds): an in-place refresh to
            # an identical byte size within the same wall-clock second
            # still invalidates the cache, which truncated whole-second
            # mtimes would miss.
            parts.append(f"{os.path.abspath(source)}:{st.st_mtime_ns}:{st.st_size}")
        except OSError:
            parts.append(f"{source}:0:0")
    digest = hashlib.sha1("|".join(parts).encode()).hexdigest()[:16]
    base = os.environ.get("XDG_CACHE_HOME") or tempfile.gettempdir()
    # Namespace the cache dir by user so that on a shared HPC node one
    # user cannot collide with (or be unable to overwrite) another's
    # cache files, and an unpickled object always came from this user.
    return Path(base) / f"virushanter2_taxdump_{_cache_user()}" / f"{kind}_{digest}.pkl"


def _cache_load(cache_path: Path):
    """Return the unpickled cache object, or None on any miss/corruption.

    Catches broadly on purpose: a truncated, version-mismatched or
    otherwise unreadable pickle can raise more than the usual
    Unpickling/EOF/Value errors (e.g. AttributeError, ModuleNotFoundError,
    MemoryError). Any failure simply falls back to a fresh parse rather
    than aborting the rule.
    """
    try:
        with open(cache_path, "rb") as fh:
            return pickle.load(fh)
    except Exception:  # noqa: BLE001 - any cache failure -> re-parse
        return None


def _cache_store(cache_path: Path, obj) -> None:
    """Atomically pickle `obj` to `cache_path`; silent on any failure.

    Writes to a pid-suffixed temp file then `os.replace` so concurrent
    Snakemake jobs parsing the same taxdump cannot leave a torn cache
    for a reader (replace is atomic; identical content makes the race
    benign).
    """
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache_path.with_name(f"{cache_path.name}.tmp.{os.getpid()}")
        with open(tmp, "wb") as fh:
            pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, cache_path)
    except OSError:
        pass

# ``pyfastx`` is used only by ``fastx_file_to_df`` (called from the
# wrangle_pilon rule's panel env). Lazy-import inside that function so
# the rest of the module (the BLAST / classifier canonicalisers used
# by the kaiju / bwa / blastn envs) does not pull pyfastx as a
# hard dependency.


_STRAIN_LIKE_MARKERS = (
    " type ",
    " strain ",
    " isolate ",
    " genotype ",
    " serotype ",
    " subtype ",
)


def _is_strain_like_name(name: str) -> bool:
    """Heuristic for "this taxid is a strain/type/isolate sub-entry".

    NCBI's viral taxonomy carries many "no rank" taxa below the species
    level that exist purely to register a sequenced isolate or a
    serologically distinct type (e.g. ``Human herpesvirus 4 type 2``,
    ``Human alphaherpesvirus 1 strain F``). Detect them by looking for
    the marker words separated by spaces so we do not match substrings
    in the middle of legitimate species names. Case-insensitive.
    """
    if not name:
        return False
    padded = " " + name.lower() + " "
    return any(marker in padded for marker in _STRAIN_LIKE_MARKERS)


# Name categories from names.dmp surfaced as "aliases" alongside the
# scientific name. Scientists need the old / common / acronym forms
# to recognise a species that NCBI / ICTV has recently renamed.
_ALIAS_NAME_CLASSES = (
    "acronym",
    "common name",
    "genbank common name",
    "equivalent name",
    "synonym",
    "genbank synonym",
    "genbank acronym",
)


def parse_nodes_dmp(path: str | Path) -> dict[int, tuple[int, str]]:
    """Parse NCBI ``nodes.dmp`` into ``{tax_id: (parent_tax_id, rank)}``.

    The NCBI ``.dmp`` files use ``\\t|\\t`` as the inter-field delimiter
    and ``\\t|`` at the end of each row. The first three fields are
    ``tax_id``, ``parent_tax_id``, ``rank``; everything after column 3
    is ignored. Returns an empty dict if ``path`` does not exist.
    """
    out: dict[int, tuple[int, str]] = {}
    if not Path(path).exists():
        return out
    cache_path = _taxdump_cache_path("nodes", str(path))
    cached = _cache_load(cache_path)
    if cached is not None:
        return cached
    with open(path) as fh:
        for line in fh:
            stripped = line.rstrip("\n").rstrip("\t|").rstrip()
            parts = stripped.split("\t|\t")
            if len(parts) < 3:
                continue
            try:
                tid = int(parts[0].strip())
                parent = int(parts[1].strip())
            except ValueError:
                continue
            out[tid] = (parent, parts[2].strip())
    _cache_store(cache_path, out)
    return out


def find_genus_taxid(
    tid: int,
    nodes: dict[int, tuple[int, str]],
    *,
    depth_limit: int = 20,
) -> int:
    """Walk the parent chain from ``tid`` and return the first ancestor
    whose rank is ``genus``. Returns ``0`` if no genus ancestor is found
    within ``depth_limit`` steps or if the chain breaks.

    The depth limit guards against retired-taxid cycles that occasionally
    appear in NCBI dumps. The chain also stops at the root (parent == tid
    for the synthetic root node 1).
    """
    if tid <= 0 or not nodes:
        return 0
    seen: set[int] = set()
    current = tid
    for _ in range(depth_limit):
        if current in seen:
            return 0
        seen.add(current)
        node = nodes.get(current)
        if node is None:
            return 0
        parent, rank = node
        if rank == "genus":
            return current
        if parent == current:
            return 0
        current = parent
    return 0


def find_species_taxid(
    tid: int,
    nodes: dict[int, tuple[int, str]],
    *,
    depth_limit: int = 20,
) -> int:
    """Walk the parent chain from ``tid`` and return the first ancestor
    whose rank is ``species``. Returns ``0`` if no species ancestor is
    found within ``depth_limit`` steps.

    A taxid that is itself at species rank returns its own value, so
    the caller can call this unconditionally and use the result as the
    canonical species taxid.
    """
    if tid <= 0 or not nodes:
        return 0
    seen: set[int] = set()
    current = tid
    for _ in range(depth_limit):
        if current in seen:
            return 0
        seen.add(current)
        node = nodes.get(current)
        if node is None:
            return 0
        parent, rank = node
        if rank == "species":
            return current
        if parent == current:
            return 0
        current = parent
    return 0


def _load_taxdump_for_species_walkup(
    nodes_dmp: str | None, names_dmp: str | None
) -> tuple[dict[int, tuple[int, str]], dict[int, str], dict[int, list[str]]] | None:
    """Parse nodes.dmp + names.dmp for the species walkup helpers.

    Returns ``(node_info, sci_name, alias_name)`` where ``alias_name``
    maps each tax_id to the deduplicated list of its non-scientific
    names from ``names.dmp`` (acronym, common name, genbank common
    name, equivalent name, synonym, genbank synonym, genbank acronym).
    Returns ``None`` if either path is missing; callers degrade to a
    no-op in that case rather than raising.
    """
    if (
        not nodes_dmp
        or not names_dmp
        or not Path(nodes_dmp).is_file()
        or not Path(names_dmp).is_file()
    ):
        return None
    cache_path = _taxdump_cache_path("walkup", str(nodes_dmp), str(names_dmp))
    cached = _cache_load(cache_path)
    if cached is not None:
        return cached
    node_info = parse_nodes_dmp(nodes_dmp)
    sci_name: dict[int, str] = {}
    alias_name: dict[int, list[str]] = {}
    alias_seen: dict[int, set[str]] = {}
    with open(names_dmp) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            name_class = parts[3]
            try:
                tid = int(parts[0])
            except ValueError:
                continue
            value = parts[1]
            if name_class == "scientific name":
                sci_name[tid] = value
            elif name_class in _ALIAS_NAME_CLASSES:
                bucket = alias_seen.setdefault(tid, set())
                if value not in bucket:
                    bucket.add(value)
                    alias_name.setdefault(tid, []).append(value)
    result = (node_info, sci_name, alias_name)
    _cache_store(cache_path, result)
    return result


def _format_aliases(
    tid: int,
    species_tid: int,
    sci_name: dict[int, str],
    alias_name: dict[int, list[str]],
) -> str:
    """Return a deduplicated "name1; name2" alias string for ``tid``.

    Collects aliases from both the row's own tax_id and its species
    ancestor (when distinct), then drops the canonical species
    scientific name itself so it does not repeat alongside the
    binomial. Order: row-own aliases first, then species-ancestor
    aliases; case-insensitive dedupe.
    """
    seen_lower: set[str] = set()
    out: list[str] = []
    canonical = sci_name.get(species_tid or tid, "")
    if canonical:
        seen_lower.add(canonical.lower())
    # Pull both the row's own scientific name (the legacy NCBI name
    # scientists still recognise — e.g. ``human gammaherpesvirus 4``
    # for taxid 10376 when its species ancestor is the binomial
    # ``Lymphocryptovirus humangamma4``) and every non-scientific
    # alias category from both the row taxid and the species
    # ancestor. Case-insensitive dedupe; preserve discovery order.
    for source_tid in (tid, species_tid):
        if not source_tid:
            continue
        own_sci = sci_name.get(source_tid, "")
        if own_sci and own_sci.lower() not in seen_lower:
            seen_lower.add(own_sci.lower())
            out.append(own_sci)
        for value in alias_name.get(source_tid, []):
            key = value.lower()
            if key in seen_lower:
                continue
            seen_lower.add(key)
            out.append(value)
    return "; ".join(out)


def canonicalise_taxon_names(
    df: pd.DataFrame,
    *,
    taxid_col: str,
    name_col: str,
    nodes_dmp: str | None,
    names_dmp: str | None,
    raw_suffix: str = "_raw",
) -> pd.DataFrame:
    """Rewrite ``name_col`` with the ICTV-binomial species name for
    every row's ``taxid_col``.

    Walks the parent chain to the first species-rank ancestor (per
    nodes.dmp) and substitutes that ancestor's scientific name (per
    names.dmp). The original value is preserved in
    ``<name_col><raw_suffix>``. Rows whose taxid is missing from the
    taxdump pass through unchanged.

    Generic over the classifier: pass ``("taxonomy_id", "name")`` for
    the Kraken report or ``("taxon_id", "taxon_name")`` for Kaiju's
    table. A missing taxdump degrades to a no-op (raw column added,
    no rewriting).
    """
    if df.empty or name_col not in df.columns or taxid_col not in df.columns:
        return df

    out = df.copy()
    raw_col = f"{name_col}{raw_suffix}"
    out[raw_col] = out[name_col].astype(str)

    parsed = _load_taxdump_for_species_walkup(nodes_dmp, names_dmp)
    if parsed is None:
        return out
    node_info, sci_name, alias_name = parsed

    species_cache: dict[int, int] = {}

    def _species_taxid(tid: int) -> int:
        if tid in species_cache:
            return species_cache[tid]
        seen: set[int] = set()
        current = tid
        for _ in range(64):
            if current in seen:
                break
            seen.add(current)
            node = node_info.get(current)
            if node is None:
                break
            parent_tid, rank = node
            if rank == "species":
                species_cache[tid] = current
                return current
            if parent_tid == current or parent_tid == 1:
                break
            current = parent_tid
        species_cache[tid] = 0
        return 0

    new_names: list[str] = []
    aliases: list[str] = []
    for raw_name, raw_taxid in zip(
        out[name_col].astype(str), out[taxid_col], strict=False
    ):
        try:
            tid = int(raw_taxid)
        except (TypeError, ValueError):
            new_names.append(raw_name)
            aliases.append("")
            continue
        species_tid = _species_taxid(tid)
        canonical = sci_name.get(species_tid, "") if species_tid else ""
        new_names.append(canonical if canonical else raw_name)
        aliases.append(_format_aliases(tid, species_tid, sci_name, alias_name))
    out[name_col] = new_names
    out["aliases"] = aliases
    return out


def dummy_contig_sentinel(
    merged: pd.DataFrame,
    checkv_df: pd.DataFrame,
    assembler: str,
) -> pd.DataFrame:
    """Carry a DUMMY_CONTIG sentinel row through a failed assembly.

    When an assembler produces no usable contigs it writes a single
    ``DUMMY_CONTIG`` sequence (Pilon polishes it to
    ``DUMMY_CONTIG_pilon``). CheckV still reports that contig, but
    BLASTN drops it because it matches no database entry, so the
    BLAST/CheckV inner join in ``merge_checkv_blastn`` is empty. An
    empty merge is indistinguishable from a genuine negative (real
    contigs assembled, none viral), so ``per_virus_metrics`` cannot
    flag the silent failure.

    When the merge is empty *and* every contig CheckV reported is a
    DUMMY_CONTIG sentinel, return a one-row frame carrying that name
    (with the merge's own columns preserved so downstream
    concatenation stays clean) so the empty-assembly note can fire.
    Otherwise return ``merged`` unchanged.
    """
    if not merged.empty or checkv_df.empty or "name" not in checkv_df.columns:
        return merged
    names = checkv_df["name"].dropna().astype(str)
    if len(names) == 0 or not names.str.startswith("DUMMY_CONTIG").all():
        return merged
    row = {col: pd.NA for col in merged.columns}
    row["name"] = names.iloc[0]
    if "assembler" in row:
        row["assembler"] = assembler
    return pd.DataFrame([row], columns=list(merged.columns))


def canonicalise_blast_match_name(
    blastn_df: pd.DataFrame,
    parquet_df: pd.DataFrame,
    nodes_dmp: str | None,
    names_dmp: str | None,
) -> pd.DataFrame:
    """Rewrite the BLAST ``match_name`` column with the canonical
    NCBI scientific name of the lowest non-strain-like ancestor.

    NCBI's RefSeq keeps multiple records per ICTV species when a
    sequenced isolate or serologically distinct type was registered
    as its own taxon: ``NC_007605`` (EBV-1, taxid 10376, rank S1,
    name ``human gammaherpesvirus 4``) and ``NC_009334`` (EBV-2,
    taxid 12509, rank S2, name ``Human herpesvirus 4 type 2``) are
    the archetypal pair. The two records appear as separate bars in
    the Assembly classification chart even though they are the same
    ICTV species; this function collapses them by walking the
    ``nodes.dmp`` parent chain until reaching the first ancestor
    whose rank is ``species``, and using that taxid's scientific
    name as the canonical label. For the EBV pair both records walk
    up to taxid 3050299 (rank species, ICTV-binomial name
    ``Lymphocryptovirus humangamma4``).

    Behaviour:
    - Lookup the BLAST hit's accession in ``parquet_df`` -> ``tax_id``.
      Try the raw accession first, then the version-stripped form.
    - Walk up via ``nodes.dmp`` parent pointers until an ancestor
      whose rank is ``species`` is reached. Fall back to the
      strain-marker heuristic for taxids that have no species-rank
      ancestor (older NCBI viral taxonomy entries that NCBI has not
      yet promoted under the binomial scheme).
    - Replace ``match_name`` with the canonical scientific name; keep
      the raw BLAST title in a new ``match_name_raw`` column for
      audit.
    - Rows whose accession is absent from the parquet, or whose taxid
      is missing from the taxdump, pass through unchanged with
      ``match_name_raw`` set equal to ``match_name``.

    A missing ``nodes.dmp`` / ``names.dmp`` path (config did not
    point at a taxdump) degrades to a no-op: the merged CSV gets a
    ``match_name_raw`` column duplicating ``match_name``, but no
    rewriting happens.
    """
    if blastn_df.empty or "match_name" not in blastn_df.columns:
        return blastn_df

    out = blastn_df.copy()
    out["match_name_raw"] = out["match_name"].astype(str)

    if (
        not nodes_dmp
        or not names_dmp
        or not Path(nodes_dmp).is_file()
        or not Path(names_dmp).is_file()
    ):
        return out

    acc_to_tax = parquet_accession_to_taxid(parquet_df)

    parsed = _load_taxdump_for_species_walkup(nodes_dmp, names_dmp)
    if parsed is None:
        return out
    node_info, sci_name, alias_name = parsed

    def _canonical_taxid(tid: int) -> int:
        """Walk up to the first species-rank ancestor.

        If no species-rank ancestor exists within the depth cap (rare
        — typically only for old NCBI viral genera whose children
        have no species row yet), fall back to the legacy strain-name
        heuristic so we never make the label worse than the raw input.
        """
        seen: set[int] = set()
        current = tid
        for _ in range(64):
            if current in seen:
                break
            seen.add(current)
            node = node_info.get(current)
            if node is None:
                break
            parent_tid, rank = node
            if rank == "species":
                return current
            if parent_tid == current or parent_tid == 1:
                break
            current = parent_tid

        # Fallback: legacy strain-marker heuristic. Walk up while the
        # name carries a strain / type / isolate marker. Preserves
        # behaviour for taxids whose species ancestor is missing from
        # the local taxdump snapshot.
        seen.clear()
        current = tid
        for _ in range(64):
            if current in seen:
                break
            seen.add(current)
            name = sci_name.get(current, "")
            if not _is_strain_like_name(name):
                return current
            node = node_info.get(current)
            if node is None:
                return current
            parent_tid, _rank = node
            if parent_tid == current or parent_tid == 1:
                return current
            current = parent_tid
        return current

    canonical_cache: dict[int, tuple[str, str]] = {}

    def _name_and_aliases_for_accession(accession: str) -> tuple[str, str] | None:
        if not accession:
            return None
        # Try versioned first, then version-stripped.
        tid = acc_to_tax.get(accession) or acc_to_tax.get(
            accession.split(".")[0]
        )
        if tid is None:
            return None
        if tid in canonical_cache:
            return canonical_cache[tid]
        cid = _canonical_taxid(tid)
        name = sci_name.get(cid, "")
        aliases = _format_aliases(tid, cid, sci_name, alias_name)
        canonical_cache[tid] = (name, aliases)
        return name, aliases

    canonical_names: list[str] = []
    canonical_aliases: list[str] = []
    for raw_match, accession in zip(
        out["match_name"].astype(str),
        out.get("accession", pd.Series([""] * len(out))).astype(str),
        strict=False,
    ):
        result = _name_and_aliases_for_accession(accession)
        if result and result[0]:
            canonical_names.append(result[0])
            canonical_aliases.append(result[1])
        else:
            canonical_names.append(raw_match)
            canonical_aliases.append("")
    out["match_name"] = canonical_names
    out["aliases"] = canonical_aliases

    return out


def parquet_accession_to_taxid(parquet_df: pd.DataFrame) -> dict[str, int]:
    """Bucket VIRUS_PARQUET rows into a base-accession -> tax_id map.

    The parquet's ``name`` column carries the original FASTA header
    (e.g. ``NC_007605.1 Human gammaherpesvirus 4 ...``); BLAST hits
    surface either the versioned (``NC_007605.1``) or unversioned
    (``NC_007605``) accession depending on the BLAST output format.
    Index both so a downstream lookup always lands the same tax_id
    regardless of which form the caller has.

    Skips rows with ``tax_id == 0`` (the build_virus_parquet sentinel
    for accessions that the accession2taxid mapping could not resolve)
    so they do not contaminate the lookup.
    """
    out: dict[str, int] = {}
    if "name" not in parquet_df.columns or "tax_id" not in parquet_df.columns:
        return out
    for row in parquet_df.itertuples():
        try:
            tid = int(row.tax_id)
        except (ValueError, TypeError):
            continue
        if tid == 0:
            continue
        first_token = str(row.name).split()[0]
        out[first_token] = tid
        out[first_token.split(".")[0]] = tid
    return out


# Assembler helpers live in the dependency-free ``assembler_utils`` module
# so the assembler scripts can import them inside their minimal per-rule
# conda envs (which lack pandas/numpy). Re-exported here for driver-env
# callers and the test-suite that import them from ``scripts.functions``.
from scripts.assembler_utils import (  # noqa: E402
    assembler_max_attempts,  # noqa: F401  - intentional re-export
    write_dummy_contig,  # noqa: F401  - intentional re-export
)


def read_file_as_blob(file_path: str) -> str:
    """Return the file content as a lowercase hex string.

    Used by `aggregate_run_information` to embed per-sample HTML reports as
    a column in the run summary CSV (matches the original virusHanter
    behavior).
    """
    with open(file_path, "rb") as f:
        return f.read().hex()


def common_suffix(folder: str) -> str:
    """Find the longest filename suffix shared by every sequencing file in
    `folder`. Used to derive the read1/read2 file extension at workflow
    setup time.
    """
    samples = sorted(
        file.name
        for file in Path(folder).iterdir()
        if re.search(r"\.(fq|fastq|fa|fasta|fna)(\.gz)?$", file.name)
    )
    if not samples:
        return ""

    test_sample = samples[0]
    suffix = ""
    for i in range(1, len(test_sample) + 1):
        index = -i
        if any(sample[index] != test_sample[index] for sample in samples):
            break
        suffix += test_sample[index]
    return suffix[::-1]


def paired_reads(folder: str) -> list:
    """Return the common-prefix sample names for paired-end FASTQ files in
    `folder`. Assumes files are sorted by name and pair up two-by-two.
    """
    def common_prefix(a: str, b: str) -> str:
        out = ""
        for ca, cb in zip(a, b):
            if ca != cb:
                break
            out += ca
        return out

    samples = sorted(
        x.stem
        for x in Path(folder).iterdir()
        if re.search(r"\.(fq|fastq|fa|fasta|fna)(\.gz)?$", x.name)
    )

    if len(samples) % 2 != 0:
        raise ValueError(
            f"Found an odd number of sequencing files ({len(samples)}) in "
            f"{folder}; paired-end discovery needs an R1/R2 mate for every "
            f"sample. Check for a missing mate or a stray FASTQ. Files: "
            f"{samples}"
        )

    prefixes = []
    for i in range(0, len(samples), 2):
        prefixes.append(common_prefix(samples[i], samples[i + 1]))
    return prefixes


def kaiju_db_files(kaiju_db: str) -> tuple:
    """Locate the `.fmi`, `names.dmp`, and `nodes.dmp` files inside a Kaiju
    database directory.

    Returns placeholder paths under the given directory when the directory
    itself does not yet exist. This lets workflow construction (including
    `snakemake --lint` and `snakemake -n`) proceed even before the database
    has been materialized; the rule that consumes these paths will still
    fail loudly at run time if the files are not present.
    """
    db_path = Path(kaiju_db)
    if not db_path.is_dir():
        return (
            db_path / "kaiju_db.fmi",
            db_path / "names.dmp",
            db_path / "nodes.dmp",
        )

    files = [x for x in db_path.iterdir() if x.is_file()]
    fmi = next((x for x in files if x.suffix == ".fmi"), db_path / "kaiju_db.fmi")
    names = next((x for x in files if x.name == "names.dmp"), db_path / "names.dmp")
    nodes = next((x for x in files if x.name == "nodes.dmp"), db_path / "nodes.dmp")
    return fmi, names, nodes


def fastx_file_to_df(fastx_file: str) -> pd.DataFrame:
    """Read a FASTA/FASTQ file into a DataFrame sorted by sequence length.

    pyfastx 2.x yields tuples from `Fastx`; the first two fields are
    (name, sequence). Earlier versions returned attribute-bearing objects;
    callers in the original virusHanter targeted that older API.
    """
    import pyfastx  # noqa: PLC0415 - lazy import; see module-level note

    fastx = pyfastx.Fastx(fastx_file)
    rows = [(record[0], record[1]) for record in fastx]
    if not rows:
        return pd.DataFrame(columns=["name", "sequence", "read_len"])

    names, seqs = zip(*rows)
    return (
        pd.DataFrame({"name": list(names), "sequence": list(seqs)})
        .assign(read_len=lambda x: x.sequence.str.len())
        .sort_values("read_len", ascending=False)
    )


def wrangle_kraken(kraken_file: str) -> pd.DataFrame:
    """Parse a Kraken2 report TSV into a DataFrame with an explicit `domain`
    column carried down from the nearest D/U/R/R1 parent row.

    Kraken2's pluspf database tags 'Viruses' as ``tax_lvl='D'``, so the
    original D/U/R anchor set was sufficient. The smaller viral-only
    databases (e.g. ``k2_viral_*``) place 'Viruses' at ``tax_lvl='R1'``
    because there is no other superkingdom alongside it to require a
    domain level. Including ``R1`` as an anchor makes the ``domain``
    column carry "Viruses" in both DB shapes. On pluspf the only extra
    R1 row is "cellular organisms", which is immediately overridden by
    the next D row (Bacteria), so the column values for every
    D-and-below row stay identical and the parity invariant holds.
    """
    kraken = (
        pd.read_csv(
            kraken_file,
            sep="\t",
            header=None,
            names=["percent", "count_clades", "count", "tax_lvl", "taxonomy_id", "name"],
        )
        .assign(name=lambda x: x.name.str.strip())
        .assign(
            domain=lambda x: np.select(
                [x.tax_lvl.isin(["D", "U", "R", "R1"])],
                [x.name],
                default=pd.NA,
            )
        )
    )
    kraken["domain"] = kraken["domain"].ffill()
    return kraken


def run_blastn(contigs_csv: str, db: str, temp_file: str, threads: int) -> pd.DataFrame:
    """Run blastn (megablast) on each row of a contigs CSV one at a time and
    return the input table joined with the best hit per contig.
    """
    # `BLASTDB` should point at the directory that holds the database files
    # (and any auxiliary `taxdb.*` files), not at the database prefix
    # itself. The original virusHanter pulled this from a dedicated
    # `BLASTDB_ENVIRON_VARIABLE` config entry; we derive it from the db
    # path so we do not need a second config key.
    os.environ["BLASTDB"] = str(Path(db).parent)
    df = pd.read_csv(contigs_csv)
    if df.empty:
        return df

    matches = []
    for contig in df.itertuples():
        with open(temp_file, "w") as f:
            f.write(f">{contig.name}\n{contig.sequence}\n")
        command = [
            "blastn", "-num_threads", str(threads), "-task", "megablast",
            "-query", temp_file, "-db", db, "-max_target_seqs", "1",
            "-outfmt", "6 stitle sacc pident slen",
        ]
        # capture_output so a non-zero blastn includes stderr in the raised
        # CalledProcessError rather than printing to the main process stdout.
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
        )
        matches.append(result.stdout.strip())

    df = df.assign(matches=matches).loc[lambda x: x.matches != ""]
    if df.empty:
        return df

    df[["match_name", "accession", "percent_identity", "sequence_len"]] = (
        df.matches.str.split("\t", expand=True).iloc[:, :4]
    )
    df = df.assign(sequence_len=lambda x: x.sequence_len.str.split("\n").str[0])
    return df
