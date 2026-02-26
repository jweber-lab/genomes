#!/usr/bin/env python3
"""
setup_blast_db.py — Build and manage local BLAST databases from downloaded
genome assemblies, annotations, and protein sequences.

Example
-------
    python setup_blast_db.py --csv bioprojects.csv build --alias
    python setup_blast_db.py --csv bioprojects.csv --dry-run build
    python setup_blast_db.py --csv bioprojects.csv update --alias
    python setup_blast_db.py --csv bioprojects.csv alias
    python setup_blast_db.py --csv bioprojects.csv list

Subcommands
-----------
  build   Build nucleotide and protein BLAST databases for every non-skipped
          row in the CSV. Optionally create combined alias databases with
          --alias.
  update  Identical to build — rebuilds all databases from current downloads.
  alias   Create alias databases (all_nucl, all_prot) from existing per-label
          databases without rebuilding anything.
  list    Print the contents of blast_db/registry.yaml.

Workflow (build)
----------------
  1. Read bioprojects.csv; skip rows where skip=true.
  2. For each label, locate genome FASTA, GFF3, and protein FASTA files in the
     downloads directory. Searches WormBase ParaSite, NCBI, and ENA
     subdirectories in that order, so WormBase's curated files take priority
     when present. Prefers manifest.json paths if available.
  3. Concatenate genome FASTAs and run makeblastdb to create a nucleotide
     database (blast_db/<label>/nucl/). Copy the first GFF3 alongside it.
  4. If protein FASTAs are found, concatenate and run makeblastdb to create a
     protein database (blast_db/<label>/prot/).
  5. If --alias is set, create combined alias databases (all_nucl, all_prot)
     using blastdb_aliastool so all genomes can be searched in a single query.
  6. Write blast_db/registry.yaml with paths, metadata, and supported BLAST
     programs for each database. Downstream pipelines read this registry.

Builds can run in parallel across labels (set build_jobs in config.yml).

Dependencies: BLAST+ suite (makeblastdb, blastdbcmd, blastdb_aliastool),
PyYAML (optional, falls back to JSON for registry).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:
    yaml = None

LOG_FILE = "blast_db.log"
DEFAULT_BLAST_DB_ROOT = "blast_db"


def detect_project_source(bioproject: str) -> str:
    """Detect project source from BioProject prefix. Returns 'ENA', 'DDBJ', or 'NCBI'."""
    bioproject_upper = bioproject.upper()
    if bioproject_upper.startswith(("PRJEB", "PRJEA")):
        return "ENA"
    elif bioproject_upper.startswith("PRJDB"):
        return "DDBJ"
    else:
        return "NCBI"
DEFAULT_DOWNLOAD_ROOT = "downloads"
DEFAULT_CSV = "bioprojects.csv"
MANIFEST_NAME = "manifest.json"
NCBI_DATA_DIR = "ncbi_dataset"
DATA_SUBDIR = "data"


def find_ncbi_data_dir(label_dir: Path) -> Path | None:
    """Find NCBI dataset data directory."""
    d = label_dir / NCBI_DATA_DIR / DATA_SUBDIR
    return d if d.is_dir() else None


def find_ena_data_dir(label_dir: Path) -> Path | None:
    """Find ENA download directory. enaDataGet creates a dir named after the accession."""
    if not label_dir.exists() or not label_dir.is_dir():
        return None
    # enaDataGet creates a directory named after the accession (e.g., GCA_000469805.3)
    # Look for common ENA directory patterns
    try:
        for item in label_dir.iterdir():
            if item.is_dir() and (item.name.startswith("GCA_") or item.name.startswith("GCF_")):
                # Check if it contains sequence files
                if any(item.rglob("*.fa*")) or any(item.rglob("*.embl")):
                    return item
    except OSError:
        return None
    return None


def find_wormbase_data_dir(label_dir: Path) -> Path | None:
    """Return the wormbase_parasite subdirectory if it exists."""
    d = label_dir / "wormbase_parasite"
    return d if d.is_dir() else None


REGISTRY_NAME = "registry.yaml"
STAGING_SUFFIX = ".staging"
MAKEBLASTDB = "makeblastdb"
BLASTDBCMD = "blastdbcmd"
BLASTDB_ALIASTOOL = "blastdb_aliastool"
NUCL_DIR = "nucl"
PROT_DIR = "prot"
SEQUENCES_NUCL = "sequences.fna"
SEQUENCES_PROT = "sequences.faa"
ANNOTATIONS_GFF = "annotations.gff"
DB_BASENAME_NUCL = "sequences"
DB_BASENAME_PROT = "sequences"


def load_config(config_path: Path | None) -> dict[str, Any]:
    out = {
        "download_root": DEFAULT_DOWNLOAD_ROOT,
        "blast_db_root": DEFAULT_BLAST_DB_ROOT,
        "build_jobs": 1,
        "max_file_sz": None,
    }
    if config_path is None or not config_path.exists():
        return out
    if yaml is None:
        return out
    try:
        with open(config_path, "r") as f:
            data = yaml.safe_load(f) or {}
        out["download_root"] = data.get("download_root", out["download_root"])
        out["blast_db_root"] = data.get("blast_db_root", out["blast_db_root"])
        out["build_jobs"] = int(data.get("build_jobs", out["build_jobs"]) or 1)
        out["max_file_sz"] = data.get("max_file_sz")
    except Exception as e:
        logging.warning("Failed to load config.yml: %s; using defaults", e)
    return out


def read_csv(csv_path: Path) -> list[dict[str, Any]]:
    """Read bioprojects CSV; normalize column names and defaults. Filter skipped rows."""
    rows = []
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = {k.strip().lower().replace(" ", "_"): v.strip() if isinstance(v, str) else v for k, v in row.items()}
            bioproject = (r.get("bioproject") or "").strip()
            if not bioproject:
                continue
            skip_val = (r.get("skip") or "false").strip().lower()
            if skip_val in ("true", "1", "yes"):
                continue
            label = (r.get("label") or bioproject).strip()
            r["bioproject"] = bioproject
            r["label"] = label
            r["source"] = detect_project_source(bioproject)
            inc_prot = (r.get("include_protein") or "true").strip().lower()
            r["include_protein"] = inc_prot in ("true", "1", "yes")
            rows.append(r)
    return rows


def get_manifest_paths(download_root: Path, label: str) -> dict[str, Any] | None:
    manifest_path = download_root / label / MANIFEST_NAME
    if not manifest_path.exists():
        return None
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)


def discover_paths(download_root: Path, label: str, include_protein: bool) -> dict[str, list[Path]]:
    """Discover genome FASTA, GFF, and optionally protein FASTA.
    Supports NCBI, ENA, and WormBase ParaSite download structures."""
    label_dir = download_root / label
    out: dict[str, list[Path]] = {"genome_fasta": [], "gff": [], "protein_fasta": []}
    seen: set[Path] = set()

    def _add(key: str, path: Path) -> None:
        if path not in seen:
            seen.add(path)
            out[key].append(path)

    # WormBase first so its curated annotations/proteins are preferred (appear first in lists)
    for d in (find_wormbase_data_dir(label_dir), find_ncbi_data_dir(label_dir), find_ena_data_dir(label_dir)):
        if d is None:
            continue
        for f in d.rglob("*.fna"):
            if f.is_file() and f.stat().st_size > 0:
                _add("genome_fasta", f)
        for f in d.rglob("*.fa"):
            if f.is_file() and f.stat().st_size > 0:
                if "protein" in f.name.lower():
                    if include_protein:
                        _add("protein_fasta", f)
                else:
                    _add("genome_fasta", f)
        for f in d.rglob("*.fasta"):
            if f.is_file() and f.stat().st_size > 0 and "protein" not in f.name.lower():
                _add("genome_fasta", f)
        for f in d.rglob("*.faa"):
            if f.is_file() and f.stat().st_size > 0 and include_protein:
                _add("protein_fasta", f)
        for f in d.rglob("*.gff*"):
            if f.is_file() and f.suffix in (".gff", ".gff3") and f.stat().st_size > 0:
                _add("gff", f)

    return out


def paths_from_manifest(manifest: dict, download_root: Path, label: str) -> dict[str, list[Path]]:
    """Convert manifest (paths as strings) to Path lists. Manifest keys: genome_fasta, gff, protein_fasta."""
    out: dict[str, list[Path]] = {"genome_fasta": [], "gff": [], "protein_fasta": []}
    for key in out:
        for p in manifest.get(key) or []:
            path = Path(p)
            if not path.is_absolute():
                path = download_root / label / path
            if path.exists():
                out[key].append(path)
    return out


def run_cmd(cmd: list[str], dry_run: bool, log: logging.Logger, cwd: Path | None = None) -> bool:
    log.info("Run: %s", " ".join(str(x) for x in cmd))
    if dry_run:
        return True
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
        return True
    except subprocess.CalledProcessError as e:
        log.error("Command failed: %s", e)
        return False
    except FileNotFoundError:
        log.error("Command not found: %s (install BLAST+)", cmd[0])
        return False


def atomic_replace_staging(staging: Path, out_dir: Path, log: logging.Logger) -> bool:
    """Atomically replace out_dir with staging. Returns True on success, False on failure."""
    if out_dir.exists():
        try:
            shutil.rmtree(out_dir)
        except OSError as e:
            log.error("Failed to remove existing DB dir %s: %s", out_dir, e)
            shutil.rmtree(staging, ignore_errors=True)
            return False
    try:
        staging.rename(out_dir)
        return True
    except OSError as e:
        log.error("Failed to rename staging to %s: %s", out_dir, e)
        shutil.rmtree(staging, ignore_errors=True)
        return False


def copy_metadata_fields(source: dict[str, Any], target: dict[str, Any], fields: list[str]) -> None:
    """Copy specified metadata fields from source to target if present and non-empty."""
    for field in fields:
        if field in source and source[field]:
            target[field] = source[field]


def concatenate_fasta(sources: list[Path], dest: Path, dry_run: bool) -> None:
    if dry_run:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as out:
        for i, f in enumerate(sorted(sources)):
            with open(f, "rb") as inp:
                data = inp.read()
            out.write(data)
            if not data.endswith(b"\n"):
                out.write(b"\n")


def build_nucl_db(
    combined_fasta: Path,
    gff_path: Path | None,
    blast_db_root: Path,
    label_name: str,
    max_file_sz: str | None,
    dry_run: bool,
    log: logging.Logger,
) -> tuple[bool, str | None, str | None]:
    """Build nucleotide DB. Returns (success, db_path, gff_path). Uses staging then atomic replace."""
    out_dir = blast_db_root / label_name / NUCL_DIR
    staging = out_dir.parent / (NUCL_DIR + STAGING_SUFFIX)
    db_basename = DB_BASENAME_NUCL
    db_path = out_dir / db_basename
    out_gff = out_dir / ANNOTATIONS_GFF
    if not dry_run:
        staging.mkdir(parents=True, exist_ok=True)
        fasta_in_staging = staging / SEQUENCES_NUCL
        shutil.copy2(combined_fasta, fasta_in_staging)
        if gff_path and gff_path.exists():
            shutil.copy2(gff_path, staging / ANNOTATIONS_GFF)
    else:
        fasta_in_staging = combined_fasta
    cmd = [
        MAKEBLASTDB,
        "-in", str(staging / SEQUENCES_NUCL) if not dry_run else str(combined_fasta),
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", str(staging / db_basename) if not dry_run else str(db_path),
        "-title", f"Genomes {label_name} nucleotide",
    ]
    if max_file_sz:
        cmd.extend(["-max_file_sz", max_file_sz])
    if not run_cmd(cmd, dry_run, log):
        return False, None, None
    if dry_run:
        return True, str(db_path), str(out_gff) if gff_path else None
    # Validate
    check_path = staging / db_basename
    val_cmd = [BLASTDBCMD, "-db", str(check_path), "-info"]
    if not run_cmd(val_cmd, False, log):
        shutil.rmtree(staging, ignore_errors=True)
        return False, None, None
    # Atomic replace
    if not atomic_replace_staging(staging, out_dir, log):
        return False, None, None
    return True, str(db_path), str(out_gff) if (out_dir / ANNOTATIONS_GFF).exists() else None


def build_prot_db(
    label_name: str,
    combined_fasta: Path,
    blast_db_root: Path,
    max_file_sz: str | None,
    dry_run: bool,
    log: logging.Logger,
) -> tuple[bool, str | None]:
    """Build protein DB. Returns (success, db_path). Uses staging then atomic replace."""
    out_dir = blast_db_root / label_name / PROT_DIR
    staging = out_dir.parent / (PROT_DIR + STAGING_SUFFIX)
    db_basename = DB_BASENAME_PROT
    db_path = out_dir / db_basename
    if not dry_run:
        staging.mkdir(parents=True, exist_ok=True)
        shutil.copy2(combined_fasta, staging / SEQUENCES_PROT)
    cmd = [
        MAKEBLASTDB,
        "-in", str(staging / SEQUENCES_PROT) if not dry_run else str(combined_fasta),
        "-dbtype", "prot",
        "-parse_seqids",
        "-out", str(staging / db_basename) if not dry_run else str(db_path),
        "-title", f"Genomes {label_name} protein",
    ]
    if max_file_sz:
        cmd.extend(["-max_file_sz", max_file_sz])
    if not run_cmd(cmd, dry_run, log):
        return False, None
    if dry_run:
        return True, str(db_path)
    check_path = staging / db_basename
    if not run_cmd([BLASTDBCMD, "-db", str(check_path), "-info"], False, log):
        shutil.rmtree(staging, ignore_errors=True)
        return False, None
    # Atomic replace
    if not atomic_replace_staging(staging, out_dir, log):
        return False, None
    return True, str(db_path)


METADATA_FIELDS = ["bioproject", "taxon", "taxon_id", "assembly_id", "assembly_level", "source", "genome_size", "reference"]


def build_one_label(
    label_name: str,
    download_root: Path,
    blast_db_root: Path,
    include_protein: bool,
    max_file_sz: str | None,
    dry_run: bool,
    log: logging.Logger,
    csv_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build nucl (and optionally prot) DB for one label. Returns dict with success, entries for registry, summary."""
    t0 = time.time()
    manifest = get_manifest_paths(download_root, label_name)
    # Collect metadata from manifest (preferred) or CSV fallback
    metadata = {}
    if manifest:
        paths = paths_from_manifest(manifest, download_root, label_name)
        copy_metadata_fields(manifest, metadata, METADATA_FIELDS)
    else:
        paths = discover_paths(download_root, label_name, include_protein)
    # CSV fallback for metadata if manifest missing fields
    if csv_metadata:
        copy_metadata_fields(csv_metadata, metadata, METADATA_FIELDS)
    if not paths["genome_fasta"]:
        log.warning("%s: no genome FASTA found", label_name)
        return {"success": False, "label": label_name, "nucl_ok": False, "prot_ok": False, "entries": [], "summary": f"{label_name} FAIL (no FASTA)"}
    # Combined FASTA for nucl
    with tempfile.NamedTemporaryFile(mode="wb", suffix=".fna", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        concatenate_fasta(paths["genome_fasta"], tmp_path, dry_run)
        gff_first = paths["gff"][0] if paths["gff"] else None
        nucl_ok, db_path_nucl, gff_path = build_nucl_db(
            tmp_path,
            gff_first,
            blast_db_root,
            label_name,
            max_file_sz,
            dry_run,
            log,
        )
    finally:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass  # Ignore cleanup errors
    entries = []
    if nucl_ok and db_path_nucl:
        entry = {
            "name": f"{label_name}_nucl",
            "db_path": db_path_nucl,
            "db_type": "nucl",
            "gff_path": gff_path,
            "blast_programs": ["blastn", "tblastn"],
        }
        # Add metadata if available
        copy_metadata_fields(metadata, entry, METADATA_FIELDS)
        entries.append(entry)
    prot_ok = False
    db_path_prot = None
    if include_protein and paths["protein_fasta"]:
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".faa", delete=False) as tmp:
            tmp_prot = Path(tmp.name)
        try:
            concatenate_fasta(paths["protein_fasta"], tmp_prot, dry_run)
            prot_ok, db_path_prot = build_prot_db(
                label_name, tmp_prot, blast_db_root, max_file_sz, dry_run, log
            )
            if prot_ok and db_path_prot:
                entry = {
                    "name": f"{label_name}_prot",
                    "db_path": db_path_prot,
                    "db_type": "prot",
                    "gff_path": gff_path,
                    "blast_programs": ["blastp", "blastx"],
                }
                # Add metadata if available
                copy_metadata_fields(metadata, entry, METADATA_FIELDS)
                entries.append(entry)
        finally:
            try:
                if tmp_prot.exists():
                    tmp_prot.unlink()
            except OSError:
                pass  # Ignore cleanup errors
    else:
        prot_ok = True  # skip is ok
    elapsed = time.time() - t0
    summary = f"{label_name} nucl={'OK' if nucl_ok else 'FAIL'}, prot={'OK' if prot_ok else 'FAIL'} ({elapsed:.1f}s)"
    log.info(summary)
    return {"success": nucl_ok, "label": label_name, "nucl_ok": nucl_ok, "prot_ok": prot_ok, "entries": entries, "summary": summary}


def build_all(
    csv_path: Path,
    config: dict[str, Any],
    dry_run: bool,
    log: logging.Logger,
    alias: bool,
) -> int:
    download_root = Path(config["download_root"])
    blast_db_root = Path(config["blast_db_root"])
    jobs = max(1, config.get("build_jobs", 1))
    max_file_sz = config.get("max_file_sz")
    rows = read_csv(csv_path)
    if not rows:
        log.warning("No rows in CSV")
        return 0
    blast_db_root.mkdir(parents=True, exist_ok=True)
    all_entries = []
    if jobs <= 1:
        for row in rows:
            label_name = row["label"]
            r = build_one_label(
                label_name,
                download_root,
                blast_db_root,
                row["include_protein"],
                max_file_sz,
                dry_run,
                log,
                row,  # Pass CSV row as metadata
            )
            all_entries.extend(r.get("entries") or [])
            print(r.get("summary", label_name))
    else:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futures = {
                ex.submit(
                    build_one_label,
                    row["label"],
                    download_root,
                    blast_db_root,
                    row["include_protein"],
                    max_file_sz,
                    dry_run,
                    log,
                    row,  # Pass CSV row as metadata
                ): row
                for row in rows
            }
            for fut in as_completed(futures):
                try:
                    r = fut.result()
                    all_entries.extend(r.get("entries") or [])
                    print(r.get("summary", ""))
                except Exception as e:
                    log.exception("%s", e)
                    return 1
    # Alias: create all_nucl (and optionally all_prot) if multiple labels
    if alias and len(rows) > 1 and not dry_run and all_entries:
        nucl_paths = [e["db_path"] for e in all_entries if e.get("db_type") == "nucl" and e.get("db_path")]
        prot_paths = [e["db_path"] for e in all_entries if e.get("db_type") == "prot" and e.get("db_path")]
        if nucl_paths:
            dblist = " ".join(nucl_paths)
            alias_cmd = [
                BLASTDB_ALIASTOOL,
                "-dblist", dblist,
                "-dbtype", "nucl",
                "-title", "All genomes nucleotide",
                "-out", str(blast_db_root / "all_nucl"),
            ]
            if run_cmd(alias_cmd, False, log):
                all_entries.append({
                    "name": "all_nucl",
                    "db_path": str(blast_db_root / "all_nucl"),
                    "db_type": "nucl",
                    "gff_path": None,
                    "blast_programs": ["blastn", "tblastn"],
                })
        if prot_paths:
            dblist = " ".join(prot_paths)
            alias_cmd = [
                BLASTDB_ALIASTOOL,
                "-dblist", dblist,
                "-dbtype", "prot",
                "-title", "All genomes protein",
                "-out", str(blast_db_root / "all_prot"),
            ]
            if run_cmd(alias_cmd, False, log):
                all_entries.append({
                    "name": "all_prot",
                    "db_path": str(blast_db_root / "all_prot"),
                    "db_type": "prot",
                    "gff_path": None,
                    "blast_programs": ["blastp", "blastx"],
                })
    # Write registry
    registry_path = blast_db_root / REGISTRY_NAME
    registry = {"databases": all_entries}
    log.info("Write registry: %s (%d entries)", registry_path, len(all_entries))
    if not dry_run and all_entries:
        try:
            with open(registry_path, "w", encoding="utf-8") as f:
                if yaml:
                    yaml.dump(registry, f, default_flow_style=False, sort_keys=False)
                else:
                    json.dump(registry, f, indent=2)
        except (OSError, yaml.YAMLError if yaml else None, json.JSONEncodeError) as e:
            log.error("Failed to write registry %s: %s", registry_path, e)
            return 1
        except (yaml.YAMLError if yaml else type(None), json.JSONEncodeError) as e:
            log.error("Failed to serialize registry %s: %s", registry_path, e)
            return 1
    return 0


def load_registry(registry_path: Path, log: logging.Logger) -> dict[str, Any] | None:
    """Load registry from file. Returns dict or None on error."""
    if not registry_path.exists():
        log.warning("No registry at %s", registry_path)
        return None
    try:
        with open(registry_path, "r", encoding="utf-8") as f:
            if yaml:
                data = yaml.safe_load(f)
            else:
                data = json.load(f)
        return data
    except OSError as e:
        log.error("Failed to read registry %s: %s", registry_path, e)
        return None
    except Exception as e:
        if yaml and isinstance(e, yaml.YAMLError):
            log.error("Failed to parse registry %s (YAML error): %s", registry_path, e)
        elif isinstance(e, json.JSONDecodeError):
            log.error("Failed to parse registry %s (JSON error): %s", registry_path, e)
        else:
            log.error("Failed to read registry %s: %s", registry_path, e)
        return None


def cmd_list(config: dict[str, Any], log: logging.Logger) -> int:
    blast_db_root = Path(config["blast_db_root"])
    registry_path = blast_db_root / REGISTRY_NAME
    data = load_registry(registry_path, log)
    if data is None:
        return 0
    for db in data.get("databases") or []:
        name = db.get("name", "")
        db_path = db.get("db_path", "")
        db_type = db.get("db_type", "")
        progs = db.get("blast_programs") or []
        log.info("%s\t%s\t%s\t%s", name, db_path, db_type, ",".join(progs))
        print(f"{name}\t{db_path}\t{db_type}\t{','.join(progs)}")
    return 0


def cmd_alias_only(config: dict[str, Any], dry_run: bool, log: logging.Logger) -> int:
    """Create alias databases from existing per-label databases in registry without rebuilding."""
    blast_db_root = Path(config["blast_db_root"])
    registry_path = blast_db_root / REGISTRY_NAME
    data = load_registry(registry_path, log)
    if data is None:
        log.error("Cannot create alias: no registry found. Run 'build' first.")
        return 1
    
    existing_entries = data.get("databases") or []
    # Filter out existing alias entries (all_nucl, all_prot) to avoid duplicates
    per_label_entries = [e for e in existing_entries if e.get("name") not in ("all_nucl", "all_prot")]
    
    if not per_label_entries:
        log.warning("No per-label databases found in registry")
        return 0
    
    # Collect nucl and prot DB paths from per-label entries
    nucl_paths = [e["db_path"] for e in per_label_entries if e.get("db_type") == "nucl" and e.get("db_path")]
    prot_paths = [e["db_path"] for e in per_label_entries if e.get("db_type") == "prot" and e.get("db_path")]
    
    if len(nucl_paths) < 2 and len(prot_paths) < 2:
        log.warning("Need at least 2 databases to create alias. Found %d nucl, %d prot", len(nucl_paths), len(prot_paths))
        return 0
    
    new_entries = list(per_label_entries)  # Start with existing per-label entries
    
    # Create all_nucl alias
    if nucl_paths and len(nucl_paths) >= 2:
        dblist = " ".join(nucl_paths)
        alias_cmd = [
            BLASTDB_ALIASTOOL,
            "-dblist", dblist,
            "-dbtype", "nucl",
            "-title", "All genomes nucleotide",
            "-out", str(blast_db_root / "all_nucl"),
        ]
        log.info("Creating all_nucl alias from %d databases", len(nucl_paths))
        if dry_run:
            log.info("Would run: %s", " ".join(alias_cmd))
            new_entries.append({
                "name": "all_nucl",
                "db_path": str(blast_db_root / "all_nucl"),
                "db_type": "nucl",
                "gff_path": None,
                "blast_programs": ["blastn", "tblastn"],
            })
        elif run_cmd(alias_cmd, False, log):
            new_entries.append({
                "name": "all_nucl",
                "db_path": str(blast_db_root / "all_nucl"),
                "db_type": "nucl",
                "gff_path": None,
                "blast_programs": ["blastn", "tblastn"],
            })
            log.info("Created all_nucl alias")
        else:
            log.error("Failed to create all_nucl alias")
            return 1
    
    # Create all_prot alias
    if prot_paths and len(prot_paths) >= 2:
        dblist = " ".join(prot_paths)
        alias_cmd = [
            BLASTDB_ALIASTOOL,
            "-dblist", dblist,
            "-dbtype", "prot",
            "-title", "All genomes protein",
            "-out", str(blast_db_root / "all_prot"),
        ]
        log.info("Creating all_prot alias from %d databases", len(prot_paths))
        if dry_run:
            log.info("Would run: %s", " ".join(alias_cmd))
            new_entries.append({
                "name": "all_prot",
                "db_path": str(blast_db_root / "all_prot"),
                "db_type": "prot",
                "gff_path": None,
                "blast_programs": ["blastp", "blastx"],
            })
        elif run_cmd(alias_cmd, False, log):
            new_entries.append({
                "name": "all_prot",
                "db_path": str(blast_db_root / "all_prot"),
                "db_type": "prot",
                "gff_path": None,
                "blast_programs": ["blastp", "blastx"],
            })
            log.info("Created all_prot alias")
        else:
            log.error("Failed to create all_prot alias")
            return 1
    
    # Update registry
    registry = {"databases": new_entries}
    log.info("Updating registry: %s (%d entries)", registry_path, len(new_entries))
    if not dry_run:
        try:
            with open(registry_path, "w", encoding="utf-8") as f:
                if yaml:
                    yaml.dump(registry, f, default_flow_style=False, sort_keys=False)
                else:
                    json.dump(registry, f, indent=2)
            log.info("Registry updated successfully")
        except (OSError, yaml.YAMLError if yaml else None, json.JSONDecodeError) as e:
            log.error("Failed to write registry %s: %s", registry_path, e)
            return 1
        except (yaml.YAMLError if yaml else type(None), json.JSONDecodeError) as e:
            log.error("Failed to serialize registry %s: %s", registry_path, e)
            return 1
    return 0


def setup_logging(log_path: Path | None, verbose: bool) -> logging.Logger:
    log = logging.getLogger("blast_db")
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    if log_path:
        fh = logging.FileHandler(log_path, encoding="utf-8")
        fh.setFormatter(fmt)
        log.addHandler(fh)
    ch = logging.StreamHandler(sys.stderr)
    ch.setFormatter(fmt)
    log.addHandler(ch)
    return log


def main() -> int:
    ap = argparse.ArgumentParser(description="Set up and manage BLAST databases from downloaded genomes")
    ap.add_argument("--csv", type=Path, default=Path(DEFAULT_CSV), help="Bioprojects CSV")
    ap.add_argument("--config", type=Path, default=Path("config.yml"), help="Config YAML")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--log-file", type=Path, default=Path(LOG_FILE), help="Log file")
    ap.add_argument("--no-log-file", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    sub = ap.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build", help="Build nucl (and optional prot) DBs per label")
    build_parser.add_argument("--alias", action="store_true", help="Create alias DB(s) combining all labels")
    update_parser = sub.add_parser("update", help="Rebuild DBs from current downloads")
    update_parser.add_argument("--alias", action="store_true", help="Create alias DB(s) combining all labels")
    sub.add_parser("alias", help="Create alias DB(s) from existing per-label databases (no rebuild)")
    sub.add_parser("list", help="List DBs in registry")
    args = ap.parse_args()
    config = load_config(args.config)
    log = setup_logging(None if args.no_log_file else args.log_file, args.verbose)
    if args.command == "list":
        return cmd_list(config, log)
    if args.command == "alias":
        return cmd_alias_only(config, args.dry_run, log)
    if args.command in ("build", "update"):
        alias_flag = getattr(args, "alias", False)
        return build_all(args.csv, config, args.dry_run, log, alias_flag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
