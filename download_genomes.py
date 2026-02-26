#!/usr/bin/env python3
"""
download_genomes.py — Download and manage genome assemblies, annotations, and
protein sequences for a set of BioProjects defined in a CSV manifest.

Example
-------
    python download_genomes.py --csv bioprojects.csv download
    python download_genomes.py --csv bioprojects.csv download --force
    python download_genomes.py --csv bioprojects.csv --dry-run download --dehydrated
    python download_genomes.py --csv bioprojects.csv status
    python download_genomes.py --csv bioprojects.csv repair

Subcommands
-----------
  status    Report per-row download state (not_started / dehydrated_only /
            rehydrated / complete).
  download  Fetch genome data for every non-skipped row. Writes a manifest.json
            in each label directory on success. Already-complete labels are
            skipped unless --force is given. Prints a summary table at the end.
  repair    Rehydrate any dehydrated datasets and re-validate (idempotent).

Workflow (download)
-------------------
  1. Read bioprojects.csv; skip rows where skip=true.
  2. For each row, check existing state:
     - If already complete and assembly_id hasn't changed, skip (unless --force).
     - If the manifest's assembly_id doesn't match the CSV, warn and skip
       (prevents silently mixing old and new data; use --force to override).
     - Remove any stale .zip from a prior interrupted run.
  3. Resolve the data source:
     a. If wormbase_species is set, try WormBase ParaSite first — download the
        genome FASTA, GFF3 annotations, and protein FASTA from the WormBase FTP
        (release controlled by wbps_version in config.yml, default WBPS19).
        Checks that all expected files are present, not just that the directory
        is non-empty (handles partial downloads from prior interrupted runs).
     b. If WormBase provides a complete set, use it as the primary source.
     c. Otherwise fall back to NCBI Datasets CLI (ncbi-datasets-cli). If NCBI
        fails and an assembly_id is available, attempt an ENA download via
        enaBrowserTools, then NCBI GenBank FTP.
     d. After NCBI/ENA download, if wormbase_species is set and WormBase data
        was not already fetched, supplement with WormBase annotations/proteins.
  4. Extract and (if dehydrated) rehydrate NCBI ZIP archives.
  5. Validate that required files are present (genome FASTA at minimum).
  6. Write downloads/<label>/manifest.json with discovered file paths and
     metadata (taxon, assembly_id, source, etc.).
  7. Warn about orphaned download directories not in the current CSV.
  8. Print a per-row summary table to stdout.

Data source detection is automatic: the bioproject prefix determines the
archive of origin (PRJNA → NCBI, PRJEB/PRJEA → ENA, PRJDB → DDBJ). No
explicit "source" column is needed in the CSV.

Dependencies: ncbi-datasets-cli, enaBrowserTools (optional), PyYAML (optional).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import re
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from shutil import which
from typing import Any
from urllib.request import urlopen, urlretrieve
from urllib.error import URLError, HTTPError

try:
    import yaml
except ImportError:
    yaml = None

LOG_FILE = "download.log"
DEFAULT_DOWNLOAD_ROOT = "downloads"
DEFAULT_CSV = "bioprojects.csv"
DATASETS_CMD = "datasets"
ENA_DATA_GET = "enaDataGet"
NCBI_DATA_DIR = "ncbi_dataset"
DATA_SUBDIR = "data"
MANIFEST_NAME = "manifest.json"
FETCH_TXT = "fetch.txt"
ENA_FTP_BASE = "ftp://ftp.ebi.ac.uk/pub/databases/ena/assembly"
ENA_FTP_ASSEMBLY_BASE = "https://www.ebi.ac.uk/ena/portal/api/filereport"
WBPS_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases"
WBPS_DEFAULT_VERSION = "WBPS19"
WBPS_DATA_DIR = "wormbase_parasite"

WBPS_FILE_TYPES = {
    "genome": "genomic.fa",
    "genome_softmasked": "genomic_softmasked.fa",
    "genome_masked": "genomic_masked.fa",
    "gff3": "annotations.gff3",
    "protein": "protein.fa",
    "cds_transcripts": "CDS_transcripts.fa",
    "mrna_transcripts": "mRNA_transcripts.fa",
    "gtf": "canonical_geneset.gtf",
}


def load_config(config_path: Path | None) -> dict[str, Any]:
    """Load config.yml if present. Returns dict with download_root, rehydrate_workers, etc."""
    out = {
        "download_root": DEFAULT_DOWNLOAD_ROOT,
        "rehydrate_workers": None,
        "assembly_level": None,
        "wbps_version": WBPS_DEFAULT_VERSION,
    }
    if config_path is None or not config_path.exists():
        return out
    if yaml is None:
        logging.warning("PyYAML not installed; config.yml ignored")
        return out
    try:
        with open(config_path, "r") as f:
            data = yaml.safe_load(f) or {}
        out["download_root"] = data.get("download_root", out["download_root"])
        out["rehydrate_workers"] = data.get("rehydrate_workers")
        out["assembly_level"] = data.get("assembly_level")
        out["wbps_version"] = data.get("wbps_version", out["wbps_version"])
    except Exception as e:
        logging.warning("Failed to load config.yml: %s; using defaults", e)
    return out


def wbps_file_url(species: str, bioproject: str, filetype: str, version: str = WBPS_DEFAULT_VERSION) -> str:
    """Construct WormBase ParaSite FTP URL for a specific file type.
    ``species`` is lowercase underscore-separated (e.g. ``hymenolepis_microstoma``).
    ``filetype`` is the suffix before ``.gz`` (e.g. ``annotations.gff3``, ``protein.fa``).
    """
    return (
        f"{WBPS_FTP_BASE}/{version}/species/{species}/{bioproject}/"
        f"{species}.{bioproject}.{version}.{filetype}.gz"
    )



def find_wormbase_data_dir(label_dir: Path) -> Path | None:
    """Return the wormbase_parasite subdirectory if it exists."""
    d = label_dir / WBPS_DATA_DIR
    return d if d.is_dir() else None


def _decompress_gz(gz_path: Path, log: logging.Logger) -> Path:
    """Decompress a .gz file in place, returning the path to the decompressed file."""
    out_path = gz_path.with_suffix("")  # strip .gz
    log.debug("Decompressing %s -> %s", gz_path.name, out_path.name)
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()
    return out_path


def download_from_wormbase(
    wormbase_species: str,
    bioproject: str,
    label_dir: Path,
    file_types: list[str],
    dry_run: bool,
    log: logging.Logger,
    wbps_version: str = WBPS_DEFAULT_VERSION,
) -> dict[str, Path]:
    """Download specified file types from WormBase ParaSite FTP.

    ``file_types`` are keys from WBPS_FILE_TYPES (e.g. ``["gff3", "protein"]``).
    Uses the row's own bioproject as the WormBase bioproject (annotations must match the assembly).
    Returns a dict mapping successfully downloaded file type keys to their local paths.
    """
    wbps_dir = label_dir / WBPS_DATA_DIR
    log.info("Attempting WormBase ParaSite download for %s/%s (version %s)",
             wormbase_species, bioproject, wbps_version)

    if dry_run:
        for ft in file_types:
            suffix = WBPS_FILE_TYPES.get(ft)
            if suffix:
                url = wbps_file_url(wormbase_species, bioproject, suffix, wbps_version)
                log.info("Would download: %s", url)
        return {}

    wbps_dir.mkdir(parents=True, exist_ok=True)
    downloaded: dict[str, Path] = {}

    for ft in file_types:
        suffix = WBPS_FILE_TYPES.get(ft)
        if suffix is None:
            log.warning("Unknown WormBase file type: %s (valid: %s)", ft, ", ".join(WBPS_FILE_TYPES))
            continue
        url = wbps_file_url(wormbase_species, bioproject, suffix, wbps_version)
        gz_filename = f"{wormbase_species}.{bioproject}.{wbps_version}.{suffix}.gz"
        gz_path = wbps_dir / gz_filename
        log.info("Downloading %s from WormBase ParaSite...", suffix)
        log.debug("URL: %s", url)
        try:
            urlretrieve(url, gz_path)
            if gz_path.exists() and gz_path.stat().st_size > 0:
                out_path = _decompress_gz(gz_path, log)
                downloaded[ft] = out_path
                log.info("Downloaded and decompressed: %s (%d bytes)", out_path.name, out_path.stat().st_size)
            else:
                log.warning("Downloaded file is empty: %s", gz_path)
                gz_path.unlink(missing_ok=True)
        except (URLError, HTTPError, OSError) as e:
            log.warning("Failed to download %s from WormBase: %s (URL: %s)", suffix, e, url)
            gz_path.unlink(missing_ok=True)

    if downloaded:
        log.info("WormBase ParaSite: downloaded %d/%d requested file(s)", len(downloaded), len(file_types))
    else:
        log.warning("WormBase ParaSite: no files downloaded for %s/%s", wormbase_species, bioproject)
    return downloaded


def read_csv(csv_path: Path) -> list[dict[str, str]]:
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
            r["assembly_level"] = (r.get("assembly_level") or "").strip() or None
            inc_ann = (r.get("include_annotation") or "true").strip().lower()
            r["include_annotation"] = inc_ann in ("true", "1", "yes")
            inc_prot = (r.get("include_protein") or "true").strip().lower()
            r["include_protein"] = inc_prot in ("true", "1", "yes")
            r["wormbase_species"] = (r.get("wormbase_species") or "").strip()
            rows.append(r)
    return rows


def get_label_dir(download_root: Path, label: str) -> Path:
    return download_root / label


def has_fetch_txt(label_dir: Path) -> bool:
    return (label_dir / FETCH_TXT).exists()


def detect_project_source(bioproject: str) -> str:
    """Detect project source from BioProject prefix. Returns 'ENA', 'DDBJ', or 'NCBI'."""
    bioproject_upper = bioproject.upper()
    if bioproject_upper.startswith(("PRJEB", "PRJEA")):
        return "ENA"
    elif bioproject_upper.startswith("PRJDB"):
        return "DDBJ"
    else:
        return "NCBI"


def find_ncbi_data_dir(label_dir: Path) -> Path | None:
    d = label_dir / NCBI_DATA_DIR / DATA_SUBDIR
    return d if d.is_dir() else None


def find_ena_data_dir(label_dir: Path) -> Path | None:
    """Find ENA download directory. enaDataGet creates a dir named after the accession."""
    # enaDataGet creates a directory named after the accession (e.g., GCA_000469805.3)
    # Look for common ENA directory patterns
    for item in label_dir.iterdir():
        if item.is_dir() and (item.name.startswith("GCA_") or item.name.startswith("GCF_")):
            # Check if it contains sequence files
            if any(item.rglob("*.fa*")) or any(item.rglob("*.embl")):
                return item
    return None




def _collect_genome_fasta(dirs: list[Path | None]) -> list[Path]:
    """Collect genome FASTA files from NCBI, ENA, and WormBase data directories."""
    fnas: list[Path] = []
    for d in dirs:
        if d is None:
            continue
        fnas.extend(list(d.rglob("*.fna")))
        fnas.extend([f for f in d.rglob("*.fa") if "protein" not in f.name.lower()])
        fnas.extend([f for f in d.rglob("*.fasta") if "protein" not in f.name.lower()])
    return fnas


def _collect_gff(dirs: list[Path | None]) -> list[Path]:
    """Collect GFF/GFF3 files from all data directories."""
    gffs: list[Path] = []
    for d in dirs:
        if d is None:
            continue
        gffs.extend(f for f in d.rglob("*.gff*") if f.suffix in (".gff", ".gff3"))
    return gffs


def _collect_protein(dirs: list[Path | None]) -> list[Path]:
    """Collect protein FASTA files from all data directories."""
    faas: list[Path] = []
    for d in dirs:
        if d is None:
            continue
        faas.extend(list(d.rglob("*.faa")))
        faas.extend([f for f in d.rglob("*.fa") if "protein" in f.name.lower()])
    return faas


def get_status(label_dir: Path, include_annotation: bool, include_protein: bool) -> dict[str, Any]:
    """Check download state for a label directory.

    Returns a dict with keys:
        status:    not_started | dehydrated | incomplete | complete
        has_genome, has_gff, has_protein:  bool
        needs_rehydration:  bool  (fetch.txt present)
        sources:  list[str]  e.g. ["ncbi", "wormbase"]
        missing:  list[str]  human-readable list of what's absent
    """
    result: dict[str, Any] = {
        "status": "not_started",
        "has_genome": False,
        "has_gff": False,
        "has_protein": False,
        "needs_rehydration": False,
        "sources": [],
        "missing": [],
    }
    if not label_dir.exists():
        return result

    data_dir = find_ncbi_data_dir(label_dir)
    ena_dir = find_ena_data_dir(label_dir)
    wbps_dir = find_wormbase_data_dir(label_dir)
    if data_dir:
        result["sources"].append("ncbi")
    if ena_dir:
        result["sources"].append("ena")
    if wbps_dir:
        result["sources"].append("wormbase")

    result["needs_rehydration"] = has_fetch_txt(label_dir)

    dirs: list[Path | None] = [data_dir, ena_dir, wbps_dir]
    fnas = _collect_genome_fasta(dirs)
    if fnas and all(f.stat().st_size > 0 for f in fnas):
        result["has_genome"] = True
    gffs = _collect_gff(dirs)
    if gffs and all(f.stat().st_size > 0 for f in gffs):
        result["has_gff"] = True
    faas = _collect_protein(dirs)
    if faas and all(f.stat().st_size > 0 for f in faas):
        result["has_protein"] = True

    if not result["has_genome"]:
        result["missing"].append("genome")
        if result["needs_rehydration"]:
            result["status"] = "dehydrated"
        elif result["sources"]:
            result["status"] = "incomplete"
        return result

    if include_annotation and not result["has_gff"]:
        result["missing"].append("gff3")
    if include_protein and not result["has_protein"]:
        result["missing"].append("protein")

    if result["missing"]:
        result["status"] = "incomplete"
    else:
        result["status"] = "complete"
    return result


def run_cmd(cmd: list[str], dry_run: bool, log: logging.Logger, max_retries: int = 3, retry_delay: int = 5) -> bool:
    """
    Run command with retry logic for network failures.
    Returns True on success, False on failure.
    """
    log.info("Run: %s", " ".join(cmd))
    if dry_run:
        return True
    for attempt in range(1, max_retries + 1):
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                return True
            # Command failed - check if it's a network error
            stderr_text = (result.stderr or "").lower()
            stdout_text = (result.stdout or "").lower()
            combined = stderr_text + stdout_text
            is_network_error = any(term in combined for term in ["timeout", "i/o timeout", "connection", "network", "gateway", "dial tcp"])
            if is_network_error and attempt < max_retries:
                log.warning("Network error (attempt %d/%d). Retrying in %d seconds...", attempt, max_retries, retry_delay)
                if stderr_text:
                    log.debug("Error details: %s", stderr_text[:200])
                time.sleep(retry_delay)
                continue
            # Not a retryable error or out of retries
            error_msg = f"Command returned {result.returncode}"
            if stderr_text:
                error_msg += f"\nstderr: {stderr_text[:500]}"
            log.error("Command failed: %s", error_msg)
            return False
        except FileNotFoundError:
            log.error("Command not found: %s (install ncbi-datasets-cli)", cmd[0])
            return False
    return False


def build_include(include_annotation: bool, include_protein: bool) -> str:
    parts = ["genome"]
    if include_annotation:
        parts.append("gff3")
    if include_protein:
        parts.append("protein")
    return ",".join(parts)


def _download_url(url: str, dest: Path, log: logging.Logger) -> bool:
    """Download a single URL to dest. Returns True on success."""
    try:
        urlretrieve(url, dest)
        if dest.exists() and dest.stat().st_size > 0:
            log.info("Downloaded: %s", dest)
            return True
        dest.unlink(missing_ok=True)
    except (URLError, OSError) as e:
        log.debug("Download failed for %s: %s", url, e)
        dest.unlink(missing_ok=True)
    return False


def _parse_ena_filereport(content: str, log: logging.Logger) -> list[str]:
    """Extract FTP/HTTP URLs from an ENA Portal API filereport response."""
    lines = content.strip().split("\n")
    if len(lines) < 2 or not lines[1].strip():
        return []
    headers = lines[0].split("\t")
    data = lines[1].split("\t")
    urls: list[str] = []
    for col in ("generated_ftp", "submitted_ftp", "generated_bytes", "submitted_bytes"):
        try:
            idx = headers.index(col)
        except ValueError:
            continue
        if "bytes" in col:
            continue
        if idx < len(data) and data[idx].strip():
            for u in data[idx].strip().split(";"):
                u = u.strip()
                if u:
                    if not u.startswith(("ftp://", "http://", "https://")):
                        u = "https://" + u
                    urls.append(u)
    if not urls:
        log.debug("ENA filereport columns: %s", ", ".join(headers))
    return urls


def _ncbi_ftp_url(assembly_id: str) -> str:
    """Construct the NCBI GenBank FTP directory URL for a GCA/GCF accession.

    Path structure: ftp.ncbi.nlm.nih.gov/genomes/all/GCA/021/556/725/GCA_021556725.1_<name>/
    We can't predict <name>, so we return the parent directory for listing.
    """
    parts = assembly_id.split("_", 1)
    prefix = parts[0]
    number = parts[1].split(".")[0] if len(parts) > 1 else ""
    # Pad to 9 digits and split into 3-digit groups
    number = number.zfill(9)
    chunks = [number[i:i+3] for i in range(0, 9, 3)]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{'/'.join(chunks)}/"


def download_from_ena_ftp(
    assembly_id: str,
    label_dir: Path,
    log: logging.Logger,
) -> bool:
    """Download assembly FASTA from ENA Portal API or direct FTP.
    Returns True on success, False on failure.
    """
    log.info("Attempting direct FTP download for %s", assembly_id)
    label_dir.mkdir(parents=True, exist_ok=True)

    # --- ENA Portal API (correct field names: generated_ftp, submitted_ftp) ---
    api_url = (
        f"{ENA_FTP_ASSEMBLY_BASE}?accession={assembly_id}"
        f"&result=assembly&fields=generated_ftp,submitted_ftp"
    )
    try:
        log.debug("Querying ENA API: %s", api_url)
        with urlopen(api_url, timeout=30) as response:
            content = response.read().decode("utf-8")
        urls = _parse_ena_filereport(content, log)
        if urls:
            log.debug("ENA API returned %d URL(s) for %s", len(urls), assembly_id)
            for url in urls:
                filename = url.split("/")[-1]
                if _download_url(url, label_dir / filename, log):
                    return True
        else:
            log.debug("ENA API returned no FTP URLs for %s", assembly_id)
    except HTTPError as e:
        body = ""
        try:
            body = e.read().decode("utf-8")[:300]
        except Exception:
            pass
        log.debug("ENA API HTTP %d for %s: %s", e.code, assembly_id, body)
    except (URLError, TimeoutError, ValueError) as e:
        log.debug("ENA API query failed for %s: %s", assembly_id, e)

    # --- NCBI GenBank FTP (works even when datasets CLI rejects the accession) ---
    ncbi_dir_url = _ncbi_ftp_url(assembly_id)
    log.debug("Trying NCBI GenBank FTP listing: %s", ncbi_dir_url)
    try:
        with urlopen(ncbi_dir_url, timeout=30) as response:
            listing = response.read().decode("utf-8")
        # HTML directory listing — look for subdirectory matching the assembly_id
        # The actual files live inside a subdirectory named like GCA_021556725.1_ASM2155672v1/
        # Match href="GCA_021556725.1_..." links
        pattern = re.compile(rf'href="({re.escape(assembly_id)}[^"/]*)/?"', re.IGNORECASE)
        matches = pattern.findall(listing)
        for subdir in matches:
            for suffix in (f"{subdir}_genomic.fna.gz", f"{subdir}_protein.faa.gz",
                           f"{subdir}_genomic.gff.gz"):
                file_url = f"{ncbi_dir_url}{subdir}/{suffix}"
                dest = label_dir / suffix
                log.info("Trying NCBI FTP: %s", file_url)
                _download_url(file_url, dest, log)
            # Check if we got at least a genome FASTA
            if any(label_dir.rglob("*.fna*")):
                log.info("Downloaded from NCBI GenBank FTP for %s", assembly_id)
                return True
    except (HTTPError, URLError, TimeoutError) as e:
        log.debug("NCBI GenBank FTP failed for %s: %s", assembly_id, e)

    log.warning("FTP download failed for %s (tried ENA API and NCBI GenBank FTP)", assembly_id)
    log.info("Check manually: https://www.ebi.ac.uk/ena/browser/view/%s", assembly_id)
    log.info("Or: https://www.ncbi.nlm.nih.gov/datasets/genome/%s/", assembly_id)
    return False


def download_from_ena(
    assembly_id: str,
    label_dir: Path,
    include_annotation: bool,
    include_protein: bool,
    dry_run: bool,
    log: logging.Logger,
) -> bool:
    """
    Download assembly from ENA using enaBrowserTools (enaDataGet), with FTP fallback.
    Returns True on success, False on failure.
    """
    log.info("Attempting ENA download for %s", assembly_id)
    if dry_run:
        log.info("Would try: %s -f fasta -d %s %s (or direct FTP)", ENA_DATA_GET, label_dir, assembly_id)
        return True
    
    # Try enaDataGet first if available
    ena_cmd_path = which(ENA_DATA_GET)
    if ena_cmd_path:
        log.debug("Found %s at: %s", ENA_DATA_GET, ena_cmd_path)
        label_dir.mkdir(parents=True, exist_ok=True)
        cmd = [ENA_DATA_GET, "-f", "fasta", "-d", str(label_dir), assembly_id]
        if run_cmd(cmd, dry_run, log):
            # Check if download succeeded
            ena_dir = find_ena_data_dir(label_dir)
            if ena_dir and any(ena_dir.rglob("*.fa*")):
                log.info("ENA download succeeded via enaDataGet: %s", ena_dir)
                return True
        log.info("enaDataGet failed or incomplete, trying direct FTP download...")
    
    # Fallback to direct FTP download
    if download_from_ena_ftp(assembly_id, label_dir, log):
        # Check if download succeeded
        if any(label_dir.rglob("*.fa*")):
            log.info("ENA download succeeded via direct FTP")
            return True
    
    log.warning("ENA download failed for %s (tried both enaDataGet and direct FTP)", assembly_id)
    log.info("You may need to download manually from ENA FTP or fix enaBrowserTools Python compatibility")
    return False


def discover_manifest_paths(label_dir: Path, include_protein: bool) -> dict[str, list[str]]:
    """Discover paths to genome FASTA, GFF, and optionally protein FASTA.
    Supports NCBI, ENA, and WormBase ParaSite download structures."""
    out: dict[str, list[str]] = {"genome_fasta": [], "gff": [], "protein_fasta": []}
    seen: set[str] = set()

    def _add(key: str, path: Path) -> None:
        resolved = str(path.resolve())
        if resolved not in seen:
            seen.add(resolved)
            out[key].append(resolved)

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


def write_manifest(label_dir: Path, manifest_data: dict, dry_run: bool, log: logging.Logger) -> None:
    path = label_dir / MANIFEST_NAME
    log.info("Write manifest: %s", path)
    if dry_run:
        return
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(manifest_data, f, indent=2)
    except (OSError, json.JSONEncodeError) as e:
        log.error("Failed to write manifest %s: %s", path, e)
        raise


def validate_download(label_dir: Path, include_annotation: bool, include_protein: bool, log: logging.Logger) -> tuple[bool, bool, bool]:
    """
    Check expected files exist and have size > 0. Returns (valid, missing_gff, missing_protein).
    Genome FASTA is required; GFF and protein are optional (warn if missing but requested).
    Supports NCBI, ENA, and WormBase ParaSite download structures.
    """
    dirs: list[Path | None] = [find_ncbi_data_dir(label_dir), find_ena_data_dir(label_dir), find_wormbase_data_dir(label_dir)]
    if all(d is None for d in dirs):
        log.warning("No ncbi_dataset/data, ENA, or WormBase download dir under %s", label_dir)
        return False, False, False
    fnas = _collect_genome_fasta(dirs)
    if not fnas:
        log.warning("No genome FASTA files under %s", label_dir)
        return False, False, False
    for f in fnas:
        if f.stat().st_size == 0:
            log.warning("Empty file: %s", f)
            return False, False, False
    missing_gff = False
    if include_annotation:
        gffs = _collect_gff(dirs)
        if not gffs:
            missing_gff = True
            log.warning("No .gff/.gff3 under %s (requested but not available for this assembly)", label_dir)
        else:
            for f in gffs:
                if f.stat().st_size == 0:
                    log.warning("Empty file: %s", f)
    missing_protein = False
    if include_protein:
        faas = _collect_protein(dirs)
        if not faas:
            missing_protein = True
            log.warning("No protein FASTA under %s (requested but not available for this assembly)", label_dir)
        else:
            for f in faas:
                if f.stat().st_size == 0:
                    log.warning("Empty file: %s", f)
    return True, missing_gff, missing_protein  # Valid if genome FASTA exists and is non-empty


def _yn(val: bool) -> str:
    return "yes" if val else "-"


def cmd_status(
    csv_path: Path,
    config: dict[str, Any],
    dry_run: bool,
    log: logging.Logger,
) -> int:
    download_root = Path(config["download_root"])
    rows = read_csv(csv_path)
    header = f"{'label':<16}{'bioproject':<16}{'status':<14}{'genome':<8}{'gff3':<8}{'protein':<8}{'sources'}"
    print(header)
    print("-" * len(header))
    for row in rows:
        label = row["label"]
        label_dir = get_label_dir(download_root, label)
        st = get_status(label_dir, row["include_annotation"], row["include_protein"])
        status_str = st["status"]
        if st["missing"] and st["status"] == "incomplete":
            status_str = f"incomplete"
        sources = ",".join(st["sources"]) if st["sources"] else "-"
        line = (f"{label:<16}{row['bioproject']:<16}{status_str:<14}"
                f"{_yn(st['has_genome']):<8}{_yn(st['has_gff']):<8}{_yn(st['has_protein']):<8}"
                f"{sources}")
        print(line)
        log.debug("%s", line)
    # Orphaned directories
    active = {row["label"] for row in rows}
    if download_root.is_dir():
        for child in sorted(download_root.iterdir()):
            if child.is_dir() and child.name not in active:
                print(f"{child.name:<16}{'':<16}{'orphaned':<14}{'-':<8}{'-':<8}{'-':<8}-")
    return 0


def _read_manifest(label_dir: Path) -> dict | None:
    """Read an existing manifest.json if present."""
    path = label_dir / MANIFEST_NAME
    if not path.exists():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def _wbps_expected_files(wbps_dir: Path, species: str, bioproject: str,
                         version: str, file_types: list[str]) -> list[str]:
    """Return basenames expected in the WormBase directory for the given file types."""
    names: list[str] = []
    for ft in file_types:
        suffix = WBPS_FILE_TYPES.get(ft)
        if suffix:
            names.append(f"{species}.{bioproject}.{version}.{suffix}")
    return names


def _wbps_download_complete(wbps_dir: Path | None, species: str, bioproject: str,
                            version: str, file_types: list[str]) -> bool:
    """Check whether all expected WormBase files are present and non-empty."""
    if wbps_dir is None or not wbps_dir.is_dir():
        return False
    for name in _wbps_expected_files(wbps_dir, species, bioproject, version, file_types):
        f = wbps_dir / name
        if not f.exists() or f.stat().st_size == 0:
            return False
    return True


def cmd_download(
    csv_path: Path,
    config: dict[str, Any],
    dry_run: bool,
    log: logging.Logger,
    use_dehydrated: bool = False,
    force: bool = False,
) -> int:
    download_root = Path(config["download_root"])
    rows = read_csv(csv_path)
    rehydrate_workers = config.get("rehydrate_workers")
    assembly_level = config.get("assembly_level")
    active_labels: set[str] = set()
    summaries: list[str] = []

    for row in rows:
        label = row["label"]
        bioproject = row["bioproject"]
        label_dir = get_label_dir(download_root, label)
        assembly_id = row.get("assembly_id", "").strip()
        wb_species = row.get("wormbase_species", "")
        wbps_version = config.get("wbps_version", WBPS_DEFAULT_VERSION)
        skip_wb = config.get("skip_wormbase", False)
        project_source = detect_project_source(bioproject)
        active_labels.add(label)

        # --- Check existing state before downloading ---
        st = get_status(label_dir, row["include_annotation"], row["include_protein"])

        if st["status"] == "complete" and not force:
            existing = _read_manifest(label_dir)
            old_asm = (existing.get("assembly_id") or "").strip() if existing else ""
            if assembly_id and old_asm and assembly_id != old_asm:
                log.warning("%s: assembly changed (%s -> %s); re-run with --force to replace",
                            label, old_asm, assembly_id)
                summaries.append(f"{label}\t{bioproject}\tSTALE (assembly mismatch: {old_asm} != {assembly_id})")
                continue
            log.info("%s: already complete; skipping (use --force to re-download)", label)
            summaries.append(f"{label}\t{bioproject}\tskipped (complete)")
            continue

        # --- Clean up stale ZIP from a previous interrupted run ---
        zip_path = download_root / f"{label}.zip"
        if zip_path.exists():
            log.warning("%s: removing stale ZIP from previous run: %s", label, zip_path)
            if not dry_run:
                zip_path.unlink()

        download_root.mkdir(parents=True, exist_ok=True)
        wormbase_primary = False

        # --- WormBase ParaSite as primary source when available ---
        if wb_species and not skip_wb:
            wbps_types = ["genome"]
            if row["include_annotation"]:
                wbps_types.append("gff3")
            if row["include_protein"]:
                wbps_types.append("protein")
            wbps_dir = find_wormbase_data_dir(label_dir)
            if _wbps_download_complete(wbps_dir, wb_species, bioproject, wbps_version, wbps_types):
                log.info("%s: WormBase files already present", label)
            else:
                log.info("%s: using WormBase ParaSite as primary source for: %s", label, ", ".join(wbps_types))
                label_dir.mkdir(parents=True, exist_ok=True)
                download_from_wormbase(wb_species, bioproject, label_dir, wbps_types, dry_run, log, wbps_version)
            if not dry_run:
                valid, _, _ = validate_download(label_dir, row["include_annotation"], row["include_protein"], log)
                if valid:
                    wormbase_primary = True

        # --- NCBI / ENA fallback (when WormBase unavailable or incomplete) ---
        if not wormbase_primary:
            include = build_include(row["include_annotation"], row["include_protein"])
            zip_name = f"{label}.zip"
            zip_path = download_root / zip_name
            accession = assembly_id or bioproject
            cmd = [
                DATASETS_CMD, "download", "genome", "accession", accession,
                "--filename", str(zip_path), "--include", include,
            ]
            row_assembly = row.get("assembly_level") or assembly_level
            if row_assembly:
                cmd.extend(["--assembly-level", row_assembly])
            if use_dehydrated:
                cmd.append("--dehydrated")
            ncbi_success = run_cmd(cmd, dry_run, log)
            ena_downloaded = False

            if not ncbi_success and assembly_id:
                log.info("%s: NCBI download failed, trying ENA download for %s", label, assembly_id)
                if download_from_ena(assembly_id, label_dir, row["include_annotation"], row["include_protein"], dry_run, log):
                    log.info("%s: Successfully downloaded from ENA", label)
                    ena_downloaded = True
                else:
                    log.error("Failed to download %s (%s) from all sources; continuing with next row", label, bioproject)
                    summaries.append(f"{label}\t{bioproject}\tFAILED (all sources)")
                    continue
            elif not ncbi_success:
                log.error("Failed to download %s (%s); continuing with next row", label, bioproject)
                summaries.append(f"{label}\t{bioproject}\tFAILED")
                continue

            if not ena_downloaded and not dry_run and zip_path.exists():
                label_dir.mkdir(parents=True, exist_ok=True)
                try:
                    with zipfile.ZipFile(zip_path, "r") as z:
                        z.extractall(label_dir)
                    zip_path.unlink()
                except (zipfile.BadZipFile, zipfile.LargeZipFile, OSError) as e:
                    log.error("Failed to extract %s: %s", zip_path, e)
                    return 1
            if not ena_downloaded and not dry_run and label_dir.exists() and has_fetch_txt(label_dir):
                rehydrate_cmd = [DATASETS_CMD, "rehydrate", "--directory", str(label_dir)]
                if rehydrate_workers is not None:
                    rehydrate_cmd.extend(["--max-workers", str(rehydrate_workers)])
                if not run_cmd(rehydrate_cmd, dry_run=False, log=log):
                    return 1

            # Supplement with WormBase annotations/protein after NCBI/ENA
            if wb_species and not skip_wb and not dry_run:
                wbps_types_supp: list[str] = []
                if row["include_annotation"]:
                    wbps_types_supp.append("gff3")
                if row["include_protein"]:
                    wbps_types_supp.append("protein")
                if wbps_types_supp:
                    wbps_dir = find_wormbase_data_dir(label_dir)
                    if not _wbps_download_complete(wbps_dir, wb_species, bioproject, wbps_version, wbps_types_supp):
                        log.info("%s: supplementing NCBI/ENA download with WormBase annotations", label)
                        download_from_wormbase(wb_species, bioproject, label_dir, wbps_types_supp, dry_run, log, wbps_version)

        # --- Validate and write manifest ---
        if not dry_run and label_dir.exists():
            valid, missing_gff, missing_protein = validate_download(label_dir, row["include_annotation"], row["include_protein"], log)
            if valid:
                manifest_data = discover_manifest_paths(label_dir, row["include_protein"])
                manifest_data["bioproject"] = bioproject
                manifest_data["label"] = label
                metadata_fields = ["taxon", "taxon_id", "assembly_id", "assembly_level", "source", "genome_size", "reference",
                                   "wormbase_species"]
                for field in metadata_fields:
                    if field in row and row[field]:
                        manifest_data[field] = row[field]
                missing_parts = []
                if missing_gff:
                    missing_parts.append("gff3")
                if missing_protein:
                    missing_parts.append("protein")
                if missing_parts:
                    if wb_species:
                        log.info("%s: some files still missing after WormBase download", label)
                    elif project_source in ("ENA", "DDBJ"):
                        log.info("%s: missing files may be available from WormBase ParaSite or ENA FTP", label)
                write_manifest(label_dir, manifest_data, dry_run=False, log=log)
                if missing_parts:
                    summaries.append(f"{label}\t{bioproject}\tdownloaded (missing {', '.join(missing_parts)})")
                else:
                    summaries.append(f"{label}\t{bioproject}\tdownloaded (complete)")
            else:
                log.warning("%s: validation failed (missing required genome FASTA); skipping manifest", label)
                summaries.append(f"{label}\t{bioproject}\tFAILED (no genome FASTA)")
        elif dry_run:
            summaries.append(f"{label}\t{bioproject}\tdry-run")

    # --- Warn about orphaned directories ---
    if download_root.is_dir():
        for child in sorted(download_root.iterdir()):
            if child.is_dir() and child.name not in active_labels:
                log.warning("Orphaned download directory (not in CSV): %s", child)

    # --- Print summary ---
    if summaries:
        print()
        print("label\tbioproject\tresult")
        for s in summaries:
            print(s)

    return 0


def cmd_repair(
    csv_path: Path,
    config: dict[str, Any],
    dry_run: bool,
    log: logging.Logger,
) -> int:
    """Rehydrate dehydrated downloads and supplement with WormBase where needed."""
    download_root = Path(config["download_root"])
    rows = read_csv(csv_path)
    rehydrate_workers = config.get("rehydrate_workers")
    summaries: list[str] = []

    for row in rows:
        label = row["label"]
        bioproject = row["bioproject"]
        label_dir = get_label_dir(download_root, label)
        st = get_status(label_dir, row["include_annotation"], row["include_protein"])

        if st["status"] == "not_started":
            log.info("%s: not_started; run download first", label)
            summaries.append(f"{label}\t{bioproject}\tskipped (not_started)")
            continue
        if st["status"] == "complete":
            log.info("%s: already complete; nothing to repair", label)
            summaries.append(f"{label}\t{bioproject}\tskipped (complete)")
            continue

        repaired = False

        # Rehydrate if needed
        if st["needs_rehydration"]:
            rehydrate_cmd = [DATASETS_CMD, "rehydrate", "--directory", str(label_dir)]
            if rehydrate_workers is not None:
                rehydrate_cmd.extend(["--max-workers", str(rehydrate_workers)])
            if not run_cmd(rehydrate_cmd, dry_run, log):
                log.error("Failed to rehydrate %s; continuing with next row", label)
                summaries.append(f"{label}\t{bioproject}\tFAILED (rehydration)")
                continue
            repaired = True

        # Supplement with WormBase if incomplete
        wb_species = row.get("wormbase_species", "")
        wbps_version = config.get("wbps_version", WBPS_DEFAULT_VERSION)
        skip_wb = config.get("skip_wormbase", False)
        if wb_species and not skip_wb and not dry_run:
            wbps_types: list[str] = []
            if row["include_annotation"] and not st["has_gff"]:
                wbps_types.append("gff3")
            if row["include_protein"] and not st["has_protein"]:
                wbps_types.append("protein")
            if not st["has_genome"]:
                wbps_types.insert(0, "genome")
            if wbps_types:
                wbps_dir = find_wormbase_data_dir(label_dir)
                if not _wbps_download_complete(wbps_dir, wb_species, bioproject, wbps_version, wbps_types):
                    log.info("%s: fetching WormBase ParaSite data for: %s", label, ", ".join(wbps_types))
                    download_from_wormbase(wb_species, bioproject, label_dir, wbps_types, dry_run, log, wbps_version)
                    repaired = True

        # Re-validate and write manifest
        if not dry_run:
            valid, missing_gff, missing_protein = validate_download(label_dir, row["include_annotation"], row["include_protein"], log)
            if valid:
                manifest_data = discover_manifest_paths(label_dir, row["include_protein"])
                manifest_data["bioproject"] = bioproject
                manifest_data["label"] = label
                metadata_fields = ["taxon", "taxon_id", "assembly_id", "assembly_level", "source", "genome_size", "reference",
                                   "wormbase_species"]
                for field in metadata_fields:
                    if field in row and row[field]:
                        manifest_data[field] = row[field]
                write_manifest(label_dir, manifest_data, dry_run=False, log=log)
                missing = []
                if missing_gff:
                    missing.append("gff3")
                if missing_protein:
                    missing.append("protein")
                action = "repaired" if repaired else "unchanged"
                if missing:
                    summaries.append(f"{label}\t{bioproject}\t{action} (missing {', '.join(missing)})")
                else:
                    summaries.append(f"{label}\t{bioproject}\t{action} (complete)")
            else:
                summaries.append(f"{label}\t{bioproject}\tFAILED (no genome FASTA)")

    if summaries:
        print()
        print("label\tbioproject\tresult")
        for s in summaries:
            print(s)

    return 0


def setup_logging(log_path: Path | None, verbose: bool) -> logging.Logger:
    log = logging.getLogger("download_genomes")
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
    ap = argparse.ArgumentParser(description="Download genomes from NCBI BioProjects (CSV input)")
    ap.add_argument("--csv", type=Path, default=Path(DEFAULT_CSV), help="Input bioprojects CSV")
    ap.add_argument("--config", type=Path, default=Path("config.yml"), help="Config YAML (optional)")
    ap.add_argument("--dry-run", action="store_true", help="Print actions, no writes or downloads")
    ap.add_argument("--log-file", type=Path, default=Path(LOG_FILE), help="Log file path")
    ap.add_argument("--no-log-file", action="store_true", help="Do not write log file")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--skip-wormbase", action="store_true",
                     help="Skip WormBase ParaSite supplement for missing annotations/proteins")
    sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("status", help="Report status per row (not_started | dehydrated_only | rehydrated | complete)")
    dl = sub.add_parser("download", help="Download and optionally rehydrate")
    dl.add_argument("--dehydrated", action="store_true", help="Download dehydrated (then rehydrate)")
    dl.add_argument("--force", action="store_true",
                    help="Re-download even if label is already complete or assembly changed")
    sub.add_parser("repair", help="Rehydrate where needed (idempotent)")
    args = ap.parse_args()
    config = load_config(args.config)
    if getattr(args, "skip_wormbase", False):
        config["skip_wormbase"] = True
    log = setup_logging(None if args.no_log_file else args.log_file, args.verbose)
    if args.command == "status":
        return cmd_status(args.csv, config, args.dry_run, log)
    if args.command == "download":
        return cmd_download(args.csv, config, args.dry_run, log,
                            use_dehydrated=getattr(args, "dehydrated", False),
                            force=getattr(args, "force", False))
    if args.command == "repair":
        return cmd_repair(args.csv, config, args.dry_run, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
