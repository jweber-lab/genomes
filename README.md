# Genomes: download and BLAST databases

This repository is the single source of truth for **which genomes we have** and **how to BLAST against them** (with annotations). It provides:

1. **Download script** — Check status, download, and repair genome data from NCBI BioProjects (input CSV). Automatically falls back to ENA download for ENA/DDBJ projects when NCBI download fails. Supplements missing annotations and proteins from [WormBase ParaSite](https://parasite.wormbase.org/) when configured.
2. **BLAST DB script** — Build and manage local BLAST databases (nucleotide and optional protein) and a **registry** for downstream pipelines.

Downstream scripts (e.g. poolseq annotation) read the registry and run BLAST against these DBs; they adapt to this repo’s interface.

---

## Install

### Option 1: Conda (recommended)

Create a conda environment with all dependencies:

```bash
conda env create -f environment.yml
conda activate genomes
```

This installs:
- Python 3.8+
- PyYAML (for config/registry YAML files)
- NCBI Datasets CLI (for downloads)
- BLAST+ (for building and querying databases)
- enaBrowserTools (optional, for ENA fallback downloads)

### Option 2: Manual installation

**Python dependencies:**
```bash
pip install -r requirements.txt
```

**CLI tools** (install separately):
- **NCBI Datasets CLI** (for downloads):
  - [Download and install](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
  - Or: `conda install -c conda-forge ncbi-datasets-cli`
- **enaBrowserTools** (optional, for ENA fallback downloads):
  - `conda install -c bioconda enabrowsertools`
  - Used automatically when NCBI download fails for ENA/DDBJ projects
- **BLAST+** (for building and querying DBs):
  - [BLAST+ executables](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
  - Or: `conda install -c bioconda blast`

---

## Input CSV: `bioprojects.csv`

One row per NCBI BioProject. Columns:

| Column               | Required | Description |
|----------------------|----------|-------------|
| `bioproject`         | Yes      | NCBI BioProject accession (e.g. `PRJNA123`, `PRJEB33226`) |
| `label`              | No       | Short name for dirs and registry (default: bioproject) |
| `taxon_id`           | No       | NCBI taxonomy ID (e.g. `94845`) — stored in manifest and registry |
| `skip`               | No       | `true`/`false`; skip this row (default: false) — rows with `skip=true` are ignored |
| `taxon`              | No       | Scientific name (e.g. `Ligula intestinalis`) — stored in manifest and registry |
| `wormbase_species`   | No       | WormBase ParaSite species key (e.g. `hymenolepis_microstoma`). When set, the script automatically downloads curated annotations and proteins from WormBase ParaSite FTP for this bioproject |
| `doi`                | No       | Publication DOI (e.g. `10.1234/example`) |
| `assembly_id`        | No       | Specific assembly accession (e.g. `GCA_036362985.1`) — if provided, used for download instead of bioproject (more specific) |
| `assembly_level`     | No       | Filter: `Complete`, `Chromosome`, `Scaffold`, `Contig` |
| `genome_size`        | No       | Genome size in bp (auto-filled by `fetch_assembly_metadata.py` if available) |

**Note:** Column order doesn't matter — scripts read by column name. Example:

```csv
bioproject,label,taxon_id,skip,taxon,wormbase_species,doi,assembly_id,assembly_level,genome_size
PRJEB124,hmic,85433,false,Hymenolepis microstoma,hymenolepis_microstoma,10.1038/nature12031,GCA_000469805.3,scaffold,
```

**Note:** 
- If `assembly_id` is provided, the download script uses it instead of `bioproject` for a more specific download.
- Metadata fields (`taxon`, `taxon_id`, `assembly_id`, `source`, `genome_size`, `reference`, `wormbase_species`) are stored in `manifest.json` and included in `registry.yaml` entries for downstream use. The `source` field is inferred from the bioproject prefix (`PRJNA` → NCBI, `PRJEB`/`PRJEA` → ENA, `PRJDB` → DDBJ).
- Rows with `skip=true` are ignored by both scripts.
- Rows are processed in CSV order.
- **Multiple assemblies per species:** Add each bioproject as a separate row. Set `skip=true` for alternative assemblies you don't actively use. This keeps a record of available options.
- Use `fetch_assembly_metadata.py` to attempt auto-filling `genome_size` and `reference` from NCBI (may require manual updates if API doesn't return data).
- **ENA/DDBJ projects**: BioProjects starting with `PRJEB` (ENA) or `PRJDB` (DDBJ) are automatically detected. If NCBI download fails for these projects and `assembly_id` is provided, the script automatically attempts an ENA download using `enaBrowserTools`. ENA downloads may not include GFF/protein files; check ENA FTP if needed.
- **WormBase ParaSite**: When `wormbase_species` is set, WormBase ParaSite is the **primary data source**. The script downloads the genome assembly, curated annotations (GFF3), and proteins from the [WormBase ParaSite FTP](https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/) using the row's own `bioproject`. NCBI/ENA is only used as a fallback if WormBase doesn't provide a complete set. Only set `wormbase_species` for bioprojects that exist on WormBase ParaSite — if the bioproject doesn't match, the annotations won't correspond to the assembly. Use `--skip-wormbase` to disable. Set `wbps_version` in `config.yml` to target a specific WormBase release (default: `WBPS19`).

---

## Quick start

**1. Download genomes**

```bash
# Check status
python download_genomes.py --csv bioprojects.csv status

# Download (and rehydrate if dehydrated)
python download_genomes.py --csv bioprojects.csv download

# Repair incomplete/dehydrated
python download_genomes.py --csv bioprojects.csv repair
```

Optional: use `config.yml` for `download_root`, `rehydrate_workers`, `wbps_version`, etc. Use `--dry-run` to print actions without downloading. Use `--skip-wormbase` to disable the automatic WormBase ParaSite supplement.

**2. Build BLAST databases**

```bash
# Build nucl (and optional prot) DB per label
python setup_blast_db.py --csv bioprojects.csv build

# Optional: create alias DB(s) combining all labels (e.g. all_nucl)
python setup_blast_db.py --csv bioprojects.csv build --alias

# Or create alias DBs later (without rebuilding per-label DBs)
python setup_blast_db.py alias

# List DBs in registry
python setup_blast_db.py list
```

Use `--config config.yml` for `blast_db_root`, `build_jobs`, etc. Use `--dry-run` to print actions without building.

---

## Directory layout

- `downloads/<label>/` — NCBI dataset (unzipped); `manifest.json` written after download/repair.
- `downloads/<label>/wormbase_parasite/` — Supplementary files from WormBase ParaSite (annotations, proteins).
- `blast_db/<label>/nucl/` — Nucleotide DB (blastn, tblastn) and `annotations.gff`.
- `blast_db/<label>/prot/` — Protein DB (blastp, blastx).
- `blast_db/registry.yaml` — **Consumer contract**: list of DBs and how to use them.

---

## Registry format

Downstream pipelines should read `blast_db/registry.yaml`. Schema:

```yaml
databases:
  - name: <label>_nucl
    db_path: /path/to/blast_db/<label>/nucl/sequences
    db_type: nucl
    gff_path: /path/to/blast_db/<label>/nucl/annotations.gff
    blast_programs: [blastn, tblastn]
    bioproject: PRJNA123
    taxon: Species name
    taxon_id: 12345
    assembly_id: GCA_000123456.1
    assembly_level: scaffold
  - name: <label>_prot
    db_path: /path/to/blast_db/<label>/prot/sequences
    db_type: prot
    gff_path: /path/optional
    blast_programs: [blastp, blastx]
    bioproject: PRJNA123
    taxon: Species name
    taxon_id: 12345
    assembly_id: GCA_000123456.1
```

- **db_path** — BLAST DB prefix (no extension); use e.g. `blastn -db <db_path>`.
- **db_type** — `nucl` or `prot`.
- **gff_path** — GFF for mapping hit IDs to genes/features (mainly for nucl).
- **blast_programs** — Which programs to use: `blastn`, `tblastn`, `blastp`, `blastx`.

Filter by `blast_programs` (e.g. use the row where `blastp` is in `blast_programs` for protein-vs-protein search).

---

## Using the databases

### Command-line usage

**List available databases:**
```bash
python setup_blast_db.py list
```

**Query a specific database:**
```bash
# Nucleotide search (blastn)
blastn -db blast_db/lint/nucl/sequences -query my_sequences.fasta -out results.txt

# Protein search (blastp)
blastp -db blast_db/spro/prot/sequences -query my_proteins.faa -out results.txt

# Translated search (blastx: translate query, search nucl DB)
blastx -db blast_db/lint/nucl/sequences -query my_proteins.faa -out results.txt

# Translated search (tblastn: translate DB, search with nucl query)
tblastn -db blast_db/lint/nucl/sequences -query my_proteins.faa -out results.txt
```

**Query all genomes at once (using alias DB):**
```bash
# Search all nucleotide databases
blastn -db blast_db/all_nucl -query my_sequences.fasta -out results.txt

# Search all protein databases
blastp -db blast_db/all_prot -query my_proteins.faa -out results.txt
```

**Get database info:**
```bash
blastdbcmd -db blast_db/lint/nucl/sequences -info
```

### Using in Python scripts

**Read the registry:**
```python
import yaml
from pathlib import Path

# Load registry
registry_path = Path("blast_db/registry.yaml")
with open(registry_path) as f:
    registry = yaml.safe_load(f)

# Find a specific database
for db in registry["databases"]:
    if db["name"] == "lint_nucl":
        db_path = db["db_path"]
        gff_path = db.get("gff_path")
        print(f"Database: {db_path}")
        print(f"GFF: {gff_path}")
        break

# Find all nucleotide databases
nucl_dbs = [db for db in registry["databases"] 
            if db["db_type"] == "nucl" and "blastn" in db["blast_programs"]]

# Find all databases for a specific taxon
taenia_dbs = [db for db in registry["databases"] 
              if db.get("taxon", "").startswith("Taenia")]

# Use alias database (all genomes)
all_nucl = next((db for db in registry["databases"] if db["name"] == "all_nucl"), None)
if all_nucl:
    db_path = all_nucl["db_path"]
    # Use db_path with blastn, blastx, tblastn
```

**Run BLAST from Python:**
```python
import subprocess
from pathlib import Path

# Get database path from registry
db_path = "blast_db/lint/nucl/sequences"
query_file = "my_sequences.fasta"
output_file = "results.txt"

# Run blastn
cmd = [
    "blastn",
    "-db", db_path,
    "-query", query_file,
    "-out", output_file,
    "-outfmt", "6",  # Tabular format
    "-max_target_seqs", "10"
]
subprocess.run(cmd, check=True)
```

**Map BLAST hits to annotations (using GFF):**
```python
from pathlib import Path
import gffutils  # or use BioPython, pybedtools, etc.

# Get GFF path from registry
gff_path = Path("blast_db/lint/nucl/annotations.gff")

# Parse GFF to map sequence IDs to genes/features
# Example using gffutils:
import gffutils
db = gffutils.create_db(str(gff_path), dbfn=":memory:", force=True)

# After BLAST, map hit IDs to features
hit_id = "scaffold_123_456"
features = list(db.features_of_type("gene", limit=(hit_id, hit_id)))
for feat in features:
    print(f"Gene: {feat.id}, Location: {feat.start}-{feat.end}")
```

### Program selection guide

| Query type | Database type | Program | Example use case |
|------------|---------------|---------|------------------|
| Nucleotide | nucl | `blastn` | Find similar DNA sequences |
| Protein | prot | `blastp` | Find similar protein sequences |
| Protein | nucl | `blastx` | Translate query protein, search nucleotide DB |
| Nucleotide | nucl | `tblastn` | Translate nucleotide DB, search with protein query |

**Note:** For alias databases (`all_nucl`, `all_prot`), use the same programs as for individual databases. The alias simply combines multiple databases into one query.

---

## Alias DBs

With `--alias` (or the `alias` subcommand), the script creates combined alias DBs (e.g. `blast_db/all_nucl`, `blast_db/all_prot`) via `blastdb_aliastool`. They appear as extra entries in the registry. Use them to query all genomes in one go without merging FASTA.

**Create alias databases after building:**
```bash
# If you built without --alias, add aliases later
python setup_blast_db.py alias
```

---

## Logging and options

- **Download script**: logs to `download.log` (override with `--log-file`; disable with `--no-log-file`).
- **BLAST DB script**: logs to `blast_db.log`.
- **Dry-run**: both scripts support `--dry-run` (no writes or downloads).
- **Validation**: download script checks files exist and size > 0 after download; BLAST script runs `blastdbcmd -db <path> -info` after each build.

---

## Memory and threads

- **makeblastdb** is single-threaded. Parallelization is **across labels** (set `build_jobs` in `config.yml`).
- Very large genomes may need substantial memory (BLAST+ 2.10+ uses LMDB for version 5 DBs). See [BLAST+ documentation](https://www.ncbi.nlm.nih.gov/books/NBK279671/) and system resources.

---

## Additional useful CSV fields

The scripts preserve any extra columns in the CSV. Consider adding these fields for better organization:

| Field              | Description |
|--------------------|-------------|
| `version`         | Assembly version or release date (e.g. `2023-01-15`) |
| `notes`            | Free-text notes or description |
| `status`           | Status (e.g. `reference`, `representative`, `latest`) |
| `group`            | Category/group name for organizing genomes |
| `download_url`     | Alternative download URL (if not using NCBI datasets) |
| `include_annotation` | `true`/`false`; request GFF (default: true) — only needed to disable annotations for a specific row |
| `include_protein`    | `true`/`false`; request protein FASTA (default: true) — only needed to disable proteins for a specific row |

**Note:** Currently, these fields are preserved in the CSV row dict but not automatically stored in `manifest.json` or `registry.yaml`. To include them, modify the scripts to copy them (similar to how `taxon`, `taxon_id`, `assembly_id` are handled).
