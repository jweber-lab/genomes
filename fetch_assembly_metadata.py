#!/usr/bin/env python3
"""
Fetch genome_size and reference (publication) from NCBI for assemblies in bioprojects.csv.
Uses Entrez API or direct HTTP requests to get assembly metadata.
"""

import argparse
import csv
import json
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET

try:
    from Bio import Entrez
    ENTREZ_AVAILABLE = True
except ImportError:
    ENTREZ_AVAILABLE = False


def fetch_with_entrez_api(assembly_id: str) -> dict[str, Any] | None:
    """Try to fetch assembly metadata using Entrez API."""
    if not ENTREZ_AVAILABLE:
        return None
    try:
        Entrez.email = "genomes@example.com"  # NCBI requires email
        # Search for assembly
        handle = Entrez.esearch(db="assembly", term=f"{assembly_id}[Assembly Accession]", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if not record.get("IdList"):
            return None
        uid = record["IdList"][0]
        # Get summary
        handle = Entrez.esummary(db="assembly", id=uid)
        summary = Entrez.read(handle)
        handle.close()
        if not summary or "DocumentSummarySet" not in summary:
            return None
        assm = summary["DocumentSummarySet"]["DocumentSummary"][0]
        out = {}
        # Get genome size
        if "Biosample" in assm:
            biosample = assm["Biosample"]
            if isinstance(biosample, dict) and "genome_size" in biosample:
                out["genome_size"] = str(biosample["genome_size"])
        # Get publication
        if "Publication" in assm and assm["Publication"]:
            pubs = assm["Publication"]
            if isinstance(pubs, list) and pubs:
                pub = pubs[0]
                if isinstance(pub, dict):
                    if "PubmedId" in pub:
                        out["reference"] = f"PMID:{pub['PubmedId']}"
                    elif "Doi" in pub:
                        out["reference"] = f"DOI:{pub['Doi']}"
        return out if out else None
    except Exception as e:
        print(f"Entrez API error for {assembly_id}: {e}", file=sys.stderr)
        return None


def fetch_with_http(assembly_id: str) -> dict[str, Any] | None:
    """Try to fetch assembly metadata using direct HTTP request to NCBI."""
    try:
        # First search for the assembly to get its UID
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={urllib.parse.quote(assembly_id)}[Assembly+Accession]&retmode=xml&retmax=1"
        with urllib.request.urlopen(search_url, timeout=10) as response:
            search_xml = response.read()
        search_root = ET.fromstring(search_xml)
        id_list = search_root.find(".//IdList")
        if id_list is None or len(id_list) == 0:
            return None
        uid = id_list[0].text
        if not uid:
            return None
        time.sleep(0.35)  # Rate limiting
        # Now get summary using UID
        summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={uid}&retmode=xml"
        with urllib.request.urlopen(summary_url, timeout=10) as response:
            xml_data = response.read()
        root = ET.fromstring(xml_data)
        out = {}
        # Parse XML for genome size and publication
        docsum = root.find(".//DocumentSummary")
        if docsum is not None:
            # Genome size - check multiple possible locations
            biosample = docsum.find("Biosample")
            if biosample is not None:
                genome_size_elem = biosample.find("genome_size")
                if genome_size_elem is not None and genome_size_elem.text:
                    out["genome_size"] = genome_size_elem.text.strip()
            # Also check Stats
            stats = docsum.find("Stats")
            if stats is not None and "genome_size" not in out:
                genome_size_elem = stats.find("genome_size")
                if genome_size_elem is not None and genome_size_elem.text:
                    out["genome_size"] = genome_size_elem.text.strip()
            # Publication
            pub_elem = docsum.find("Publication")
            if pub_elem is not None:
                pubmed_id = pub_elem.find("PubmedId")
                if pubmed_id is not None and pubmed_id.text:
                    out["reference"] = f"PMID:{pubmed_id.text.strip()}"
                else:
                    doi_elem = pub_elem.find("Doi")
                    if doi_elem is not None and doi_elem.text:
                        out["reference"] = f"DOI:{doi_elem.text.strip()}"
        return out if out else None
    except Exception as e:
        print(f"HTTP fetch error for {assembly_id}: {e}", file=sys.stderr)
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch genome_size and reference from NCBI for assemblies")
    ap.add_argument("--csv", type=Path, default=Path("bioprojects.csv"), help="Input CSV")
    ap.add_argument("--dry-run", action="store_true", help="Print what would be fetched")
    args = ap.parse_args()
    rows = []
    updated = 0
    with open(args.csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        for row in reader:
            assembly_id = (row.get("assembly_id") or "").strip()
            if not assembly_id:
                rows.append(row)
                continue
            # Skip if already has both
            if row.get("genome_size") and row.get("reference"):
                rows.append(row)
                continue
            if args.dry_run:
                print(f"Would fetch: {assembly_id}")
                rows.append(row)
                continue
            # Try Entrez API first, then HTTP fallback
            print(f"Fetching metadata for {assembly_id}...", file=sys.stderr)
            metadata = fetch_with_entrez_api(assembly_id) or fetch_with_http(assembly_id)
            if metadata:
                if "genome_size" in metadata and not row.get("genome_size"):
                    row["genome_size"] = metadata["genome_size"]
                    print(f"  Found genome_size: {metadata['genome_size']}", file=sys.stderr)
                if "reference" in metadata and not row.get("reference"):
                    row["reference"] = metadata["reference"]
                    print(f"  Found reference: {metadata['reference']}", file=sys.stderr)
                updated += 1
            else:
                print(f"  No metadata found", file=sys.stderr)
            rows.append(row)
            # Rate limiting: wait 0.5 seconds between requests
            if not args.dry_run:
                time.sleep(0.5)
    if not args.dry_run and updated > 0:
        with open(args.csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Updated {updated} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
