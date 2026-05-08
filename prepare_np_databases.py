#!/usr/bin/env python3
"""
Prepare NP Atlas and COCONUT databases for virtual screening with AutoDock Vina.

Downloads both databases, cleans SMILES (removes salts/fragments, invalid entries),
and creates:
  1. Clean .smi files for docking (one SMILES per line with compound ID)
  2. Metadata TSV files mapping compound ID to name, source organism, MW, etc.

Run on Expanse: python3 prepare_np_databases.py

Output directory: ~/final_project/np_databases/
"""
import csv
import os
import subprocess
import sys
import zipfile

OUTPUT_DIR = os.path.expanduser("~/final_project/np_databases")

NPATLAS_URL = "https://www.npatlas.org/static/downloads/NPAtlas_download.tsv"
COCONUT_URL = "https://coconut.s3.uni-jena.de/prod/downloads/2026-05/coconut_csv-05-2026.zip"


def download_file(url, dest):
    if os.path.exists(dest):
        print(f"  Already downloaded: {dest}")
        return True
    print(f"  Downloading {url}...")
    result = subprocess.run(["curl", "-sL", url, "-o", dest], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"  ERROR: Download failed: {result.stderr.decode()}")
        return False
    print(f"  Saved: {dest} ({os.path.getsize(dest) / 1024 / 1024:.1f} MB)")
    return True


def clean_smiles(smi):
    """Remove salts/fragments (keep largest component), skip invalid SMILES."""
    if not smi or smi.strip() == "":
        return None
    smi = smi.strip()
    if "." in smi:
        parts = smi.split(".")
        smi = max(parts, key=len)
    if len(smi) < 3:
        return None
    bad_patterns = ["[Na+]", "[K+]", "[Ca+", "[Cl-]", "[Br-]", "[I-]", "[OH-]"]
    if any(bp in smi for bp in bad_patterns):
        remaining = smi
        for bp in bad_patterns:
            remaining = remaining.replace(bp, "")
        if len(remaining.strip()) < 3:
            return None
    return smi


def process_npatlas(raw_path):
    """Process NP Atlas TSV into clean SMILES + metadata files."""
    print("\n--- Processing NP Atlas ---")
    smi_path = os.path.join(OUTPUT_DIR, "npatlas.smi")
    meta_path = os.path.join(OUTPUT_DIR, "npatlas_metadata.tsv")

    total = 0
    kept = 0
    with open(raw_path, "r", encoding="utf-8") as fin, \
         open(smi_path, "w") as fsmi, \
         open(meta_path, "w") as fmeta:

        reader = csv.DictReader(fin, delimiter="\t")
        fmeta.write("compound_id\tname\tsmiles\tmw\torigin_type\tgenus\tspecies\tnpaid\n")

        for row in reader:
            total += 1
            smiles = clean_smiles(row.get("compound_smiles", ""))
            if smiles is None:
                continue

            cid = row.get("npaid", f"NPA{total}")
            name = row.get("compound_name", "").replace("\t", " ")
            mw = row.get("compound_molecular_weight", "")
            origin = row.get("origin_type", "")
            genus = row.get("genus", "")
            species = row.get("origin_species", "")

            fsmi.write(f"{smiles} {cid}\n")
            fmeta.write(f"{cid}\t{name}\t{smiles}\t{mw}\t{origin}\t{genus}\t{species}\t{cid}\n")
            kept += 1

    print(f"  Total compounds: {total}")
    print(f"  Clean compounds: {kept} ({kept*100//total}%)")
    print(f"  SMILES file: {smi_path}")
    print(f"  Metadata file: {meta_path}")
    return kept


def process_coconut(raw_path):
    """Process COCONUT CSV into clean SMILES + metadata files."""
    print("\n--- Processing COCONUT ---")
    smi_path = os.path.join(OUTPUT_DIR, "coconut.smi")
    meta_path = os.path.join(OUTPUT_DIR, "coconut_metadata.tsv")

    total = 0
    kept = 0
    with open(raw_path, "r", encoding="utf-8", errors="replace") as fin, \
         open(smi_path, "w") as fsmi, \
         open(meta_path, "w") as fmeta:

        reader = csv.DictReader(fin)
        fmeta.write("compound_id\tname\tsmiles\tmw\torganisms\tchemical_class\tnp_class\tidentifier\n")

        for row in reader:
            total += 1
            smiles = clean_smiles(row.get("canonical_smiles", ""))
            if smiles is None:
                continue

            cid = row.get("identifier", f"CNP{total}")
            name = row.get("name", "").replace("\t", " ")
            mw = row.get("molecular_weight", "")
            organisms = row.get("organisms", "").replace("\t", " ")
            chem_class = row.get("chemical_class", "").replace("\t", " ")
            np_class = row.get("np_classifier_class", "").replace("\t", " ")

            fsmi.write(f"{smiles} {cid}\n")
            fmeta.write(f"{cid}\t{name}\t{smiles}\t{mw}\t{organisms}\t{chem_class}\t{np_class}\t{cid}\n")
            kept += 1

    print(f"  Total compounds: {total}")
    print(f"  Clean compounds: {kept} ({kept*100//total}%)")
    print(f"  SMILES file: {smi_path}")
    print(f"  Metadata file: {meta_path}")
    return kept


def main():
    print("=" * 60)
    print("NP Database Preparation for Virtual Screening")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    npatlas_raw = os.path.join(OUTPUT_DIR, "NPAtlas_download.tsv")
    coconut_zip = os.path.join(OUTPUT_DIR, "coconut_csv.zip")
    coconut_raw = os.path.join(OUTPUT_DIR, "coconut_csv-05-2026.csv")

    print("\n1. Downloading NP Atlas...")
    if not download_file(NPATLAS_URL, npatlas_raw):
        print("ERROR: Could not download NP Atlas")
        sys.exit(1)

    print("\n2. Downloading COCONUT...")
    if not download_file(COCONUT_URL, coconut_zip):
        print("ERROR: Could not download COCONUT")
        sys.exit(1)

    if not os.path.exists(coconut_raw):
        print("  Extracting COCONUT CSV...")
        with zipfile.ZipFile(coconut_zip, "r") as z:
            z.extractall(OUTPUT_DIR)
        print("  Extracted")

    print("\n3. Processing databases...")
    npa_count = process_npatlas(npatlas_raw)
    coc_count = process_coconut(coconut_raw)

    print("\n" + "=" * 60)
    print("DONE")
    print(f"NP Atlas: {npa_count} compounds -> {OUTPUT_DIR}/npatlas.smi")
    print(f"COCONUT:  {coc_count} compounds -> {OUTPUT_DIR}/coconut.smi")
    print(f"Metadata files in same directory for compound identification")
    print("")
    print("Combined total: " + str(npa_count + coc_count) + " compounds")
    print("=" * 60)


if __name__ == "__main__":
    main()
