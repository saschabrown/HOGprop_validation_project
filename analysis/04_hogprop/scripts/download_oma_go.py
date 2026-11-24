"""Fetch May 2026 OMA or use --database PATH; export GO assignments and embedded OBO."""
import argparse
import bisect
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

import h5py

URL = "https://omabrowser.org/All.May2026/OmaServer.h5"


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def value(item):
    return item.decode() if isinstance(item, bytes) else str(item)


def export(db, folder):
    genomes = sorted((int(g["EntryOff"]), int(g["TotEntries"]), value(g["UniProtSpeciesCode"]))
                     for g in db["Genome"])
    offsets = [g[0] for g in genomes]
    ontology = db["Ontologies/GO"][:].tobytes()
    terms = {line[4:].decode() for line in ontology.splitlines() if line.startswith(b"id: GO:")}
    if not terms or b"[Term]" not in ontology:
        raise ValueError("Embedded ontology is not recognizable GO OBO")
    (folder / "go.obo").write_bytes(ontology)
    table = db["Annotations/GeneOntology"]
    fields = table.dtype.names
    count = 0
    missing = set()
    with gzip.open(folder / "oma_go.tsv.gz", "wt", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["OMA_ID", "GO_ID", *fields])
        for start in range(0, len(table), 250000):
            for row in table[start:start + 250000]:
                entry = int(row["EntryNr"])
                idx = bisect.bisect_left(offsets, entry) - 1
                offset, size, species = genomes[idx]
                if idx < 0 or not 1 <= entry - offset <= size:
                    raise ValueError(f"Unmapped OMA entry {entry}")
                term = f"GO:{int(row['TermNr']):07d}"
                if term not in terms:
                    missing.add(term)
                writer.writerow([f"{species}{entry-offset:05d}", term, *(value(row[k]) for k in fields)])
                count += 1
    if not count or count != len(table):
        raise ValueError("Incomplete or empty GO assignment export")
    return dict(annotation_rows=count, ontology_terms=len(terms), absent_primary_ids=sorted(missing))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, help="existing May 2026 OmaServer.h5 (read only)")
    args = parser.parse_args()
    data = Path(os.environ["WP1_ROOT"]) / "HOGprop_validation_project/analysis/04_hogprop/data"
    data.mkdir(parents=True, exist_ok=True)
    output = data / "oma_go_May2026"
    if output.exists():
        raise SystemExit(f"Output already exists: {output}")
    with tempfile.TemporaryDirectory(dir=data) as tmp:
        work = Path(tmp)
        db_path = args.database
        if db_path is None:
            db_path = work / "OmaServer.h5"
            print(f"Downloading full OMA database from {URL}", flush=True)
            with urlopen(URL, timeout=180) as response, db_path.open("xb") as out:
                shutil.copyfileobj(response, out, 8 << 20)
        bundle = work / "export"
        bundle.mkdir()
        with h5py.File(db_path, "r") as db:
            release = value(db.attrs.get("oma_version", ""))
            if "may2026" not in "".join(release.lower().split()):
                raise ValueError(f"Could not verify May 2026 from oma_version={release!r}")
            record = export(db, bundle)
        record.update(release=release, source=str(args.database.resolve()) if args.database else URL,
                      source_sha256=digest(db_path), h5py=h5py.__version__,
                      generated_utc=datetime.now(timezone.utc).isoformat(),
                      output_sha256={p.name: digest(p) for p in bundle.iterdir()})
        (bundle / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
        bundle.rename(output)
    print(output)


if __name__ == "__main__":
    main()
