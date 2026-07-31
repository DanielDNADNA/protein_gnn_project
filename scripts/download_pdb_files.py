import httpx
import csv
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent

output_dir = project_root / "data" / "benchmark_raw"
#creates directory if not there, but ok if it already is
output_dir.mkdir(parents=True, exist_ok=True)

manifest_path = project_root / "data" / "benchmark_raw" / "benchmark_manifest.csv"

with open(manifest_path, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader: 
        pdb_id = row["pdb_id"]
        min_atoms = int(row["min_atoms"])
        max_atoms = int(row["max_atoms"])
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        file_path = output_dir / f"{pdb_id}.cif"

        if not file_path.exists():

            response = httpx.get(url)
            #Python throws exception if error code is returned
            response.raise_for_status()
            file_path.write_bytes(response.content)


