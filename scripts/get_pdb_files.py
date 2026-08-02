from rcsbapi.search import search_attributes as attrs
from pathlib import Path
import csv

output_file = Path(__file__).resolve().parent.parent / "data" / "benchmark_raw" / "benchmark_manifest.csv"
output_file.parent.mkdir(parents=True, exist_ok=True)

def query_atoms(min_atoms, max_atoms, num_entries=20):
    q1 = attrs.rcsb_entry_info.deposited_atom_count >= min_atoms
    q2 = attrs.rcsb_entry_info.deposited_atom_count <= max_atoms
    query = (q1 & q2)

    results = query(return_type="entry", rows=num_entries, results_verbosity="compact")

    return results

def main():

    sizes = [(1000, 5000), (5000, 10000), (10000, 50000), (50000, 100000), (100000, 1000000), (1000000, 3000000)]

    with open(output_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow([
            "pdb_id", "min_atoms", "max_atoms"
        ])
        for size in sizes: 

            print(f"===== {size[0]}-{size[1]} atoms =====")
            results = query_atoms(size[0], size[1])

            count = 0
                
            for pdb_id in results: 
                print(pdb_id)
                count += 1
                print (count)
                writer.writerow([pdb_id, size[0], size[1]])

                if count >= 20:
                    break
    
if __name__ == "__main__":
    main()

