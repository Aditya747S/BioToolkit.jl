import sys
import os
import time
import subprocess
import tempfile
import xml.etree.ElementTree as ET

print("Generating sequences for BLAST+ comparison...")
# Target: 100,000 bp
target_base = "ACGT" * 25000
target = target_base[:50000] + "GGCCGATCGATCGATCGATCGA" + target_base[50022:]

# Query: 1,000 bp
query = ("T" * 500) + "GGCCGATCGATCGATCGATCGA" + ("A" * 478)

with tempfile.TemporaryDirectory() as tmpdir:
    query_file = os.path.join(tmpdir, "query.fa")
    target_file = os.path.join(tmpdir, "target.fa")
    db_name = os.path.join(tmpdir, "target_db")
    out_file = os.path.join(tmpdir, "blast_out.xml")

    with open(query_file, "w") as f:
        f.write(">query\n" + query + "\n")
    
    with open(target_file, "w") as f:
        f.write(">target\n" + target + "\n")

    print("Building NCBI BLAST database...")
    subprocess.run(["makeblastdb", "-in", target_file, "-dbtype", "nucl", "-out", db_name], 
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    print("Running NCBI BLASTN binary...")
    
    # Warmup
    subprocess.run(["blastn", "-query", query_file, "-db", db_name, "-outfmt", "5", "-out", out_file], check=True)

    reps = 10
    start_time = time.perf_counter()
    for _ in range(reps):
        subprocess.run(["blastn", "-query", query_file, "-db", db_name, "-outfmt", "5", "-out", out_file], check=True)
    end_time = time.perf_counter()
    
    ms_per_rep = ((end_time - start_time) * 1000) / reps

    print(f"\nNCBI BLAST+ benchmark")
    print(f"  reps={reps}")
    print(f"  ncbi_blast_ms={ms_per_rep:.4f}")
    
    # Verify it found it via simple XML parsing
    tree = ET.parse(out_file)
    root = tree.getroot()
    hsp = root.find(".//Hsp_hit-from")
    if hsp is not None:
        print(f"  Found hit! Target Start: {hsp.text}")
    else:
        print("  WARNING: No hit found.")
