import os
import shutil
import sys

missing_file = sys.argv[1] if len(sys.argv) > 1 else  "/home/anirudh/genomes/scripts/missing_faas.txt"

missing_list = []
# --- CONFIG ---
with open(missing_file, 'r') as file:
    lines = file.readlines()
    for line in lines:
        missing = line.strip("\n")
        missing_list.append(missing)


print(len(missing_list))

source_dir = "/home/anirudh/genomes/selected_genomes/prokka_results"
destination_dir = "/home/anirudh/genomes/missing_faa"

# Create destination dir if it doesn't exist
os.makedirs(destination_dir, exist_ok=True)

# Copy matching files
copied = 0
for root, dirs, files in os.walk(source_dir):
    for file in files:
        # Check if it's a .faa file and matches any missing name
        for name in missing_list:
            if file.startswith(name) and file.endswith(".faa"):
                src_path = os.path.join(root, file)
                dst_path = os.path.join(destination_dir, file)
                shutil.copy2(src_path, dst_path)
                print(f"✅ Copied {file}")
                copied += 1

print(f"\nTotal missing files: {len(missing_list)}")
print(f"\nTotal files copied: {copied}")

