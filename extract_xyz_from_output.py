import os
import re

# Source and destination folders
source_folder = "./por-out"       # Your .out files location
destination_folder = "./por-geom-xyz"  # Your .xyz output location
os.makedirs(destination_folder, exist_ok=True)

# Mapping of atomic numbers or indices to element symbols (optional, just using atom names)
def extract_xyz_from_out(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    # Extract the geometry block
    match = re.search(r"Standard Nuclear Orientation \(Angstroms\)(.*?)Nuclear Repulsion Energy", content, re.DOTALL)
    if not match:
        print(f"Geometry not found in {file_path}")
        return

    geometry_block = match.group(1)

    # Extract lines with atom info
    lines = geometry_block.strip().splitlines()
    atom_lines = [line for line in lines if re.match(r"\s*\d+\s+[A-Z][a-z]?\s", line)]

    # Build .xyz content
    xyz_lines = []
    xyz_lines.append(str(len(atom_lines)))               # Number of atoms
    xyz_lines.append("porphyrin dimer")                  # Comment line

    for line in atom_lines:
        parts = line.split()
        symbol = parts[1]
        x, y, z = parts[2:5]
        xyz_lines.append(f"{symbol:2} {float(x):>12.6f} {float(y):>12.6f} {float(z):>12.6f}")

    # Save to .xyz
    base_name = os.path.basename(file_path).replace(".out", ".xyz")
    xyz_path = os.path.join(destination_folder, base_name)
    with open(xyz_path, 'w') as out_file:
        out_file.write("\n".join(xyz_lines))

    print(f"Saved: {xyz_path}")


# Loop over all .out files
for filename in os.listdir(source_folder):
    if filename.endswith(".out"):
        extract_xyz_from_out(os.path.join(source_folder, filename))
