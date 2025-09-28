import pandas as pd
import numpy as np
from ase.io import read

def _best_fit_plane_normal(points):
    # Subtract centroid, SVD → last right-singular vector = normal
    centered = points - points.mean(axis=0)
    _, _, vh = np.linalg.svd(centered)
    n = vh[-1]
    n /= np.linalg.norm(n)
    return n

def _principal_inplane_axis(points, normal):
    # PCA of anchor atoms projected into plane → first principal axis
    centered = points - points.mean(axis=0)
    # Project to plane
    proj = centered - np.outer(centered @ normal, normal)
    # SVD of projected coords
    u, s, vh = np.linalg.svd(proj, full_matrices=False)
    e1 = vh[0]
    # Ensure e1 is orthogonal to normal and unit-length
    e1 -= (e1 @ normal) * normal
    e1 /= np.linalg.norm(e1)
    return e1

def compute_slippage_general(xyz_path, mon1_idx=range(0,37), mon2_idx=range(37,74), anchor_symbol="N"):
    try:
        atoms = read(xyz_path)
        pos = atoms.get_positions()
        sym = atoms.get_chemical_symbols()

        # Fe positions
        fe1 = np.array([pos[i] for i in mon1_idx if sym[i] == "Fe"])[0]
        fe2 = np.array([pos[i] for i in mon2_idx if sym[i] == "Fe"])[0]

        # Anchor atoms (e.g., the 4 porphyrin N's) for each monomer
        n1 = np.array([pos[i] for i in mon1_idx if sym[i] == anchor_symbol])
        n2 = np.array([pos[i] for i in mon2_idx if sym[i] == anchor_symbol])

        if len(n1) < 3 or len(n2) < 3:
            raise ValueError("Not enough anchor atoms to define planes.")

        # Best-fit plane normals
        n1_norm = _best_fit_plane_normal(n1)
        n2_norm = _best_fit_plane_normal(n2)

        # Make normals co-directional (avoid near anti-parallel sum)
        if n1_norm @ n2_norm < 0:
            n2_norm = -n2_norm

        # Average normal
        n_avg = n1_norm + n2_norm
        n_avg /= np.linalg.norm(n_avg)

        # Fe–Fe vector and decomposition
        v = fe2 - fe1
        d_perp = abs(v @ n_avg)                          # interlayer
        v_para = v - (v @ n_avg) * n_avg                 # in-plane component
        s_mag = np.linalg.norm(v_para)                   # lateral slippage magnitude

        # Define in-plane axes tied to monomer 1
        e1 = _principal_inplane_axis(n1, n1_norm)        # principal axis in plane
        e2 = np.cross(n_avg, e1)                         # completes right-handed triad
        e2 /= np.linalg.norm(e2)

        # Components and direction
        s1 = v_para @ e1
        s2 = v_para @ e2
        phi = np.degrees(np.arctan2(s2, s1))             # angle in plane (deg)

        # Plane-tilt angle between monomers (deg)
        tilt = np.degrees(np.arccos(np.clip(n1_norm @ n2_norm, -1.0, 1.0)))

        return {
            "Interlayer_Dist_Ang": d_perp,
            "Lateral_Slippage_Ang": s_mag,
            "Slippage_X_Ang": s1,
            "Slippage_Y_Ang": s2,
            "Slippage_Angle_Deg": phi,
            "Plane_Tilt_Deg": tilt,
            "FeFe_Dist_Ang": np.linalg.norm(v)
        }

    except Exception as e:
        print(f"Error processing {xyz_path}: {e}")
        return {
            "Interlayer_Dist_Ang": np.nan,
            "Lateral_Slippage_Ang": np.nan,
            "Slippage_X_Ang": np.nan,
            "Slippage_Y_Ang": np.nan,
            "Slippage_Angle_Deg": np.nan,
            "Plane_Tilt_Deg": np.nan,
            "FeFe_Dist_Ang": np.nan
        }

# ---------- Batch over your 15 principal structures ----------
# CSV columns (as you pasted):
# Principal Structure Index,Coupling / eV,Fe-Fe Distance / Å,Shortest Atom-Atom Distance / Å,Parallel Displacement / Å
df = pd.read_csv("principal_structures.csv")

records = []
for idx in df["Principal Structure Index"]:
    xyz_file = f"structure_{idx}.xyz"  # your filenames
    geom = compute_slippage_general(xyz_file)
    row = {"Structure": idx, "Filename": xyz_file}
    row.update(geom)
    records.append(row)

geom_df = pd.DataFrame(records)

# Merge with your table (left join on index/Structure)
out = df.merge(geom_df, left_on="Principal Structure Index", right_on="Structure", how="left")
out.to_csv("principal_structures_with_general_slippage.csv", index=False)
print("Saved -> principal_structures_with_general_slippage.csv")
