"""Run fpocket on the trimmed A+D models and extract the descriptors of the pocket whose
alpha-sphere ('STP') atoms lie closest to OG1 of the catalytic Thr (first residue of the chain-A beta chain)."""
import argparse
import glob
import os
import subprocess
from multiprocessing import Pool

import numpy as np
import pandas as pd


def run_fpocket(pdb):
    stem = os.path.splitext(pdb)[0]
    if not os.path.isdir(stem + '_out'):
        subprocess.run(['fpocket', '-f', pdb, '-D', '3.3'], check=True, stdout=subprocess.DEVNULL)
    return pdb


def closest_pocket(out_pdb, beta_start):
    """Pocket whose alpha spheres ('STP' HETATM records) are nearest to OG1 of the catalytic Thr in chain A.
    STP records are parsed from the text: Biopython would merge them (repeated atom names per residue)
    and silently drop most alpha spheres."""
    og1, stp = [], []
    for line in open(out_pdb):
        if line.startswith('ATOM') and line[21] == 'A' and line[12:16].strip() == 'OG1' and int(line[22:26]) == beta_start:
            og1.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
        elif line.startswith('HETATM') and line[17:20] == 'STP':
            stp.append((int(line[22:26]), np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    if len(og1) != 1 or not stp:
        return None, None
    d, pocket = min((np.linalg.norm(xyz - og1[0]), res) for res, xyz in stp)
    return pocket, d


def parse_info(info_file, pocket):
    data, capture = {}, False
    with open(info_file) as f:
        for line in f:
            if f"Pocket {pocket} :" in line:
                capture = True
            elif line.startswith("Pocket") and capture:
                break
            elif capture and ":" in line:
                k, v = line.strip().split(":", 1)
                data[k.strip()] = float(v.strip())
    return data or None


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fpocket_dir", required=True)
    p.add_argument("--cleavage_csv", required=True)
    p.add_argument("--output_csv", required=True)
    p.add_argument("--n_jobs", type=int, default=4)
    a = p.parse_args()
    beta = dict(pd.read_csv(a.cleavage_csv).values)
    pdbs = sorted(glob.glob(f"{a.fpocket_dir}/*_trimmed_ad.pdb"))
    with Pool(a.n_jobs) as pool:
        pool.map(run_fpocket, pdbs)
    rows = []
    for pdb in pdbs:
        name = os.path.basename(pdb).split('_')[0]
        d = f"{a.fpocket_dir}/{name}_trimmed_ad_out"
        pocket, dist = closest_pocket(f"{d}/{name}_trimmed_ad_out.pdb", int(beta[name]))
        if pocket is None:
            print(f"{name}: catalytic Thr OG1 not found / no pocket"); continue
        data = parse_info(f"{d}/{name}_trimmed_ad_info.txt", pocket)
        if data:
            print(f"{name}: pocket {pocket}, dist {dist:.2f} A")
            rows.append({'enzyme': name, **data})
    pd.DataFrame(rows).to_csv(a.output_csv, index=False, na_rep='NaN')
    print(f"Wrote {len(rows)} rows to {a.output_csv}")
