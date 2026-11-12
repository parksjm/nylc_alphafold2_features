"""Buried surface area (BSA) of the A-B and A-D interfaces (Shrake-Rupley, Biopython defaults:
1.4 A probe, 100 points). BSA = (SASA_A + SASA_other - SASA_dimer) / 2."""
import argparse
import glob
import os
from multiprocessing import Pool

import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley


def sasa(entity):
    ShrakeRupley().compute(entity, level="A")
    return sum(a.sasa for a in entity.get_atoms() if a.sasa is not None)


def keep_chains(structure, keep):
    s = structure.copy()
    for chain in list(s[0]):
        if chain.id not in keep:
            s[0].detach_child(chain.id)
    return s


def calculate_bsa(pdb_file):
    structure = PDBParser(QUIET=True).get_structure("protein", pdb_file)
    ids = [c.id for c in structure[0]]
    if 'A' not in ids:
        raise ValueError("Chain A is not present")
    res = {"AB_BSA": None, "AD_BSA": None}
    for other, key in (('B', 'AB_BSA'), ('D', 'AD_BSA')):
        if other not in ids:
            continue
        a = sasa(keep_chains(structure, ['A'])[0]['A'])
        o = sasa(keep_chains(structure, [other])[0][other])
        d = sasa(keep_chains(structure, ['A', other])[0])
        res[key] = round((a + o - d) / 2, 2)
    return res


def work(pdb_file):
    name = os.path.basename(pdb_file).split(".")[0]
    try:
        r = calculate_bsa(pdb_file)
        return {"enzyme": name, **r}
    except Exception as e:  # keep going, but report
        print(f"Error calculating BSA for {name}: {e}")
        return {"enzyme": name, "AB_BSA": None, "AD_BSA": None}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pdb_dir", required=True)
    p.add_argument("--output_csv", required=True)
    p.add_argument("--n_jobs", type=int, default=4)
    a = p.parse_args()
    files = sorted(f for f in glob.glob(f"{a.pdb_dir}/*.pdb") if not os.path.basename(f).lower().startswith('5y0m'))
    with Pool(a.n_jobs) as pool:
        rows = pool.map(work, files)
    pd.DataFrame(rows, columns=["enzyme", "AB_BSA", "AD_BSA"]).to_csv(a.output_csv, index=False, na_rep="NaN")
    print(f"BSA results for {len(rows)} models written to {a.output_csv}")
