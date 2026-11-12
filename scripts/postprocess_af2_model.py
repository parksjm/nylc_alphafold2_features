"""Reorder chains and renumber residues of an AlphaFold2 (ColabFold) tetramer-of-heterodimers
model so it matches the NylC X-ray structure (PDB 5y0m).

Raw ColabFold models have chains A-D = alpha chains and E-H = beta chains. After processing,
each alpha/beta pair shares one chain ID (A-D, ordered as in 5y0m) with continuous residue numbering.
Run with PyMOL's python (needs `pymol`)."""
import argparse
import os
import numpy as np
from pymol import cmd
from scipy.optimize import linear_sum_assignment

parser = argparse.ArgumentParser()
parser.add_argument("--model", required=True, help="Relaxed rank-1 ColabFold model (PDB)")
parser.add_argument("--beta_start", required=True, type=int, help="Residue number of beta chain N-terminus")
parser.add_argument("--output_dir", required=True, help="Directory for processed PDB files")
parser.add_argument("--output_dir_fpocket", required=True, help="Directory for trimmed PDB files for fpocket")
parser.add_argument("--ref", required=True, help="Reference X-ray tetramer (5y0m_tetramer.pdb)")
args = parser.parse_args()

ref, ref_name = args.ref, "nylc_xray"
model = args.model
model_name = os.path.splitext(os.path.basename(model))[0].split('_')[0]
model_beta_resi_start = args.beta_start
tmp_chains = ['W', 'X', 'Y', 'Z']
print(f"Processing {model_name}...")


def load_pdb_files():
    cmd.load(ref, ref_name)
    cmd.load(model, model_name)
    # Fixed orientation: the output coordinates are rotated into this frame in write_pdb().
    # Keep it unchanged: Shrake-Rupley SASA and fpocket depend slightly on the coordinate frame,
    # so a different frame changes BSA/fpocket descriptors by ~0.5-1% (and occasionally the chosen pocket).
    cmd.set_view((
        -0.033802517,    0.158777967,    0.986735702,
        -0.999350965,   -0.017669253,   -0.031391911,
         0.012449547,   -0.987155139,    0.159273610,
        -0.000052944,    0.000145555, -306.412078857,
       -21.511169434,  -33.441307068,   42.405929565,
       253.244506836,  359.581939697,  -20.000000000))


def get_model_chain_ids():
    model_chains = cmd.get_chains(model_name)
    assert len(model_chains) == 8, f"{model_name}: expected 8 chains, found {model_chains}"
    return model_chains[:4], model_chains[-4:]


def align_model_to_ref():
    cmd.align(f"{model_name} and chain {alpha_chains[0]}", f"{ref_name} and chain A")


def assign_alpha_chains_temp_names():
    for i, chain in enumerate(alpha_chains):
        cmd.alter(f"({model_name} and chain {chain})", f'chain="{tmp_chains[i]}"')
    cmd.sort(model_name)


def best_chain_order(model_chains, resi_sel):
    """Assign each reference chain (A-D) a distinct model chain, minimizing the total RMSD.
    (Legacy: independent per-chain argmin, which could pick the same model chain twice and
    produce scrambled chain IDs, e.g. for dp171.) Identical to argmin whenever the argmins are distinct."""
    ref_chains = cmd.get_chains(ref_name)
    rms = np.zeros((len(ref_chains), len(model_chains)))
    for r, i in enumerate(ref_chains):
        isel = f"{ref_name} and chain {i} and resi {resi_sel}"
        isel2 = isel + ' & aln'
        for c, j in enumerate(model_chains):
            jsel = f"{model_name} and chain {j}"
            jsel2 = jsel + ' & aln'
            cmd.align(isel, jsel, cycles=0, transform=0, object="aln")
            rms[r, c] = cmd.rms_cur(jsel2, isel2, matchmaker=-1)
            cmd.delete("aln")
    rows, cols = linear_sum_assignment(rms)
    if not np.array_equal(cols, rms.argmin(axis=1)):
        print(f"WARNING {model_name}: per-chain argmin was not a permutation; used optimal unique assignment")
    return [model_chains[c] for c in cols]


def rename_chains(new_chains):
    for i, chain in enumerate(new_chains):
        cmd.alter(f"({model_name} and chain {chain})", f'chain="{alpha_chains[i]}"')
    cmd.sort(model_name)


def renumber_beta_chains():
    offset = model_beta_resi_start - 1
    cmd.alter(f"({model_name} and chain E+F+G+H)", f'resi=str(int(resi)+{offset})')
    cmd.sort()


def write_pdb(output_dir):
    pdbname = os.path.join(output_dir, model_name + '.pdb')
    m = cmd.get_view(1)
    ttt = [m[0], m[1], m[2], 0.0, m[3], m[4], m[5], 0.0, m[6], m[7], m[8], 0.0, 0.0, 0.0, 0.0, 1.0]
    cmd.transform_object(model_name, ttt, transpose=1)
    cmd.save(pdbname, model_name)


def write_fpocket_pdb(output_dir_fpocket):
    """Keep chains A and D; remove the ten alpha-chain residues preceding the cleavage site."""
    fp_pdbname = os.path.join(output_dir_fpocket, model_name + '_trimmed_ad.pdb')
    cmd.remove(f"{model_name} and chain B+C")
    last_alpha_resi = model_beta_resi_start - 10
    for chain in cmd.get_chains(model_name):
        cmd.remove(f"{model_name} and chain {chain} and resi {last_alpha_resi}-{model_beta_resi_start - 1}")
    cmd.save(fp_pdbname, model_name)


load_pdb_files()
alpha_chains, beta_chains = get_model_chain_ids()          # ['A'..'D'], ['E'..'H']
align_model_to_ref()
assign_alpha_chains_temp_names()                           # alpha A-D -> W-Z
rename_chains(best_chain_order(list(tmp_chains), '5-260'))  # W-Z -> A-D in 5y0m order
renumber_beta_chains()
rename_chains(best_chain_order(['E', 'F', 'G', 'H'], '267-340'))  # beta E-H -> A-D in 5y0m order
write_pdb(args.output_dir)
write_fpocket_pdb(args.output_dir_fpocket)
