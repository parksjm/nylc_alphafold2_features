"""One-time derivation of the repository's small input files from the legacy project files.
Drops the 24 extra picks (nyl096-nyl119, 'dp25' set); nylc-gyaq (the reference) is kept.
Alignment columns are NOT removed so MSA-based features stay comparable with legacy output."""
import sys, re
import pandas as pd
legacy, out = sys.argv[1], sys.argv[2]
extra = {f"nyl{i:03d}" for i in range(96, 120)}
# MSA
keep, name = [], None
for line in open(f"{legacy}/nylc_all.mafft.fasta"):
    if line.startswith('>'):
        name = line[1:].split('/')[0].lower(); on = name not in extra
    if on: keep.append(line)
open(f"{out}/nylc_all.mafft.fasta", 'w').writelines(keep)
# sequences
x = pd.read_excel(f"{legacy}/nylc_all.xlsx")[['Diversity set ID', 'AAseq']]
x.columns = ['enzyme', 'AAseq']; x['enzyme'] = x.enzyme.str.lower()
x[~x.enzyme.isin(extra)].to_csv(f"{out}/nylc_all_seqs.csv", index=False)
# activity (assay data for the original diversity panel)
b = pd.read_excel(f"{legacy}/NylC_Build_Results.xlsx", usecols=['Diversity set ID', 'PA6 L1 Activity 20231117', 'PA66 L2 Activity 20231117'])
b = b.dropna(subset=['Diversity set ID']); b['Diversity set ID'] = b['Diversity set ID'].str.lower()
b = b[b['Diversity set ID'] != 'nylc']  # wild-type NylC has no model; keep nylc-gyaq
b.columns = ['enzyme', 'PA6_L1_Activity_20231117', 'PA66_L2_Activity_20231117']
b[~b.enzyme.isin(extra)].to_csv(f"{out}/nylc_activity.csv", index=False)
