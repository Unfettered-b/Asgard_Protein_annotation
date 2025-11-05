# faced dependency issues in local runtime thus shifted to colab. This script is discontinued


from Bio import SeqIO
import torch
import esm
import biotite.structure.io as bsio
import pandas as pd
import os
from tqdm import tqdm 

input_fasta = "/home/anirudh/genomes/asCOGs/results/selected/denovo_reps_large.faa"
outfolder = "/home/anirudh/genomes/predicted/"
pdbfolder = os.path.join(outfolder, "pdbs")
os.makedirs(outfolder, exist_ok=True)
os.makedirs(pdbfolder, exist_ok=True)

pdbs = [] # empty list for storing the names of the pdb structures
score = [] # empty list for storing pLDDT values of the predicted structure of proteins
ids = []
seqs = []

pdbase = pd.DataFrame({"IDs": ids,
                        "seqs": seqs,
                        "pdb model" :pdbs,
                        "confidence": score })

pdbcsvfile = os.path.join(outfolder, "prediction_summary.csv")


model = torch.hub.load("facebookresearch/esm:main", "esmfold_v1")

model = model.eval().cuda()

print("Loaded ESMfold model")

for record in tqdm(SeqIO.parse(input_fasta, "fasta"), desc="Predicting Structures"):
    

    print(f"Predicting structure for {record.id}")
    ids.append(record.id)
    seq = str(record.seq)
    seqs.append(seq)

    # ESMfold struggles structure prediction for long proteins (>1024 aa); prediction of long proteins will be done through colabfold
    if len(seq) > 1024:
        print(f"Skipping {record.id} (sequence too long: {len(seq)} aa)")
        pdbs.append(None)
        score.append(None)
        continue
    
    try: 
        with torch.no_grad():
            pdb = model.infer_pdb(seq)

        # structure prediction and storage
        pdb_path = os.path.join(pdbfolder, f"{record.id}.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb)
        print(f"{record.id}.pdb written and saved at {pdb_path}")
        pdbs.append(pdb_path)

        # this will be the pLDDT
        try:
            struct = bsio.load_structure(os.path.join(pdbfolder, f"{record.id}.pdb"), extra_fields=["b_factor"])
            score.append(struct.b_factor.mean())
        except Exception as e:
            print(f"Warning: failed to load structure for {record.id}: {e}")
            score.append(None)

        print("pLDDT for the structure: ", score[-1]) 

        # removing cuda cache to prevent cache build up
        torch.cuda.empty_cache()

        if len(ids) % 10 == 0:
            pdbase = pd.DataFrame({"IDs": ids, "seqs": seqs, "pdb model": pdbs, "confidence": score})
            pdbase.to_csv(pdbcsvfile, index=False)

    
    except Exception as e:
        print(f"Failed for {record.id}: {e}")
        pdbs.append(None)
        score.append(None)



pdbase.to_csv(pdbcsvfile)

print(f"Prediction summary saved to {pdbcsvfile}")


