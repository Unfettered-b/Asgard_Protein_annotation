from Bio import SeqIO
import torch
import biotite.structure.io as bsio
import pandas as pd
import os
from tqdm import tqdm 
import esm
from datetime import datetime
import signal 


# define timeout for structure prediction (in seconds)
TIMEOUT = 420  # 7 minutes



# =====================
# Logging setup
# =====================
LOG_DIR = "/home/anirudh/genomes/predicted/logs/"
os.makedirs(LOG_DIR, exist_ok=True)
LOG_FILE = os.path.join(LOG_DIR, "esm_inference.log")

def log(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {message}"
    # Write to log file
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    # Print above tqdm bar
    tqdm.write(line)



# =====================
# Setting up model run and timeout
# =====================


def timeout_handler(signum, frame):
    raise TimeoutError("Took too long!")

signal.signal(signal.SIGALRM, timeout_handler)


# =====================
# Input / Output paths
# =====================
input_fasta = "/home/anirudh/genomes/asCOGs/results/selected/denovo_reps_large.faa"
outfolder = "/home/anirudh/genomes/predicted/"
pdbfolder = os.path.join(outfolder, "pdbs")
os.makedirs(outfolder, exist_ok=True)
os.makedirs(pdbfolder, exist_ok=True)
pdbcsvfile = os.path.join(outfolder, "prediction_summary.csv")

# =====================
# Load previous predictions if they exist
# =====================
ids, seqs, pdbs, score = [], [], [], []

pdbase = pd.DataFrame({"IDs": ids, "seqs": seqs, "pdb model": pdbs, "confidence": score})
if os.path.exists(pdbcsvfile):
    log(f"Loading existing prediction summary from {pdbcsvfile}")
    pdbase = pd.read_csv(pdbcsvfile)
    ids = pdbase["IDs"].tolist()
    seqs = pdbase["seqs"].tolist()
    pdbs = pdbase["pdb model"].tolist()
    score = pdbase["confidence"].tolist()
    log(f"Loaded {len(ids)} existing predictions")

# =====================
# Counters
# =====================
processed_count = 0
skipped_long_count = 0
failed_count = 0

# =====================
# Load ESMFold model
# =====================
log("Loading ESMFold model...")
model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
log("ESMFold model loaded.")

# =====================
# Structure prediction loop
# =====================
for i, record in enumerate(tqdm(SeqIO.parse(input_fasta, "fasta"), desc="Predicting Structures")):
    if record.id in ids:
        log(f"Skipping {record.id}, already predicted.")
        continue

    processed_count += 1
    seq = str(record.seq)
    ids.append(record.id)
    seqs.append(seq)

    if len(seq) > 1000:
        log(f"Skipping {record.id} (sequence too long: {len(seq)} aa)")
        skipped_long_count += 1
        pdbs.append(None)
        score.append(None)
    else:
        try:
            with torch.no_grad():
                signal.alarm(TIMEOUT)  # Start timeout countdown
                pdb = model.infer_pdb(seq)  # Warm-up run
                
                

            pdb_path = os.path.join(pdbfolder, f"{record.id}.pdb")
            with open(pdb_path, "w") as f:
                f.write(pdb)
            log(f"{record.id}.pdb written and saved at {pdb_path}")
            pdbs.append(pdb_path)

            try:
                struct = bsio.load_structure(pdb_path, extra_fields=["b_factor"])
                score.append(struct.b_factor.mean())
            except Exception as e:
                log(f"Warning: failed to load structure for {record.id}: {e}")
                score.append(None)

            log(f"pLDDT for {record.id}: {score[-1]}") 
            torch.cuda.empty_cache()

        except TimeoutError:
            failed_count += 1
            log(f"Timeout: Prediction for {record.id} exceeded {TIMEOUT} seconds.")
            pdbs.append(None)
            score.append(None)

            torch.cuda.empty_cache()

        except Exception as e:
            failed_count += 1
            log(f"Failed for {record.id}: {e}")
            pdbs.append(None)
            score.append(None)

            torch.cuda.empty_cache()
        
        finally:
            signal.alarm(0)  # Disable timeout
        


    # Save intermediate results and append summary every 10 predictions
    if len(ids) % 10 == 0:
        pdbase = pd.DataFrame({"IDs": ids, "seqs": seqs, "pdb model": pdbs, "confidence": score})
        pdbase.to_csv(pdbcsvfile, index=False)

        summary = (
            f"=== ESMFold Summary after {len(ids)} runs ===\n"
            f"Total IDs processed: {processed_count}, {len(ids)}\n"
            f"Skipped due to long sequences in this session: {skipped_long_count}\n"
            f"Failed runs in this session: {failed_count}\n"
            f"Successful predictions: {processed_count - skipped_long_count - failed_count}\n"
            f"Valid structures (pLDDT > 80): {(lambda lst: sum(1 for x in lst if x is not None and x > 80))(score)}\n"

            f"============================================"
        )
        log(summary)

# =====================
# Final save
# =====================
pdbase = pd.DataFrame({"IDs": ids, "seqs": seqs, "pdb model": pdbs, "confidence": score})
pdbase.to_csv(pdbcsvfile, index=False)
torch.cuda.empty_cache()

log(f"Prediction summary saved to {pdbcsvfile}")
log("ESMFold pipeline completed.")