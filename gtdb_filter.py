import pandas as pd
import numpy as np

archaeal_metadata = "ar122_metadata.tsv"

metadata = pd.read_csv(archaeal_metadata, sep="\t")

metadata['genus'] = metadata['gtdb_taxonomy'].str.extract(r"g__([^;]+)")

metadata['GQS'] = metadata['checkm_completeness'] - 5*metadata['checkm_contamination']

metadata = metadata.dropna(subset=['GQS', 'genus'])


representatives = (
    metadata.sort_values("GQS", ascending=False)
    .groupby('genus')
    .apply(lambda x:x.sample(1) if len(x)>1 and x['GQS'].nunique() == 1 else x.head(1))
    .reset_index(drop=True)
)

mask = (
    representatives["gtdb_taxonomy"].str.contains("p__Thermoproteota")
    | representatives["gtdb_taxonomy"].str.contains("c__Methanobacteria_B")
    | representatives["gtdb_taxonomy"].str.contains("p__Hadarchaeota")
)
mask &= ~representatives["gtdb_taxonomy"].str.contains("c__Korarchaeia")
outgroups = representatives[mask]



representatives.to_csv('archaea_genus_representatives.tsv', sep="\t", index=False)
