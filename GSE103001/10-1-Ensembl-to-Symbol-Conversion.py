import pandas as pd
import mygene
import os


# ============================================================
# Paths
# ============================================================

BASE_DIR = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\Arriba-Downstream-Analysis"
)


INPUT_FILE = os.path.join(
    BASE_DIR,
    "Expression",
    "GSE103001_norm_counts_TPM_GRCh38.p13_NCBI.tsv"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "10-1.expression_matrix_symbol_TPM.csv"
)



# ============================================================
# Read TPM matrix
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)


print("Original shape:")
print(df.shape)


print(df.head())



# ============================================================
# Rename GeneID column
# ============================================================

gene_id_col = df.columns[0]


df = df.rename(
    columns={
        gene_id_col: "GeneID"
    }
)



# ============================================================
# Clean GeneID
# ============================================================

df["GeneID_clean"] = (
    df["GeneID"]
    .astype(str)
)



# ============================================================
# Convert NCBI GeneID to Gene Symbol
# ============================================================

print()
print("Converting GeneID to Symbol...")


mg = mygene.MyGeneInfo()


mapping = mg.querymany(
    df["GeneID_clean"].tolist(),
    scopes="entrezgene",
    fields="symbol",
    species="human"
)



mapping_df = pd.DataFrame(mapping)



print("Mapping result:")
print(mapping_df.head())



# ============================================================
# Keep successful mappings
# ============================================================

mapping_df = mapping_df[
    mapping_df["symbol"].notna()
]


mapping_df = mapping_df[
    [
        "query",
        "symbol"
    ]
]


mapping_df = mapping_df.rename(
    columns={
        "query": "GeneID_clean",
        "symbol": "Symbol"
    }
)



# ============================================================
# Merge annotation
# ============================================================

df = df.merge(
    mapping_df,
    on="GeneID_clean",
    how="left"
)



# ============================================================
# Remove unmapped genes
# ============================================================

before = len(df)


df = df[
    df["Symbol"].notna()
]


after = len(df)


print()
print("Mapped genes:")
print(after, "/", before)



# ============================================================
# Sample columns
# ============================================================

sample_cols = [
    c for c in df.columns
    if c.startswith("GSM")
]



# ============================================================
# Remove duplicated symbols
# Keep row with highest total TPM
# ============================================================

df["Total_TPM"] = (
    df[sample_cols]
    .sum(axis=1)
)


df = (
    df
    .sort_values(
        by="Total_TPM",
        ascending=False
    )
    .drop_duplicates(
        subset="Symbol"
    )
)



# ============================================================
# Final expression matrix
# ============================================================

final_df = df[
    [
        "Symbol"
    ]
    +
    sample_cols
]



# ============================================================
# Save
# ============================================================

final_df.to_csv(
    OUTPUT_FILE,
    index=False
)



print()
print("==============================")
print("GeneID to Symbol conversion finished")
print("==============================")

print()

print("Final matrix shape:")
print(final_df.shape)

print()

print(final_df.head())


print()

print("Saved:")
print(OUTPUT_FILE)