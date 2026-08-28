import pandas as pd


# ===============================
# Paths
# ===============================

INPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001\GEO2R\GSE103001.top.table.tsv"
)


OUTPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\DESeq2\GSE103001_Target_Genes_DESeq2.csv"
)


# ===============================
# Target genes
# ===============================

target_genes = [

    "TIMM17B",
    "ERAS",
    "TIMM23",
    "TIMM50",
    "TOMM20",
    "ATF5",
    "PIK3CA",
    "AKT1",
    "MTOR",
    "RPS6KB1"

]


# ===============================
# Read DESeq2 result
# ===============================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)


print("Columns:")
print(df.columns.tolist())



# ===============================
# Detect gene column
# ===============================

possible_gene_cols = [

    "gene",
    "Gene",
    "symbol",
    "SYMBOL",
    "external_gene_name"

]


gene_col = None


for c in possible_gene_cols:

    if c in df.columns:

        gene_col = c

        break



if gene_col is None:

    raise ValueError(
        "Gene column not found"
    )



# ===============================
# Filter genes
# ===============================

target_df = df[
    df[gene_col]
    .isin(target_genes)
].copy()



# ===============================
# Add regulation
# ===============================

def regulation(row):

    if pd.isna(row["padj"]):

        return "Not significant"


    if row["padj"] < 0.05:

        if row["log2FoldChange"] > 0:

            return "Upregulated"

        else:

            return "Downregulated"


    return "Not significant"



target_df["Regulation"] = target_df.apply(
    regulation,
    axis=1
)



# ===============================
# Save
# ===============================

target_df.to_csv(
    OUTPUT_FILE,
    index=False
)



print()
print("==============================")
print("Target gene analysis finished")
print("==============================")

print()

print(target_df)

print()

print("Saved:")
print(OUTPUT_FILE)