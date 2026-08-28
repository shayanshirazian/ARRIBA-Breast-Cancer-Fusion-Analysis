import pandas as pd
import os


# ============================================================
# Paths
# ============================================================

BASE_DIR = r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001\Arriba-Downstream-Analysis"

INPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "0. Fusion-21-8-26.xlsx"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "1. single_gene_frequency.csv"
)



# ============================================================
# Read ARRIBA master file
# ============================================================

df = pd.read_excel(INPUT_FILE)

print("Total fusion rows:", len(df))
print(df.columns.tolist())

# ============================================================
# Check required columns
# ============================================================

required_columns = [
    "gene1",
    "gene2"
]


for col in required_columns:

    if col not in df.columns:

        raise ValueError(
            f"Missing required column: {col}"
        )



# ============================================================
# Extract genes from ARRIBA output
# ============================================================

gene1 = (
    df["gene1"]
    .dropna()
    .astype(str)
    .str.strip()
)


gene2 = (
    df["gene2"]
    .dropna()
    .astype(str)
    .str.strip()
)



# ============================================================
# Count gene frequency
# ============================================================

gene1_counts = gene1.value_counts()

gene2_counts = gene2.value_counts()



# ============================================================
# Combine all genes
# ============================================================

all_genes = sorted(
    set(gene1_counts.index) |
    set(gene2_counts.index)
)



result = pd.DataFrame(
    {
        "GENE": all_genes
    }
)



result["gene1_count"] = (
    result["GENE"]
    .map(gene1_counts)
    .fillna(0)
    .astype(int)
)



result["gene2_count"] = (
    result["GENE"]
    .map(gene2_counts)
    .fillna(0)
    .astype(int)
)



result["total_count"] = (
    result["gene1_count"] +
    result["gene2_count"]
)



# ============================================================
# Sort frequency
# ============================================================

result = result.sort_values(
    by=[
        "total_count",
        "gene1_count",
        "gene2_count"
    ],
    ascending=False
).reset_index(drop=True)



# ============================================================
# Save output
# ============================================================

result.to_csv(
    OUTPUT_FILE,
    index=False
)



# ============================================================
# Summary
# ============================================================

print()
print("==============================================")
print("Single Gene Frequency Finished")
print("==============================================")

print("Unique genes:", len(result))

print("Output:")
print(OUTPUT_FILE)

print()

print("Top 20 genes:")
print(
    result.head(20)
    .to_string(index=False)
)