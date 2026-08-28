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
    "2. pair_gene_frequency.csv"
)



# ============================================================
# Read ARRIBA master file
# ============================================================

df = pd.read_excel(INPUT_FILE)


print("Total fusion rows:", len(df))



# ============================================================
# Check columns
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
# Clean gene names
# ============================================================

df["gene1"] = (
    df["gene1"]
    .astype(str)
    .str.strip()
)


df["gene2"] = (
    df["gene2"]
    .astype(str)
    .str.strip()
)



# ============================================================
# Create fusion pair
# ============================================================

df["fusion_pair"] = (
    df["gene1"]
    + "--" +
    df["gene2"]
)



# ============================================================
# Count fusion pairs
# ============================================================

result = (
    df["fusion_pair"]
    .value_counts()
    .reset_index()
)


result.columns = [
    "fusion_pair",
    "event_count"
]



# ============================================================
# Split genes back
# ============================================================

result["gene1"] = (
    result["fusion_pair"]
    .str.split("--")
    .str[0]
)


result["gene2"] = (
    result["fusion_pair"]
    .str.split("--")
    .str[1]
)



# ============================================================
# Reorder columns
# ============================================================

result = result[
    [
        "gene1",
        "gene2",
        "fusion_pair",
        "event_count"
    ]
]



# ============================================================
# Save
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
print("PAIR GENE FREQUENCY FINISHED")
print("==============================================")

print("Unique fusion pairs:", len(result))

print("Output:")
print(OUTPUT_FILE)

print()

print("Top 20 fusion pairs:")
print(
    result.head(20)
    .to_string(index=False)
)