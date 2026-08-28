import pandas as pd
import os


# ============================================================
# Input / Output
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
    "4. multi_raw.csv"
)



# ============================================================
# Read Arriba master file
# ============================================================

df = pd.read_excel(INPUT_FILE)

print("Total Arriba rows:", len(df))



# ============================================================
# Detect sample column
# ============================================================

sample_candidates = [
    "sample-id",
    "SRR-ID",
    "SRR_ID",
    "SRR",
    "sample_id",
    "Sample",
    "sample"
]


SAMPLE_COL = None


for col in sample_candidates:
    if col in df.columns:
        SAMPLE_COL = col
        break


if SAMPLE_COL is None:
    raise ValueError(
        "Sample column not found.\n"
        f"Available columns: {list(df.columns)}"
    )


print("Sample column detected:", SAMPLE_COL)



# ============================================================
# Check required columns
# ============================================================

for col in ["gene1", "gene2"]:

    if col not in df.columns:
        raise ValueError(
            f"Required column not found: {col}"
        )



# ============================================================
# Prepare data
# ============================================================

data = df[
    [
        SAMPLE_COL,
        "gene1",
        "gene2"
    ]
].copy()


data.columns = [
    "sample_id",
    "gene1_raw",
    "gene2_raw"
]


data = data.dropna(
    subset=[
        "sample_id",
        "gene1_raw",
        "gene2_raw"
    ]
)



# ============================================================
# Clean spaces only
# Keep raw annotation unchanged
# ============================================================

for col in [
    "sample_id",
    "gene1_raw",
    "gene2_raw"
]:

    data[col] = (
        data[col]
        .astype(str)
        .str.strip()
    )



# ============================================================
# Detect multi-gene annotations
# ============================================================

multi_mask = (
    data["gene1_raw"].str.contains(",", regex=False) |
    data["gene2_raw"].str.contains(",", regex=False)
)


multi = data[multi_mask].copy()



# ============================================================
# Count event frequency
# ============================================================

event_counts = (
    multi
    .groupby(
        [
            "gene1_raw",
            "gene2_raw"
        ]
    )
    .size()
    .reset_index(name="event_count")
)



# ============================================================
# Count independent samples
# ============================================================

sample_counts = (
    multi
    .drop_duplicates(
        subset=[
            "sample_id",
            "gene1_raw",
            "gene2_raw"
        ]
    )
    .groupby(
        [
            "gene1_raw",
            "gene2_raw"
        ]
    )["sample_id"]
    .nunique()
    .reset_index(name="sample_count")
)



# ============================================================
# Merge results
# ============================================================

result = event_counts.merge(
    sample_counts,
    on=[
        "gene1_raw",
        "gene2_raw"
    ],
    how="left"
)



# ============================================================
# Sort
# ============================================================

result = result.sort_values(
    by=[
        "sample_count",
        "event_count",
        "gene1_raw",
        "gene2_raw"
    ],
    ascending=[
        False,
        False,
        True,
        True
    ]
).reset_index(drop=True)



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
print("================================================")
print("MULTI-GENE RAW ANALYSIS FINISHED")
print("================================================")

print("Multi-gene Arriba rows :", len(multi))
print("Unique raw patterns    :", len(result))
print("Output                 :", OUTPUT_FILE)

print()
print(result.head(20).to_string(index=False))