import pandas as pd
import os


# ============================================================
# Input / Output
# ============================================================

BASE_DIR = r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135\Arriba-Downstream-Analysis"


INPUT_FILE = os.path.join(
    BASE_DIR,
    "Outputs",
    "0. Fusion-13-8-26.xlsx"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "Outputs",
    "4. multi_raw.csv"
)



# ============================================================
# Read Arriba master file
# ============================================================

df = pd.read_excel(INPUT_FILE)


print("Total fusion rows:", len(df))



# ============================================================
# Detect sample column
# ============================================================

sample_candidates = [
    "SRR-ID",
    "sample-id",
    "sample_id",
    "sample"
]


SAMPLE_COL = None


for col in sample_candidates:

    if col in df.columns:

        SAMPLE_COL = col
        break


if SAMPLE_COL is None:

    raise ValueError(
        "Sample column not found"
    )


print("Sample column:", SAMPLE_COL)



# ============================================================
# Required columns
# ============================================================

required_columns = [
    "#gene1",
    "gene2"
]


for col in required_columns:

    if col not in df.columns:

        raise ValueError(
            f"Missing required column: {col}"
        )



# ============================================================
# Prepare data
# ============================================================

data = df[
    [
        SAMPLE_COL,
        "#gene1",
        "gene2"
    ]
].copy()


data.columns = [
    "sample_id",
    "gene1_raw",
    "gene2_raw"
]



# Fill sample IDs if merged cells exist

data["sample_id"] = data["sample_id"].ffill()



# Remove missing

data = data.dropna(
    subset=[
        "gene1_raw",
        "gene2_raw",
        "sample_id"
    ]
)



# Clean

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
# Select multi-gene annotations
# ============================================================

multi = data[
    data["gene1_raw"].str.contains(",", regex=False)
    |
    data["gene2_raw"].str.contains(",", regex=False)
].copy()



# ============================================================
# Event count
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
    .reset_index(
        name="event_count"
    )

)



# ============================================================
# Sample recurrence
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
    .reset_index(
        name="sample_count"
    )

)



# ============================================================
# Merge
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
        "event_count"
    ],

    ascending=False

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

print("==============================================")

print("MULTI RAW FINISHED")

print("==============================================")


print("Multi-gene rows:", len(multi))

print("Unique patterns:", len(result))

print("Output:", OUTPUT_FILE)


print()

print(result.head(20))