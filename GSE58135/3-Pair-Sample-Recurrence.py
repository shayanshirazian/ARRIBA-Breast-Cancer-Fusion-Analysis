import pandas as pd
import os


# ============================================================
# Paths
# ============================================================

INPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\0. Fusion-13-8-26.xlsx"
)


OUTPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\3. pair_sample_recurrence.csv"
)



# ============================================================
# Read input
# ============================================================

df = pd.read_excel(
    INPUT_FILE
)


print("Total fusion rows:", len(df))

print(df.columns.tolist())



# ============================================================
# Rename columns
# ============================================================

if "#gene1" in df.columns:

    df = df.rename(
        columns={
            "#gene1": "gene1"
        }
    )


if "SRR-ID" in df.columns:

    df = df.rename(
        columns={
            "SRR-ID": "sample_id"
        }
    )



# ============================================================
# Check required columns
# ============================================================

required = [

    "sample_id",
    "gene1",
    "gene2"

]


for col in required:

    if col not in df.columns:

        raise ValueError(
            f"Missing required column: {col}"
        )



# ============================================================
# Clean data
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


df["sample_id"] = (

    df["sample_id"]
    .astype(str)
    .str.strip()

)



# ============================================================
# Create fusion pair
# ============================================================

df["fusion_pair"] = (

    df["gene1"]
    +
    "--"
    +
    df["gene2"]

)



# ============================================================
# Calculate recurrence
# ============================================================

result = (

    df.groupby(
        "fusion_pair"
    )

    .agg(

        gene1=(
            "gene1",
            "first"
        ),

        gene2=(
            "gene2",
            "first"
        ),

        event_count=(
            "fusion_pair",
            "count"
        ),

        sample_count=(
            "sample_id",
            "nunique"
        )

    )

    .reset_index()

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

)



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

print("============================================")

print("PAIR SAMPLE RECURRENCE FINISHED")

print("============================================")


print(
    "Unique fusion pairs:",
    len(result)
)


print()

print(
    result.head(10)
)


print()

print("Output:")

print(OUTPUT_FILE)