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
    "5. pair_annotated_summary.csv"
)



# ============================================================
# Read Arriba master file
# ============================================================

df = pd.read_excel(INPUT_FILE)

print("Total Arriba rows:", len(df))



# ============================================================
# Required columns
# ============================================================

required_columns = [
    "sample-id",
    "gene1",
    "gene2",
    "breakpoint1",
    "breakpoint2",
    "type",
    "split_reads1",
    "split_reads2",
    "discordant_mates",
    "confidence",
    "reading_frame",
    "tags",
    "retained_protein_domains"
]


for col in required_columns:

    if col not in df.columns:

        raise ValueError(
            f"Required column not found: {col}"
        )



# ============================================================
# Fill sample IDs
# ============================================================

df["sample-id"] = df["sample-id"].ffill()



# ============================================================
# Clean fields
# ============================================================

df = df.dropna(
    subset=[
        "gene1",
        "gene2",
        "sample-id"
    ]
).copy()



for col in [
    "gene1",
    "gene2",
    "sample-id"
]:

    df[col] = (
        df[col]
        .astype(str)
        .str.strip()
    )



# ============================================================
# Helper function
# ============================================================

def unique_values(series):

    values = []

    for value in series.dropna():

        value = str(value).strip()

        if value in ["", ".", "nan"]:
            continue

        if value not in values:
            values.append(value)

    return " | ".join(values)



# ============================================================
# Confidence ranking
# ============================================================

confidence_rank = {
    "low": 1,
    "medium": 2,
    "high": 3
}


def best_confidence(series):

    values = [
        str(x).strip().lower()
        for x in series.dropna()
        if str(x).strip().lower()
        in confidence_rank
    ]


    if not values:
        return ""


    return max(
        values,
        key=lambda x: confidence_rank[x]
    )



# ============================================================
# Numeric conversion
# ============================================================

numeric_columns = [
    "split_reads1",
    "split_reads2",
    "discordant_mates"
]


for col in numeric_columns:

    df[col] = pd.to_numeric(
        df[col],
        errors="coerce"
    )



# ============================================================
# Build annotated summary
# ============================================================

rows = []


for (gene1, gene2), group in df.groupby(
    [
        "gene1",
        "gene2"
    ],
    sort=False
):


    event_count = len(group)


    sample_count = (
        group["sample-id"]
        .nunique()
    )


    samples = unique_values(
        group["sample-id"]
    )


    support = (
        group["split_reads1"].fillna(0)
        +
        group["split_reads2"].fillna(0)
        +
        group["discordant_mates"].fillna(0)
    )


    rows.append({

        "gene1": gene1,
        "gene2": gene2,

        "event_count": event_count,

        "sample_count": sample_count,

        "samples": samples,

        "best_confidence":
            best_confidence(
                group["confidence"]
            ),

        "max_split_reads1":
            group["split_reads1"].max(),

        "max_split_reads2":
            group["split_reads2"].max(),

        "max_discordant_mates":
            group["discordant_mates"].max(),

        "max_total_read_support":
            support.max(),

        "types":
            unique_values(
                group["type"]
            ),

        "reading_frames":
            unique_values(
                group["reading_frame"]
            ),

        "breakpoint1":
            unique_values(
                group["breakpoint1"]
            ),

        "breakpoint2":
            unique_values(
                group["breakpoint2"]
            ),

        "tags":
            unique_values(
                group["tags"]
            ),

        "retained_protein_domains":
            unique_values(
                group["retained_protein_domains"]
            )
    })



# ============================================================
# Create dataframe
# ============================================================

result = pd.DataFrame(rows)



# ============================================================
# Sort
# ============================================================

result = result.sort_values(
    by=[
        "sample_count",
        "event_count",
        "max_total_read_support"
    ],
    ascending=[
        False,
        False,
        False
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
print("====================================================")
print("PAIR ANNOTATED SUMMARY FINISHED")
print("====================================================")

print("Unique fusion pairs :", len(result))

print("Output:")
print(OUTPUT_FILE)


print()

print(
    result[
        [
            "gene1",
            "gene2",
            "event_count",
            "sample_count",
            "best_confidence",
            "max_total_read_support"
        ]
    ]
    .head(20)
    .to_string(index=False)
)