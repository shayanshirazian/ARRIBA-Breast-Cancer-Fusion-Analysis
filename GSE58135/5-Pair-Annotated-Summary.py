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
    "5. pair_annotated_summary.csv"
)



# ============================================================
# Read Arriba master file
# ============================================================

df = pd.read_excel(INPUT_FILE)

print("Total fusion rows:", len(df))



# ============================================================
# Normalize columns
# ============================================================

rename_map = {}

for col in df.columns:

    clean = str(col).strip().lower()

    if clean == "srr-id":
        rename_map[col] = "sample_id"

    elif clean == "#gene1":
        rename_map[col] = "gene1"

    elif clean == "gene2":
        rename_map[col] = "gene2"


df = df.rename(columns=rename_map)



# ============================================================
# Required columns
# ============================================================

required_columns = [

    "sample_id",
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
            f"Missing required column: {col}"
        )



# ============================================================
# Fill sample IDs
# ============================================================

df["sample_id"] = df["sample_id"].ffill()



# ============================================================
# Clean values
# ============================================================

df = df.dropna(
    subset=[
        "gene1",
        "gene2",
        "sample_id"
    ]
).copy()


for col in [
    "gene1",
    "gene2",
    "sample_id"
]:

    df[col] = (
        df[col]
        .astype(str)
        .str.strip()
    )



# ============================================================
# Helper functions
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



confidence_rank = {

    "low": 1,
    "medium": 2,
    "high": 3

}



def best_confidence(series):

    values = [

        str(x).lower().strip()

        for x in series.dropna()

        if str(x).lower().strip()
        in confidence_rank

    ]


    if not values:

        return ""


    return max(
        values,
        key=lambda x: confidence_rank[x]
    )



# ============================================================
# Numeric columns
# ============================================================

for col in [

    "split_reads1",
    "split_reads2",
    "discordant_mates"

]:

    df[col] = pd.to_numeric(
        df[col],
        errors="coerce"
    )



# ============================================================
# Build summary
# ============================================================

rows = []


for (gene1, gene2), group in df.groupby(
    [
        "gene1",
        "gene2"
    ],
    sort=False
):


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

        "event_count": len(group),

        "sample_count":
            group["sample_id"].nunique(),

        "samples":
            unique_values(
                group["sample_id"]
            ),

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

    ascending=False

).reset_index(drop=True)



# ============================================================
# Save
# ============================================================

result.to_csv(
    OUTPUT_FILE,
    index=False
)



print()

print("================================================")

print("PAIR ANNOTATED SUMMARY FINISHED")

print("================================================")

print("Unique fusion pairs:", len(result))

print("Output:", OUTPUT_FILE)


print()

print(result.head(20).to_string(index=False))