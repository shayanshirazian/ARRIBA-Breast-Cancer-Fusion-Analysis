import pandas as pd
import os


# ============================================================
# Input / Output
# ============================================================

BASE_DIR = r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001\Arriba-Downstream-Analysis"


INPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "6. candidate_ranking.csv"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "7. candidate_ranking_intergene.csv"
)



# ============================================================
# Read candidate ranking
# ============================================================

df = pd.read_csv(INPUT_FILE)


print("Total candidate pairs:", len(df))



# ============================================================
# Required columns
# ============================================================

required_columns = [

    "gene1",
    "gene2",
    "sample_count",
    "event_count",
    "best_confidence",
    "max_total_read_support"

]


for col in required_columns:

    if col not in df.columns:

        raise ValueError(
            f"Required column not found: {col}"
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
# Keep true inter-gene fusions
# ============================================================

intergene = df[

    (df["gene1"].str.upper() != df["gene2"].str.upper())

    &

    (~df["gene1"].str.contains(",", regex=False))

    &

    (~df["gene2"].str.contains(",", regex=False))

].copy()



# ============================================================
# Confidence score
# ============================================================

confidence_score = {

    "low": 1,
    "medium": 2,
    "high": 3

}


intergene["confidence_score"] = (

    intergene["best_confidence"]
    .astype(str)
    .str.lower()
    .map(confidence_score)
    .fillna(0)
    .astype(int)

)



# ============================================================
# Re-rank
# ============================================================

intergene = intergene.sort_values(

    by=[

        "sample_count",

        "confidence_score",

        "max_total_read_support",

        "event_count"

    ],

    ascending=[

        False,
        False,
        False,
        False

    ]

).reset_index(drop=True)



# ============================================================
# Replace rank
# ============================================================

if "candidate_rank" in intergene.columns:

    intergene = intergene.drop(
        columns=["candidate_rank"]
    )


intergene.insert(

    0,

    "candidate_rank",

    range(
        1,
        len(intergene)+1
    )

)



# ============================================================
# Save
# ============================================================

intergene.to_csv(

    OUTPUT_FILE,

    index=False

)



# ============================================================
# Summary
# ============================================================

print()

print("================================================")

print("INTER-GENE CANDIDATE RANKING FINISHED")

print("================================================")


print("Original pairs       :", len(df))

print("Inter-gene candidates:", len(intergene))

print("Output               :", OUTPUT_FILE)


print()

print("Top 30 inter-gene fusion candidates:")


show_columns = [

    "candidate_rank",

    "gene1",

    "gene2",

    "sample_count",

    "event_count",

    "best_confidence",

    "max_total_read_support"

]


if "reading_frames" in intergene.columns:

    show_columns.append(
        "reading_frames"
    )


print(

    intergene[show_columns]
    .head(30)
    .to_string(index=False)

)