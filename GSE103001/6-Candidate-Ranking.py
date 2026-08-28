import pandas as pd
import os


# ============================================================
# Input / Output
# ============================================================

BASE_DIR = r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001\Arriba-Downstream-Analysis"


INPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "5. pair_annotated_summary.csv"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "6. candidate_ranking.csv"
)



# ============================================================
# Read annotated pair summary
# ============================================================

df = pd.read_csv(INPUT_FILE)


print("Total fusion pairs:", len(df))



# ============================================================
# Confidence score
# ============================================================

confidence_score = {

    "low": 1,

    "medium": 2,

    "high": 3

}


df["confidence_score"] = (

    df["best_confidence"]
    .astype(str)
    .str.lower()
    .map(confidence_score)
    .fillna(0)
    .astype(int)

)



# ============================================================
# Descriptive flags
# ============================================================

df["recurrent"] = (
    df["sample_count"] > 1
)



df["in_frame"] = (

    df["reading_frames"]
    .fillna("")
    .astype(str)
    .str.contains(
        "in-frame",
        case=False,
        regex=False
    )

)



# ============================================================
# Rank candidates
# ============================================================

df = df.sort_values(

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
# Add ranking number
# ============================================================

df.insert(

    0,

    "candidate_rank",

    range(
        1,
        len(df)+1
    )

)



# ============================================================
# Save
# ============================================================

df.to_csv(

    OUTPUT_FILE,

    index=False

)



# ============================================================
# Summary
# ============================================================

print()

print("============================================")

print("CANDIDATE RANKING FINISHED")

print("============================================")


print("Output:")

print(OUTPUT_FILE)


print()

print("Top 30 candidates:")


columns_to_show = [

    "candidate_rank",

    "gene1",

    "gene2",

    "sample_count",

    "event_count",

    "best_confidence",

    "max_total_read_support",

    "reading_frames"

]


print(

    df[columns_to_show]
    .head(30)
    .to_string(index=False)

)