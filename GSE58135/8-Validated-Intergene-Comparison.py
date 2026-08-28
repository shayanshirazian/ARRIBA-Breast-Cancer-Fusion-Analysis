import pandas as pd


# ============================================================
# Paths
# ============================================================

DISCOVERY_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\Arriba-Downstream-Analysis\output"
    r"\7. candidate_ranking_intergene.csv"
)


VALIDATION_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\7. candidate_ranking_intergene.csv"
)


OUTPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\8. validated_intergene_comparison.csv"
)



# ============================================================
# Read files
# ============================================================

discovery = pd.read_csv(
    DISCOVERY_FILE
)


validation = pd.read_csv(
    VALIDATION_FILE
)



print("Discovery rows:", len(discovery))

print("Validation rows:", len(validation))



# ============================================================
# Required columns
# ============================================================

required = [

    "gene1",
    "gene2",
    "sample_count",
    "event_count",
    "best_confidence"

]


for col in required:

    if col not in discovery.columns:

        raise ValueError(
            f"Missing discovery column: {col}"
        )


    if col not in validation.columns:

        raise ValueError(
            f"Missing validation column: {col}"
        )



# ============================================================
# Clean gene names
# ============================================================

for df in [
    discovery,
    validation
]:

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
# Normalize fusion orientation
# ============================================================

def normalize_pair(row):

    return "--".join(
        sorted(
            [
                row["gene1"],
                row["gene2"]
            ]
        )
    )


discovery["fusion_pair"] = (
    discovery.apply(
        normalize_pair,
        axis=1
    )
)


validation["fusion_pair"] = (
    validation.apply(
        normalize_pair,
        axis=1
    )
)



# ============================================================
# Find shared fusion pairs
# ============================================================

shared = discovery.merge(

    validation,

    on="fusion_pair",

    suffixes=(
        "_GSE103001",
        "_GSE58135"
    )

)



# ============================================================
# Select output columns
# ============================================================

result = shared[

    [

        "gene1_GSE103001",
        "gene2_GSE103001",

        "sample_count_GSE103001",
        "sample_count_GSE58135",

        "event_count_GSE103001",
        "event_count_GSE58135",

        "best_confidence_GSE103001",
        "best_confidence_GSE58135"

    ]

].copy()



result["validation_status"] = (
    "Gene_pair_validated"
)



# ============================================================
# Sort
# ============================================================

result = result.sort_values(

    by=[
        "sample_count_GSE58135",
        "sample_count_GSE103001"
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



print()

print("============================================")

print("VALIDATED INTERGENE COMPARISON FINISHED")

print("============================================")

print(
    "Shared fusion pairs:",
    len(result)
)

print()

print("Output:")

print(OUTPUT_FILE)