import pandas as pd


# ============================================================
# Paths
# ============================================================

VALIDATION_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\8. validated_intergene_comparison.csv"
)


BREAKPOINT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\9. breakpoint_validation_summary.csv"
)


OUTPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\10. final_validation_summary.csv"
)



# ============================================================
# Read
# ============================================================

validation = pd.read_csv(
    VALIDATION_FILE
)


breakpoint = pd.read_csv(
    BREAKPOINT_FILE
)



print("Gene pair validation:", len(validation))

print("Breakpoint validation:", len(breakpoint))



# ============================================================
# Normalize gene names
# ============================================================

def normalize_pair(g1, g2):

    return tuple(
        sorted(
            [
                str(g1),
                str(g2)
            ]
        )
    )



# ============================================================
# Prepare validation table
# ============================================================

validation["fusion_key"] = validation.apply(

    lambda x:
    normalize_pair(
        x["gene1_GSE103001"],
        x["gene2_GSE103001"]
    ),

    axis=1

)



validation = validation.drop_duplicates(
    subset=["fusion_key"]
)



# ============================================================
# Prepare breakpoint table
# ============================================================

breakpoint["fusion_key"] = breakpoint.apply(

    lambda x:
    normalize_pair(
        x["gene1"],
        x["gene2"]
    ),

    axis=1

)



breakpoint = breakpoint.drop_duplicates(

    subset=[
        "fusion_key"
    ],

    keep="first"

)



# ============================================================
# Merge
# ============================================================

final = validation.merge(

    breakpoint[

        [

            "fusion_key",

            "breakpoint_status",

            "discovery_breakpoint1",

            "validation_breakpoint1",

            "discovery_breakpoint2",

            "validation_breakpoint2"

        ]

    ],

    on="fusion_key",

    how="left"

)



# ============================================================
# Classification
# ============================================================

def classify(row):

    if row["breakpoint_status"] == "Exact_breakpoint_match":

        return "Strong_validation"


    else:

        return "Gene_pair_validation"



final["final_validation_status"] = final.apply(

    classify,

    axis=1

)



# ============================================================
# Remove helper columns
# ============================================================

final = final.drop(

    columns=[
        "fusion_key"
    ]

)



# ============================================================
# Sort
# ============================================================

final = final.sort_values(

    by=[
        "final_validation_status",
        "sample_count_GSE58135",
        "sample_count_GSE103001"
    ],

    ascending=[
        True,
        False,
        False
    ]

)



# ============================================================
# Save
# ============================================================

final.to_csv(

    OUTPUT_FILE,

    index=False

)



print()

print("====================================")

print("FINAL VALIDATION SUMMARY COMPLETE")

print("====================================")


print()

print(
    final["final_validation_status"]
    .value_counts()
)


print()

print("Output:")

print(OUTPUT_FILE)