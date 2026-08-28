import pandas as pd


# ============================================================
# Paths
# ============================================================

SHARED_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\8. validated_intergene_comparison.csv"
)


DISCOVERY_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\Arriba-Downstream-Analysis\output"
    r"\5. pair_annotated_summary.csv"
)


VALIDATION_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\5. pair_annotated_summary.csv"
)


OUTPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135"
    r"\Arriba-Downstream-Analysis\Outputs"
    r"\9. breakpoint_validation_summary.csv"
)



# ============================================================
# Load
# ============================================================

shared = pd.read_csv(SHARED_FILE)

discovery = pd.read_csv(DISCOVERY_FILE)

validation = pd.read_csv(VALIDATION_FILE)



print("Shared:", len(shared))
print("Discovery:", len(discovery))
print("Validation:", len(validation))



# ============================================================
# Rename ARRIBA columns
# ============================================================

for df in [discovery, validation]:

    if "#gene1" in df.columns:

        df.rename(
            columns={
                "#gene1": "gene1"
            },
            inplace=True
        )



# ============================================================
# Cleaning
# ============================================================

for df in [discovery, validation]:

    for c in [
        "gene1",
        "gene2"
    ]:

        df[c] = (
            df[c]
            .astype(str)
            .str.strip()
        )



# ============================================================
# Functions
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



def parse_bp(value):

    if pd.isna(value):
        return set()

    return set(
        str(value)
        .split("|")
    )



def get_gene_bp(row):

    """
    Create gene -> breakpoint mapping
    for both possible orientations
    """

    g1 = row["gene1"]
    g2 = row["gene2"]

    bp1 = parse_bp(
        row["breakpoint1"]
    )

    bp2 = parse_bp(
        row["breakpoint2"]
    )


    return {

        g1: bp1,
        g2: bp2

    }



def compare_breakpoints(d, v):

    d_map = get_gene_bp(d)

    v_map = get_gene_bp(v)


    if set(d_map.keys()) != set(v_map.keys()):

        return False


    for gene in d_map:

        if len(
            d_map[gene]
            &
            v_map[gene]
        ) == 0:

            return False


    return True



# ============================================================
# Prepare shared pairs
# ============================================================

shared["pair_key"] = shared.apply(

    lambda x:
    normalize_pair(
        x["gene1_GSE103001"],
        x["gene2_GSE103001"]
    ),

    axis=1

)



shared = shared.drop_duplicates(
    subset=["pair_key"]
)



# ============================================================
# Validation
# ============================================================

results = []


for _, row in shared.iterrows():


    pair = row["pair_key"]



    d_match = discovery[

        discovery.apply(

            lambda x:

            normalize_pair(
                x["gene1"],
                x["gene2"]
            )

            ==
            pair,

            axis=1

        )

    ]



    v_match = validation[

        validation.apply(

            lambda x:

            normalize_pair(
                x["gene1"],
                x["gene2"]
            )

            ==
            pair,

            axis=1

        )

    ]



    if len(d_match)==0 or len(v_match)==0:

        continue



    status = "Gene_pair_only"

    best_d = d_match.iloc[0]

    best_v = v_match.iloc[0]



    for _, d in d_match.iterrows():

        for _, v in v_match.iterrows():

            if compare_breakpoints(d,v):

                status = "Exact_breakpoint_match"

                best_d = d

                best_v = v

                break


        if status == "Exact_breakpoint_match":

            break



    results.append({

        "gene1": pair[0],

        "gene2": pair[1],

        "discovery_breakpoint1":
            best_d["breakpoint1"],

        "validation_breakpoint1":
            best_v["breakpoint1"],

        "discovery_breakpoint2":
            best_d["breakpoint2"],

        "validation_breakpoint2":
            best_v["breakpoint2"],

        "breakpoint_status":
            status

    })



# ============================================================
# Final cleanup
# ============================================================

result = pd.DataFrame(results)



result = result.drop_duplicates(

    subset=[
        "gene1",
        "gene2"
    ]

)



result.to_csv(

    OUTPUT_FILE,

    index=False

)



print()
print("================================")
print("BREAKPOINT VALIDATION COMPLETE")
print("================================")


print()

print(
    result["breakpoint_status"]
    .value_counts()
)


print()

print(
    "Unique fusion pairs:",
    len(result)
)


print()

print(
    "Output:",
    OUTPUT_FILE
)