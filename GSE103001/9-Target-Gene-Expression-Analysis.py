import pandas as pd


# ============================================================
# Paths
# ============================================================

INPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\Arriba-Downstream-Analysis\GEO2R\GSE103001.top.table.tsv"
)


OUTPUT_FILE = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\Arriba-Downstream-Analysis\output\9.target_gene_expression_analysis.csv"
)



# ============================================================
# Target genes
# ============================================================

target_genes = [

    "TIMM17B",
    "ERAS",
    "TIMM23",
    "TIMM50",
    "TOMM20",
    "ATF5",
    "PIK3CA",
    "AKT1",
    "MTOR",
    "RPS6KB1"

]



# ============================================================
# Read GEO2R result
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)


print(df.columns.tolist())



# ============================================================
# Check Symbol
# ============================================================

if "Symbol" not in df.columns:

    raise ValueError(
        "Symbol column not found"
    )



# ============================================================
# Extract genes
# ============================================================

result = []


for gene in target_genes:


    match = df[
        df["Symbol"]
        .astype(str)
        .str.upper()
        ==
        gene.upper()
    ]


    if len(match) > 0:


        row = match.iloc[0].copy()


        if row["padj"] < 0.05:

            if row["log2FoldChange"] > 0:

                regulation = "Upregulated"

            else:

                regulation = "Downregulated"


        else:

            regulation = "Not significant"



        result.append({

            "Symbol": gene,

            "Detected_in_GEO2R": "Yes",

            "GeneID":
                row["GeneID"],

            "baseMean":
                row["baseMean"],

            "log2FoldChange":
                row["log2FoldChange"],

            "pvalue":
                row["pvalue"],

            "padj":
                row["padj"],

            "Regulation":
                regulation

        })


    else:


        result.append({

            "Symbol": gene,

            "Detected_in_GEO2R": "No",

            "GeneID": None,

            "baseMean": None,

            "log2FoldChange": None,

            "pvalue": None,

            "padj": None,

            "Regulation":
                "Not detected in GEO2R"

        })



# ============================================================
# Save
# ============================================================

out = pd.DataFrame(result)


out.to_csv(
    OUTPUT_FILE,
    index=False
)



print()

print("==============================")

print("Target gene analysis complete")

print("==============================")


print(out)


print()

print("Saved:")

print(OUTPUT_FILE)