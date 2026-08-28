import pandas as pd
import numpy as np
import os
from scipy.stats import mannwhitneyu


# ============================================================
# Paths
# ============================================================

BASE_DIR = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\Arriba-Downstream-Analysis"
)


INPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "10-2.ssGSEA_pathway_scores.csv"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "10-3.pathway_comparison.csv"
)



# ============================================================
# Read ssGSEA result
# ============================================================

df = pd.read_csv(
    INPUT_FILE
)


print(df.head())


# ============================================================
# Convert long format to matrix
# ============================================================

matrix = df.pivot(
    index="Term",
    columns="Name",
    values="NES"
)


print()
print("Pathway matrix:")
print(matrix.shape)



# ============================================================
# Define groups
# ============================================================

normal_samples = [
    f"GSM27523{i}"
    for i in range(50,72)
]


tumor_samples = [
    f"GSM27523{i}"
    for i in range(72,94)
]



# Check samples

normal_samples = [
    x for x in normal_samples
    if x in matrix.columns
]


tumor_samples = [
    x for x in tumor_samples
    if x in matrix.columns
]



print()
print("Normal:", len(normal_samples))
print("Tumor:", len(tumor_samples))



# ============================================================
# Statistical comparison
# ============================================================

results = []


for pathway in matrix.index:


    normal_values = (
        matrix.loc[pathway, normal_samples]
        .values
    )


    tumor_values = (
        matrix.loc[pathway, tumor_samples]
        .values
    )


    stat, pvalue = mannwhitneyu(
        normal_values,
        tumor_values,
        alternative="two-sided"
    )


    results.append({

        "Pathway": pathway,

        "Normal_mean":
            np.mean(normal_values),

        "Tumor_mean":
            np.mean(tumor_values),

        "Difference":
            np.mean(tumor_values)
            -
            np.mean(normal_values),

        "pvalue":
            pvalue

    })



result_df = pd.DataFrame(results)



# Multiple testing correction

from statsmodels.stats.multitest import multipletests


result_df["padj"] = (
    multipletests(
        result_df["pvalue"],
        method="fdr_bh"
    )[1]
)



# ============================================================
# Save
# ============================================================

result_df = result_df.sort_values(
    "padj"
)


result_df.to_csv(
    OUTPUT_FILE,
    index=False
)



print()
print("==============================")
print("Pathway comparison finished")
print("==============================")

print(result_df)


print()
print("Saved:")
print(OUTPUT_FILE)