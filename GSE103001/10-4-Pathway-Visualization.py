import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


# ============================================================
# Paths
# ============================================================

BASE_DIR = (
    r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001"
    r"\Arriba-Downstream-Analysis"
)


INPUT_SCORE = os.path.join(
    BASE_DIR,
    "output",
    "10-2.ssGSEA_pathway_scores.csv"
)


INPUT_COMPARE = os.path.join(
    BASE_DIR,
    "output",
    "10-3.pathway_comparison.csv"
)


OUT_SUMMARY = os.path.join(
    BASE_DIR,
    "output",
    "10-4.pathway_summary.csv"
)


OUT_BOX = os.path.join(
    BASE_DIR,
    "output",
    "10-4.pathway_boxplots.png"
)


OUT_HEATMAP = os.path.join(
    BASE_DIR,
    "output",
    "10-4.pathway_heatmap.png"
)



# ============================================================
# Read data
# ============================================================

scores = pd.read_csv(
    INPUT_SCORE
)


comparison = pd.read_csv(
    INPUT_COMPARE
)



# ============================================================
# Add group information
# ============================================================

normal_samples = [
    f"GSM27523{i}"
    for i in range(50,72)
]


tumor_samples = [
    f"GSM27523{i}"
    for i in range(72,94)
]


def group(sample):

    if sample in normal_samples:
        return "Normal"

    elif sample in tumor_samples:
        return "Tumor"

    else:
        return "Unknown"



scores["Group"] = scores["Name"].apply(group)



# ============================================================
# Save summary
# ============================================================

summary = comparison.copy()


summary["Significance"] = np.where(
    summary["padj"] < 0.05,
    "Significant",
    "Not significant"
)


summary.to_csv(
    OUT_SUMMARY,
    index=False
)



# ============================================================
# Boxplots
# ============================================================

plt.figure(
    figsize=(12,8)
)


sns.boxplot(
    data=scores,
    x="Term",
    y="NES",
    hue="Group"
)


plt.xticks(
    rotation=45,
    ha="right"
)


plt.title(
    "ssGSEA Pathway Activity: Tumor vs Normal"
)


plt.tight_layout()


plt.savefig(
    OUT_BOX,
    dpi=300
)


plt.close()



# ============================================================
# Heatmap
# ============================================================

matrix = scores.pivot(
    index="Term",
    columns="Name",
    values="NES"
)


plt.figure(
    figsize=(14,6)
)


sns.heatmap(
    matrix,
    cmap="RdBu_r",
    center=matrix.values.mean()
)


plt.title(
    "ssGSEA Pathway Score Heatmap"
)


plt.tight_layout()


plt.savefig(
    OUT_HEATMAP,
    dpi=300
)


plt.close()



print()
print("==============================")
print("Pathway visualization finished")
print("==============================")

print("Saved:")
print(OUT_SUMMARY)
print(OUT_BOX)
print(OUT_HEATMAP)