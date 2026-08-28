import pandas as pd
import gseapy as gp
import os
import numpy as np


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
    "10-1.expression_matrix_symbol_TPM.csv"
)


GMT_FILE = os.path.join(
    BASE_DIR,
    "GeneSets",
    "pathways.gmt"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "10-2.ssGSEA_pathway_scores.csv"
)



# ============================================================
# Read expression matrix
# ============================================================

expr = pd.read_csv(
    INPUT_FILE
)


expr = expr.set_index(
    "Symbol"
)


print("Expression matrix:")
print(expr.shape)



# ============================================================
# log2 transformation
# ============================================================

expr_log = np.log2(expr + 1)




# ============================================================
# UPPERCASE GENE SYMBOLS
# ============================================================

expr.index = expr.index.astype(str).str.upper()



# ============================================================
# ssGSEA
# ============================================================

ssgsea = gp.ssgsea(
    data=expr_log,
    gene_sets=GMT_FILE,
    sample_norm_method="rank",
    outdir=None,
    no_plot=True,
    permutation_num=0,
    min_size=2,
    max_size=500
)



# ============================================================
# Extract scores
# ============================================================

scores = ssgsea.res2d



# ============================================================
# Save
# ============================================================

scores.to_csv(
    OUTPUT_FILE,
    index=False
)



print()
print("==============================")
print("ssGSEA finished")
print("==============================")

print(scores.head())

print()

print("Saved:")
print(OUTPUT_FILE)