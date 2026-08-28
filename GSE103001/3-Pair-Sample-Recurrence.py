import pandas as pd
import os


# ============================================================
# Paths
# ============================================================

BASE_DIR = r"C:\MyFiles\Bioinfo\ARRIBA\GSE103001\Arriba-Downstream-Analysis"


INPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "0. Fusion-21-8-26.xlsx"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "output",
    "3. pair_sample_recurrence.csv"
)



# ============================================================
# Read ARRIBA master file
# ============================================================

df = pd.read_excel(INPUT_FILE)


print("Total fusion rows:", len(df))



# ============================================================
# Check columns
# ============================================================

required_columns = [
    "gene1",
    "gene2",
    "sample-id"
]


for col in required_columns:

    if col not in df.columns:

        raise ValueError(
            f"Missing required column: {col}"
        )



# ============================================================
# Clean columns
# ============================================================

for col in ["gene1", "gene2", "sample-id"]:

    df[col] = (
        df[col]
        .astype(str)
        .str.strip()
    )



# ============================================================
# Count event and independent samples
# ============================================================

result = (

    df
    .groupby(
        [
            "gene1",
            "gene2"
        ]
    )
    .agg(
        event_count=("gene1", "size"),
        sample_count=("sample-id", "nunique")
    )
    .reset_index()

)



# ============================================================
# Sort recurrence
# ============================================================

result = result.sort_values(
    by=[
        "sample_count",
        "event_count"
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



# ============================================================
# Summary
# ============================================================

print()
print("==============================================")
print("PAIR SAMPLE RECURRENCE FINISHED")
print("==============================================")

print("Fusion pairs:", len(result))

print("Output:")
print(OUTPUT_FILE)

print()

print("Top 20 recurrent fusion pairs:")
print(
    result.head(20)
    .to_string(index=False)
)