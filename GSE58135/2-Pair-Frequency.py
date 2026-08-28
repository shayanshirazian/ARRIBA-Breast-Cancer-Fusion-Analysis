import pandas as pd
import os


# ============================================================
# Input / Output
# ============================================================

BASE_DIR = r"C:\MyFiles\Bioinfo\ARRIBA\GSE58135\Arriba-Downstream-Analysis"


INPUT_FILE = os.path.join(
    BASE_DIR,
    "Outputs",
    "0. Fusion-13-8-26.xlsx"
)


OUTPUT_FILE = os.path.join(
    BASE_DIR,
    "Outputs",
    "2. pair_frequency.csv"
)



# ============================================================
# Read Arriba Master File
# ============================================================

df = pd.read_excel(INPUT_FILE)


print("Total fusion rows:", len(df))



# ============================================================
# Check columns
# ============================================================

required_columns = [
    "#gene1",
    "gene2"
]


for col in required_columns:

    if col not in df.columns:

        raise ValueError(
            f"Missing required column: {col}"
        )



# ============================================================
# Extract fusion pairs
# ============================================================

pairs = df[
    [
        "#gene1",
        "gene2"
    ]
].copy()


pairs.columns = [
    "gene1",
    "gene2"
]



# ============================================================
# Remove missing values
# ============================================================

pairs = pairs.dropna(
    subset=[
        "gene1",
        "gene2"
    ]
)



# ============================================================
# Clean names
# ============================================================

pairs["gene1"] = (
    pairs["gene1"]
    .astype(str)
    .str.strip()
)


pairs["gene2"] = (
    pairs["gene2"]
    .astype(str)
    .str.strip()
)



pairs = pairs[
    (pairs["gene1"] != "") &
    (pairs["gene2"] != "")
]



# ============================================================
# Count fusion pair frequency
# ============================================================

pair_frequency = (

    pairs
    .groupby(
        [
            "gene1",
            "gene2"
        ]
    )
    .size()
    .reset_index(
        name="event_count"
    )

)



# ============================================================
# Sort
# ============================================================

pair_frequency = pair_frequency.sort_values(

    by=[
        "event_count",
        "gene1",
        "gene2"
    ],

    ascending=[
        False,
        True,
        True
    ]

).reset_index(drop=True)



# ============================================================
# Save
# ============================================================

pair_frequency.to_csv(

    OUTPUT_FILE,

    index=False

)



# ============================================================
# Summary
# ============================================================

print()

print("================================================")

print("PAIR FREQUENCY FINISHED")

print("================================================")


print("Fusion rows analyzed:", len(pairs))

print("Unique fusion pairs:", len(pair_frequency))

print("Output:", OUTPUT_FILE)


print()

print(pair_frequency.head(20))