import os
import pandas as pd

# ============================================================
# Input / Output
# ============================================================

DISCOVERY_FILES = [
    "../ARRIBA output whole results.xlsx",
    "../ARRIBA main output (without inc-dec).xlsx"
]

VALIDATION_FILE = "../GSE58135 Fusion-13-8-26.xlsx"

OUTPUT_FILE = "../GSE58135_validated_intergene_comparison.csv"


# ============================================================
# Fusion pairs already detected in BOTH datasets
# ============================================================

TARGET_PAIRS = [
    ("THBS1", "RP11-624L4.1"),
    ("ANKRD62P1-PARP4P3", "KB-7G2.8"),
    ("MZB1", "DNAJC18"),
]


# ============================================================
# Find discovery file
# ============================================================

DISCOVERY_FILE = None

for file in DISCOVERY_FILES:
    if os.path.exists(file):
        DISCOVERY_FILE = file
        break

if DISCOVERY_FILE is None:
    raise FileNotFoundError(
        "Discovery Arriba Excel file was not found."
    )

print("Discovery file :", DISCOVERY_FILE)
print("Validation file:", VALIDATION_FILE)


# ============================================================
# Column normalization
# ============================================================

def normalize_columns(df):

    rename_map = {}

    for col in df.columns:

        clean = str(col).strip().lower().replace(" ", "_")

        if clean in ["#gene1", "gene1"]:
            rename_map[col] = "gene1"

        elif clean == "gene2":
            rename_map[col] = "gene2"

        elif clean in ["srr-id", "srr_id", "srr", "sample_id", "sample"]:
            rename_map[col] = "sample_id"

        elif clean == "breakpoint1":
            rename_map[col] = "breakpoint1"

        elif clean == "breakpoint2":
            rename_map[col] = "breakpoint2"

        elif clean == "type":
            rename_map[col] = "type"

        elif clean == "confidence":
            rename_map[col] = "confidence"

        elif clean == "reading_frame":
            rename_map[col] = "reading_frame"

        elif clean == "split_reads1":
            rename_map[col] = "split_reads1"

        elif clean == "split_reads2":
            rename_map[col] = "split_reads2"

        elif clean == "discordant_mates":
            rename_map[col] = "discordant_mates"

    return df.rename(columns=rename_map)


# ============================================================
# Read Excel:
# collect all sheets containing gene1/gene2
# ============================================================

def read_arriba_excel(path):

    sheets = pd.read_excel(
        path,
        sheet_name=None
    )

    usable = []

    for sheet_name, df in sheets.items():

        df = normalize_columns(df)

        if "gene1" in df.columns and "gene2" in df.columns:

            df = df.copy()
            df["source_sheet"] = sheet_name

            usable.append(df)

    if not usable:
        raise ValueError(
            f"No sheet containing gene1/gene2 found in {path}"
        )

    df = pd.concat(
        usable,
        ignore_index=True
    )

    # Fill sample IDs if the Excel visually writes SRR only once
    if "sample_id" in df.columns:
        df["sample_id"] = df["sample_id"].ffill()

    # Clean gene names
    df["gene1"] = df["gene1"].astype(str).str.strip()
    df["gene2"] = df["gene2"].astype(str).str.strip()

    # Numeric evidence columns
    for col in [
        "split_reads1",
        "split_reads2",
        "discordant_mates"
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col],
                errors="coerce"
            )

    return df


# ============================================================
# Load datasets
# ============================================================

discovery = read_arriba_excel(DISCOVERY_FILE)
validation = read_arriba_excel(VALIDATION_FILE)

print()
print("Discovery Arriba rows :", len(discovery))
print("Validation Arriba rows:", len(validation))


# ============================================================
# Helpers
# ============================================================

def unique_values(series):

    values = []

    for value in series.dropna():

        value = str(value).strip()

        if value in ["", ".", "nan"]:
            continue

        if value not in values:
            values.append(value)

    return " | ".join(values)


confidence_rank = {
    "low": 1,
    "medium": 2,
    "high": 3
}


def best_confidence(series):

    values = [
        str(x).strip().lower()
        for x in series.dropna()
        if str(x).strip().lower() in confidence_rank
    ]

    if not values:
        return ""

    return max(
        values,
        key=lambda x: confidence_rank[x]
    )


# ============================================================
# Extract one target fusion
#
# Both orientations are accepted.
#
# Breakpoints are normalized to:
# candidate_gene1 -> candidate_gene2
# ============================================================

def extract_pair(df, candidate_gene1, candidate_gene2):

    g1 = candidate_gene1.upper()
    g2 = candidate_gene2.upper()

    direct = df[
        (df["gene1"].str.upper() == g1) &
        (df["gene2"].str.upper() == g2)
    ].copy()

    reverse = df[
        (df["gene1"].str.upper() == g2) &
        (df["gene2"].str.upper() == g1)
    ].copy()

    direct["orientation"] = "direct"
    reverse["orientation"] = "reverse"

    matches = pd.concat(
        [direct, reverse],
        ignore_index=True
    )

    if len(matches) == 0:
        return matches

    # --------------------------------------------------------
    # Normalize breakpoints to candidate gene order
    # --------------------------------------------------------

    normalized_bp1 = []
    normalized_bp2 = []

    for _, row in matches.iterrows():

        orientation = row["orientation"]

        bp1 = row.get("breakpoint1", "")
        bp2 = row.get("breakpoint2", "")

        if orientation == "direct":
            normalized_bp1.append(bp1)
            normalized_bp2.append(bp2)

        else:
            normalized_bp1.append(bp2)
            normalized_bp2.append(bp1)

    matches["normalized_breakpoint_gene1"] = normalized_bp1
    matches["normalized_breakpoint_gene2"] = normalized_bp2

    matches["normalized_breakpoint_pair"] = (
        matches["normalized_breakpoint_gene1"].astype(str)
        + " -> " +
        matches["normalized_breakpoint_gene2"].astype(str)
    )

    # --------------------------------------------------------
    # Total read support for each event
    # --------------------------------------------------------

    support = pd.Series(
        0,
        index=matches.index,
        dtype=float
    )

    for col in [
        "split_reads1",
        "split_reads2",
        "discordant_mates"
    ]:

        if col in matches.columns:
            support += matches[col].fillna(0)

    matches["total_read_support"] = support

    return matches


# ============================================================
# Summarize a cohort
# ============================================================

def summarize(matches, prefix):

    if len(matches) == 0:

        return {
            f"{prefix}_event_count": 0,
            f"{prefix}_sample_count": 0,
            f"{prefix}_samples": "",
            f"{prefix}_best_confidence": "",
            f"{prefix}_max_read_support": 0,
            f"{prefix}_types": "",
            f"{prefix}_reading_frames": "",
            f"{prefix}_breakpoint_pairs": ""
        }

    if "sample_id" in matches.columns:

        sample_count = matches["sample_id"].dropna().nunique()
        samples = unique_values(matches["sample_id"])

    else:

        sample_count = ""
        samples = ""

    return {

        f"{prefix}_event_count":
            len(matches),

        f"{prefix}_sample_count":
            sample_count,

        f"{prefix}_samples":
            samples,

        f"{prefix}_best_confidence":
            best_confidence(matches["confidence"])
            if "confidence" in matches.columns
            else "",

        f"{prefix}_max_read_support":
            matches["total_read_support"].max(),

        f"{prefix}_types":
            unique_values(matches["type"])
            if "type" in matches.columns
            else "",

        f"{prefix}_reading_frames":
            unique_values(matches["reading_frame"])
            if "reading_frame" in matches.columns
            else "",

        f"{prefix}_breakpoint_pairs":
            unique_values(
                matches["normalized_breakpoint_pair"]
            )
    }


# ============================================================
# Compare discovery vs validation
# ============================================================

results = []

for gene1, gene2 in TARGET_PAIRS:

    discovery_matches = extract_pair(
        discovery,
        gene1,
        gene2
    )

    validation_matches = extract_pair(
        validation,
        gene1,
        gene2
    )

    discovery_summary = summarize(
        discovery_matches,
        "discovery"
    )

    validation_summary = summarize(
        validation_matches,
        "validation"
    )

    # --------------------------------------------------------
    # Exact breakpoint overlap
    # --------------------------------------------------------

    discovery_breakpoints = set(
        discovery_matches.get(
            "normalized_breakpoint_pair",
            pd.Series(dtype=str)
        )
        .dropna()
        .astype(str)
    )

    validation_breakpoints = set(
        validation_matches.get(
            "normalized_breakpoint_pair",
            pd.Series(dtype=str)
        )
        .dropna()
        .astype(str)
    )

    breakpoint_overlap = (
        discovery_breakpoints &
        validation_breakpoints
    )

    # --------------------------------------------------------
    # Type overlap
    # --------------------------------------------------------

    discovery_types = set(
        discovery_matches.get(
            "type",
            pd.Series(dtype=str)
        )
        .dropna()
        .astype(str)
    )

    validation_types = set(
        validation_matches.get(
            "type",
            pd.Series(dtype=str)
        )
        .dropna()
        .astype(str)
    )

    type_overlap = (
        discovery_types &
        validation_types
    )

    # --------------------------------------------------------
    # Reading-frame overlap
    # --------------------------------------------------------

    discovery_frames = set(
        discovery_matches.get(
            "reading_frame",
            pd.Series(dtype=str)
        )
        .dropna()
        .astype(str)
    )

    validation_frames = set(
        validation_matches.get(
            "reading_frame",
            pd.Series(dtype=str)
        )
        .dropna()
        .astype(str)
    )

    frame_overlap = (
        discovery_frames &
        validation_frames
    )

    row = {
        "gene1": gene1,
        "gene2": gene2,

        **discovery_summary,
        **validation_summary,

        "exact_breakpoint_overlap":
            len(breakpoint_overlap) > 0,

        "shared_breakpoint_pairs":
            " | ".join(sorted(breakpoint_overlap)),

        "fusion_type_overlap":
            len(type_overlap) > 0,

        "shared_fusion_types":
            " | ".join(sorted(type_overlap)),

        "reading_frame_overlap":
            len(frame_overlap) > 0,

        "shared_reading_frames":
            " | ".join(sorted(frame_overlap))
    }

    results.append(row)


# ============================================================
# Save comparison
# ============================================================

result = pd.DataFrame(results)

result.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# Terminal summary
# ============================================================

print()
print("========================================================")
print("VALIDATED INTER-GENE COMPARISON FINISHED")
print("========================================================")

print("Output:", OUTPUT_FILE)

print()

columns_to_show = [
    "gene1",
    "gene2",
    "discovery_event_count",
    "validation_event_count",
    "validation_sample_count",
    "discovery_best_confidence",
    "validation_best_confidence",
    "exact_breakpoint_overlap",
    "fusion_type_overlap",
    "reading_frame_overlap"
]

print(
    result[columns_to_show]
    .to_string(index=False)
)