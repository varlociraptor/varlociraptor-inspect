import altair as alt
import pandas as pd
import re


def phred_to_prob(phred_value):
    """Convert PHRED score to probability"""
    if phred_value is None:
        return None
    return 10 ** (-phred_value / 10)


def visualize_event_probabilities(record):
    """Visualize event probabilities from INFO column (PROB_* fields)"""
    prob_data = []

    for key, value in record.info.items():
        if key.startswith("PROB_"):
            event_name = key.replace("PROB_", "")

            # Handle tuple/list outside phred_to_prob
            if len(value) != 1:
                raise ValueError(
                    f"Unexpected number of values for {key} (expected: 1, found: {len(value)})"
                )
            value = value[0]

            probability = 0.0 if value == float("inf") else phred_to_prob(value)

            prob_data.append({"Event": event_name, "Probability": probability})

    df = pd.DataFrame(prob_data)

    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            alt.X("Event:N"),
            alt.Y("Probability:Q"),
            tooltip=[
                "Event",
                alt.Tooltip("Probability:Q", format=".6f"),
            ],
        )
        .properties(title="Event Probabilities", width=400, height=300)
    )


def visualize_allele_frequency_distribution(record, sample_name):
    """Visualize allele frequency distribution (AFD field)"""
    sample = record.samples[sample_name]

    # Assume these are always tuples/lists
    af_ml = sample["AF"][0]
    afd_entries = sample["AFD"]

    afd_data = []

    # Process all distribution points
    for entry in afd_entries:
        for part in entry.split(","):
            if "=" in part:
                freq, phred = part.split("=")
                freq = float(freq)
                prob = phred_to_prob(float(phred))

                afd_data.append(
                    {
                        "Allele Frequency": freq,
                        "Probability": prob,
                        "Type": "Distribution",
                    }
                )

    # Explicitly add ML estimate
    # Find probability at ML frequency, or use 1.0 if not found
    ml_prob = next(
        (
            d["Probability"]
            for d in afd_data
            if abs(d["Allele Frequency"] - af_ml) < 0.001
        ),
        1.0,
    )

    afd_data.append(
        {
            "Allele Frequency": af_ml,
            "Probability": ml_prob,
            "Type": "ML Estimate",
        }
    )

    df = pd.DataFrame(afd_data)

    return (
        alt.Chart(df)
        .mark_circle()
        .encode(
            alt.X("Allele Frequency:Q"),
            alt.Y("Probability:Q", axis=None),
            alt.Color(
                "Type:N",
                scale=alt.Scale(
                    domain=["Distribution", "ML Estimate"],
                    range=["blue", "red"],
                ),
            ),
            alt.Size(
                "Type:N",
                scale=alt.Scale(
                    domain=["Distribution", "ML Estimate"],
                    range=[60, 100],
                ),
                legend=None,
            ),
            alt.Opacity(
                "Type:N",
                scale=alt.Scale(
                    domain=["Distribution", "ML Estimate"],
                    range=[0.7, 1.0],
                ),
                legend=None,
            ),
            tooltip=[
                "Allele Frequency",
                alt.Tooltip("Probability:Q", format=".6f"),
                "Type",
            ],
        )
        .properties(
            title="Allele Frequency Distribution (ML Estimate in Red)",
            width=500,
            height=300,
        )
        .configure_view(strokeWidth=0)
        .configure_axis(grid=False)
    )


def visualize_observations(record, sample_name):
    """Visualize observations from OBS field"""
    sample = record.samples[sample_name]
    obs = sample["OBS"]

    if len(obs) != 1:
        raise ValueError(
            f"Unexpected number of string values for OBS (expected: 1, found: {len(obs)})"
        )
    obs_string = obs[0]
    # TODO maybe the check for "." is not needed (pysam should recognize that as missing value).
    if obs_string is None or obs_string == ".":
        obs_string = ""

    ref_observations = []
    alt_observations = []

    # Pattern: COUNT + 2-letter odds + 8 single chars
    pattern = r"(\d+)([a-zA-Z]{2})(.)(.)(.)(.)(.)(.)(.)(.)"
    matches = re.findall(pattern, obs_string)

    strand_map = {"+": "Forward strand", "-": "Reverse strand", "*": "Both strands"}
    read_pos_map = {
        "^": "Most common position",
        "*": "Other position",
        ".": "Irrelevant position",
    }
    orientation_map = {
        ">": "F1R2 orientation",
        "<": "F2R1 orientation",
        "*": "Unknown orientation",
        "!": "Non-standard orientation",
    }
    softclip_map = {"$": "Soft clipped", ".": "No soft clipping"}
    indel_map = {"*": "Contains indel", ".": "No indel"}

    for idx, match in enumerate(matches):
        count = int(match[0])
        odds_code = match[1]

        # 8 fields from the match
        edit_distance_char = match[2]
        strand = match[5]
        orientation = match[6]
        read_position = match[7]
        softclip = match[8]
        indel = match[9]

        allele_type = odds_code[0].upper()
        kass = odds_code[1]

        kr_names = {
            "N": "None",
            "E": "Equal",
            "B": "Barely",
            "P": "Positive",
            "S": "Strong",
            "V": "Very Strong",
            "n": "None",
            "e": "Equal",
            "b": "Barely",
            "p": "Positive",
            "s": "Strong",
            "v": "Very Strong",
        }

        # Convert edit distance
        edit_distance = (
            0
            if edit_distance_char == "."
            else (int(edit_distance_char) if edit_distance_char.isdigit() else 0)
        )

        obs_entry = {
            "obs_index": idx,
            "count": count,
            "Posterior Odds": kr_names.get(kass, kass),
            "Strand": strand_map.get(strand, strand),
            "Read Position": read_pos_map.get(read_position, read_position),
            "Orientation": orientation_map.get(orientation, orientation),
            "Softclip": softclip_map.get(softclip, softclip),
            "Indel": indel_map.get(indel, indel),
            "Edit Distance": edit_distance,
        }

        if allele_type == "A":
            alt_observations.append(obs_entry)
        else:
            ref_observations.append(obs_entry)

    metrics = [
        "Posterior Odds",
        "Strand",
        "Read Position",
        "Orientation",
        "Softclip",
        "Indel",
        "Edit Distance",
    ]

    odds_order = ["None", "Equal", "Barely", "Positive", "Strong", "Very Strong"]
    odds_colors = ["#AAAAAA", "#999999", "#D4EFF7", "#AFDFEE", "#6CC5E0", "#2DACD2"]

    max_count = max(
        sum(o["count"] for o in ref_observations) or 0,
        sum(o["count"] for o in alt_observations) or 0,
    )

    def create_panel(observations, allele, show_y_axis=True, show_legend=True):
        if not observations:
            return (
                alt.Chart(pd.DataFrame({"Metric": [], "Count": []}))
                .mark_bar()
                .properties(
                    width=220,
                    height=400,
                    title=f"{allele} Allele Observations (No Data)",
                )
            )

        rows = []
        for obs in observations:
            for metric in metrics:
                rows.append(
                    {
                        "Metric": metric,
                        "Category": str(obs[metric]),
                        "Count": obs["count"],
                        "obs_index": obs["obs_index"],
                    }
                )

        df = pd.DataFrame(rows)

        # Determine edit distance domain
        edit_values = (
            df[df["Metric"] == "Edit Distance"]["Category"].astype(int).unique()
        )
        edit_domain = None
        if len(edit_values) == 1:
            k = int(edit_values[0])
            edit_domain = [0, k] if k > 0 else [0, 1]

        edit_scale = (
            alt.Scale(scheme="reds", domain=edit_domain)
            if edit_domain
            else alt.Scale(scheme="reds")
        )

        base = alt.Chart(df).encode(
            alt.X("Metric:N", sort=metrics, title=None),
            alt.Y(
                "Count:Q",
                stack="zero",
                scale=alt.Scale(domain=[0, max_count]),
                title="Count" if show_y_axis else None,
                axis=None if not show_y_axis else alt.Axis(),
            ),
            alt.Order("obs_index:Q"),
            alt.Tooltip(["Metric", "Category", "Count"]),
        )

        odds_layer = (
            base.transform_filter(alt.datum.Metric == "Posterior Odds")
            .mark_bar(size=18)
            .encode(
                alt.Color(
                    "Category:N",
                    scale=alt.Scale(domain=odds_order, range=odds_colors),
                    legend=alt.Legend(title="Posterior Odds") if show_legend else None,
                )
            )
        )

        edit_layer = (
            base.transform_filter(alt.datum.Metric == "Edit Distance")
            .mark_bar(size=18)
            .encode(
                alt.Color(
                    "Category:Q",
                    scale=edit_scale,
                    legend=alt.Legend(title="Edit distance") if show_legend else None,
                )
            )
        )

        other_layer = (
            base.transform_filter(
                (alt.datum.Metric != "Posterior Odds")
                & (alt.datum.Metric != "Edit Distance")
            )
            .mark_bar(size=18)
            .encode(
                alt.Color(
                    "Category:N",
                    legend=alt.Legend(title="Category") if show_legend else None,
                )
            )
        )

        return (odds_layer + edit_layer + other_layer).properties(
            width=220, height=400, title=f"{allele} Allele Observations"
        )

    ref_chart = create_panel(ref_observations, "REF", True, False)
    alt_chart = create_panel(alt_observations, "ALT", True, True)

    return (
        alt.hconcat(ref_chart, alt_chart, spacing=10)
        .resolve_scale(y="shared")
        .configure_legend(orient="right")
        .configure_view(strokeWidth=0)
        .configure_axis(grid=False)
    )
