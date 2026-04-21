import altair as alt
import pandas as pd
import re
from dataclasses import dataclass
from typing import Sequence
from typing import Self


def phred_to_prob(phred_value: float) -> float:
    """Convert PHRED score to probability"""
    return float(10 ** (-phred_value / 10))


@dataclass
class ProbEntry:
    event: str
    probability: float


@dataclass
class ProbData:
    entries: list[ProbEntry]

    @classmethod
    def from_record(cls, record) -> "ProbData":
        """Parse PROB_* fields from a pysam VariantRecord INFO column."""
        entries = []
        for key, value in record.info.items():
            if not key.startswith("PROB_"):
                continue
            event_name = key.removeprefix("PROB_")
            if isinstance(value, Sequence) and not isinstance(value, str):
                if len(value) == 0:
                    raise ValueError(f"Empty sequence for PROB field '{event_name}'")
                value = value[0]
            if isinstance(value, str):
                try:
                    value = float(value)
                except ValueError:
                    raise ValueError(
                        f"Cannot convert PROB value '{value}' to float for event '{event_name}'"
                    )
            probability = 0.0 if value == float("inf") else float(phred_to_prob(value))
            entries.append(ProbEntry(event=event_name, probability=probability))
        return cls(entries=entries)

    @classmethod
    def from_dict(cls, prob_fields: dict[str, str]) -> "ProbData":
        """Build ProbData directly from URL parameter dict {PROB_EVENT: phred_str}."""
        entries = []
        for key, value in prob_fields.items():
            event_name = key.removeprefix("PROB_")
            try:
                phred = float(value)
            except (ValueError, TypeError):
                continue
            probability = 0.0 if phred == float("inf") else phred_to_prob(phred)
            entries.append(ProbEntry(event=event_name, probability=probability))
        return cls(entries=entries)


@dataclass
class AFDEntry:
    allele_frequency: float
    probability: float
    entry_type: str = "Distribution"


@dataclass
class AFDData:
    sample_name: str
    entries: list[AFDEntry]

    @classmethod
    def from_record(cls, record, sample_name: str) -> Self | None:
        """Parse AFD field from a pysam VariantRecord sample."""
        sample = record.samples[sample_name]
        afd_entries = sample.get("AFD")
        if afd_entries is None:
            afd_entries = []
        elif isinstance(afd_entries, str):
            afd_entries = [afd_entries]
        elif not isinstance(afd_entries, Sequence):
            afd_entries = [afd_entries]
        return cls._parse_afd_entries(sample_name, afd_entries)

    @classmethod
    def from_string(cls, sample_name: str, afd_string: str) -> Self | None:
        """Build AFDData directly from a raw AFD string (e.g. from URL params)."""
        if not afd_string or afd_string == ".":
            return None
        return cls._parse_afd_entries(sample_name, [afd_string])

    @classmethod
    def _parse_afd_entries(
        cls, sample_name: str, afd_entries: "Sequence[object]"
    ) -> Self | None:
        entries = []
        for entry in afd_entries:
            if not isinstance(entry, str):
                raise ValueError(
                    f"AFD entry must be a string, got {type(entry).__name__}: {entry!r}"
                )
            for part in entry.split(","):
                if "=" not in part:
                    raise ValueError(
                        f"Invalid AFD format in '{part}': expected 'freq=phred'"
                    )
                try:
                    freq, phred = part.split("=")
                    freq = float(freq)
                    prob = phred_to_prob(float(phred))
                except (ValueError, TypeError) as e:
                    raise ValueError(f"Failed to parse AFD entry '{part}': {e}")
        if not entries:
            return None
        ml = max(entries, key=lambda e: e.probability)
        entries.append(
            AFDEntry(
                allele_frequency=ml.allele_frequency,
                probability=ml.probability,
                entry_type="ML Estimate",
            )
        )
        return cls(sample_name=sample_name, entries=entries)


@dataclass
class ObsEntry:
    obs_index: int
    count: int
    posterior_odds: str
    mapq: str
    strand: str
    read_position: str
    orientation: str
    softclip: str
    indel: str
    edit_distance: int
    allele_type: str


@dataclass
class OBSData:
    sample_name: str
    ref_observations: list[ObsEntry]
    alt_observations: list[ObsEntry]

    @classmethod
    def from_record(cls, record, sample_name: str) -> "OBSData":
        """Parse OBS field from a pysam VariantRecord sample."""
        sample = record.samples[sample_name]
        obs = sample.get("OBS")
        if obs is None:
            obs_string = ""
        elif isinstance(obs, str):
            obs_string = obs
        elif isinstance(obs, Sequence):
            obs_string = obs[0] if obs else ""
        else:
            obs_string = str(obs)
        if obs_string == ".":
            obs_string = ""
        return cls._parse_obs_string(sample_name, obs_string)

    @classmethod
    def from_string(cls, sample_name: str, obs_string: str) -> "OBSData":
        """Build OBSData directly from a raw OBS string (e.g. from URL params)."""
        if not obs_string or obs_string == ".":
            obs_string = ""
        return cls._parse_obs_string(sample_name, obs_string)

    @classmethod
    def _parse_obs_string(cls, sample_name: str, obs_string: str) -> "OBSData":
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
        ref_observations = []
        alt_observations = []
        pattern = r"(\d+)([a-zA-Z]{2})(.)(.)(.)(.)(.)(.)(.)(.)"
        for idx, match in enumerate(re.findall(pattern, obs_string)):
            count = int(match[0])
            odds_code = match[1]
            edit_distance_char = match[2]
            strand = match[5]
            orientation = match[6]
            read_position = match[7]
            softclip = match[8]
            indel = match[9]
            allele_type = odds_code[0].upper()
            kass = odds_code[1]
            mapq = "High MAPQ" if kass.isupper() else "Low MAPQ"
            edit_distance = (
                0
                if edit_distance_char == "."
                else (int(edit_distance_char) if edit_distance_char.isdigit() else 0)
            )
            entry = ObsEntry(
                obs_index=idx,
                count=count,
                posterior_odds=kr_names.get(kass, kass.upper()),
                mapq=mapq,
                strand=strand_map.get(strand, strand),
                read_position=read_pos_map.get(read_position, read_position),
                orientation=orientation_map.get(orientation, orientation),
                softclip=softclip_map.get(softclip, softclip),
                indel=indel_map.get(indel, indel),
                edit_distance=edit_distance,
                allele_type=allele_type,
            )
            if allele_type == "A":
                alt_observations.append(entry)
            else:
                ref_observations.append(entry)
        return cls(
            sample_name=sample_name,
            ref_observations=ref_observations,
            alt_observations=alt_observations,
        )


def visualize_event_probabilities(prob_data: ProbData):
    """Visualize event probabilities."""
    df = pd.DataFrame(
        [{"Event": e.event, "Probability": e.probability} for e in prob_data.entries]
    )
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            alt.X("Event:N"),
            alt.Y("Probability:Q"),
            tooltip=["Event", alt.Tooltip("Probability:Q", format=".6f")],
        )
        .properties(title="Event Probabilities", width=400, height=300)
    )


def visualize_allele_frequency_distribution(afd_data: AFDData):
    """Visualize allele frequency distribution."""
    df = pd.DataFrame(
        [
            {
                "Allele Frequency": e.allele_frequency,
                "Probability": e.probability,
                "Type": e.entry_type,
            }
            for e in afd_data.entries
        ]
    )
    return (
        alt.Chart(df)
        .mark_circle()
        .encode(
            alt.X("Allele Frequency:Q"),
            alt.Y("Probability:Q", axis=None),
            alt.Color(
                "Type:N",
                scale=alt.Scale(
                    domain=["Distribution", "ML Estimate"], range=["blue", "red"]
                ),
            ),
            alt.Size(
                "Type:N",
                scale=alt.Scale(
                    domain=["Distribution", "ML Estimate"], range=[60, 100]
                ),
                legend=None,
            ),
            alt.Opacity(
                "Type:N",
                scale=alt.Scale(
                    domain=["Distribution", "ML Estimate"], range=[0.7, 1.0]
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


def visualize_observations(obs_data: OBSData):
    """Visualize observations from OBS field."""
    metrics = [
        "Posterior Odds",
        "MAPQ",
        "Strand",
        "Read Position",
        "Orientation",
        "Softclip",
        "Indel",
        "Edit Distance",
    ]
    odds_order = ["None", "Equal", "Barely", "Positive", "Strong", "Very Strong"]
    odds_colors = ["#AAAAAA", "#999999", "#D4EFF7", "#AFDFEE", "#6CC5E0", "#2DACD2"]
    mapq_colors = {"High MAPQ": "#4575b4", "Low MAPQ": "#d73027"}
    strand_colors = {
        "Forward strand": "#1f77b4",
        "Reverse strand": "#ff7f0e",
        "Both strands": "#2ca02c",
    }
    read_pos_colors = {
        "Most common position": "#d62728",
        "Other position": "#9467bd",
        "Irrelevant position": "#8c564b",
    }
    orientation_colors = {
        "F1R2 orientation": "#e377c2",
        "F2R1 orientation": "#7f7f7f",
        "Unknown orientation": "#bcbd22",
        "Non-standard orientation": "#17becf",
    }
    softclip_colors = {"Soft clipped": "#ff9896", "No soft clipping": "#c5b0d5"}
    indel_colors = {"Contains indel": "#c49c94", "No indel": "#f7b6d2"}

    ref_observations = obs_data.ref_observations
    alt_observations = obs_data.alt_observations

    max_count = max(
        sum(o.count for o in ref_observations) or 0,
        sum(o.count for o in alt_observations) or 0,
    )

    def create_panel(
        observations: list[ObsEntry], allele: str, show_y_axis=True, show_legend=True
    ):
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
            attr_map = {
                "Posterior Odds": obs.posterior_odds,
                "MAPQ": obs.mapq,
                "Strand": obs.strand,
                "Read Position": obs.read_position,
                "Orientation": obs.orientation,
                "Softclip": obs.softclip,
                "Indel": obs.indel,
                "Edit Distance": obs.edit_distance,
            }
            for metric in metrics:
                rows.append(
                    {
                        "Metric": metric,
                        "Category": str(attr_map[metric]),
                        "Count": obs.count,
                        "obs_index": obs.obs_index,
                    }
                )
        df = pd.DataFrame(rows)
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
        mapq_layer = (
            base.transform_filter(alt.datum.Metric == "MAPQ")
            .mark_bar(size=18)
            .encode(
                alt.Color(
                    "Category:N",
                    scale=alt.Scale(
                        domain=list(mapq_colors.keys()),
                        range=list(mapq_colors.values()),
                    ),
                    legend=alt.Legend(title="MAPQ") if show_legend else None,
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
                & (alt.datum.Metric != "MAPQ")
                & (alt.datum.Metric != "Edit Distance")
            )
            .mark_bar(size=18)
            .encode(
                alt.Color(
                    "Category:N",
                    scale=alt.Scale(
                        domain=list(strand_colors)
                        + list(read_pos_colors)
                        + list(orientation_colors)
                        + list(softclip_colors)
                        + list(indel_colors),
                        range=list(strand_colors.values())
                        + list(read_pos_colors.values())
                        + list(orientation_colors.values())
                        + list(softclip_colors.values())
                        + list(indel_colors.values()),
                    ),
                    legend=alt.Legend(title="Category") if show_legend else None,
                )
            )
        )
        return (odds_layer + mapq_layer + edit_layer + other_layer).properties(
            width=220, height=400, title=f"{allele} Allele Observations"
        )

    has_ref = bool(ref_observations)
    has_alt = bool(alt_observations)
    if has_ref and has_alt:
        ref_chart = create_panel(ref_observations, "REF", True, False)
        alt_chart = create_panel(alt_observations, "ALT", True, True)
    elif has_ref:
        ref_chart = create_panel(ref_observations, "REF", True, True)
        alt_chart = create_panel(alt_observations, "ALT", True, False)
    elif has_alt:
        ref_chart = create_panel(ref_observations, "REF", True, False)
        alt_chart = create_panel(alt_observations, "ALT", True, True)
    else:
        ref_chart = create_panel(ref_observations, "REF", True, False)
        alt_chart = create_panel(alt_observations, "ALT", True, False)

    return (
        alt.hconcat(ref_chart, alt_chart, spacing=10)
        .resolve_scale(y="shared")
        .configure_legend(orient="right")
        .configure_view(strokeWidth=0)
        .configure_axis(grid=False)
    )
