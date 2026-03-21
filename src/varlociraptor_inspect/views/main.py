import streamlit as st
import pysam
import tempfile
import os
import re
from varlociraptor_inspect import plotting


def build_vcf_from_url_params():
    """Build VCF record from URL query parameters"""
    params = st.query_params
    
    if not params:
        return None
    
    # Extract PROB_ fields
    prob_fields = {}
    afd_fields = {}
    obs_fields = {}
    
    for key in params:
        if key.startswith("PROB_"):
            prob_fields[key] = params[key]
        elif key.startswith("AFD_"):
            sample_name = key.replace("AFD_", "")
            afd_fields[sample_name] = params[key]
        elif key.startswith("OBS_"):
            sample_name = key.replace("OBS_", "")
            obs_fields[sample_name] = params[key]
    
    if not prob_fields:
        return None
    
    # Get all sample names
    sample_names = sorted(set(list(afd_fields.keys()) + list(obs_fields.keys())))
    
    if not sample_names:
        return None
    
    # Build INFO field
    info_parts = [f"{k}={v}" for k, v in prob_fields.items()]
    info_field = ";".join(info_parts)
    
    # Build sample columns
    format_field = "AF:AFD:DP:OBS"
    sample_columns = []
    
    for sample in sample_names:
        af = "0.5"  # Default
        afd = afd_fields.get(sample, "0.0=0.01")
        dp = "100"  # Default
        obs = obs_fields.get(sample, ".")
        
        sample_columns.append(f"{af}:{afd}:{dp}:{obs}")
    
    # Build complete data line
    data_line = f"chr1\t1000\t.\tA\tT\t.\t.\t{info_field}\t{format_field}\t" + "\t".join(sample_columns)
    
    return data_line


def main_view():
    st.set_page_config(
        page_title="Varlociraptor Inspect",
    )
    st.title("Varlociraptor Inspect")
    st.text("Visual inspection of Varlociraptor VCF records.")

    # Check for URL parameters and build VCF if present
    url_vcf = build_vcf_from_url_params()
    if url_vcf:
        st.info("Loaded data from URL parameters")

    # Load record from text input
    record_text = st.text_area(
        "Paste your Varlociraptor VCF record here (including header lines starting with #)",
        value=url_vcf or "",
        height=200,
    )

    if record_text:
        try:
            # Check if header is included
            if not record_text.startswith("##fileformat"):
                # No header - generate from the data line
                lines = record_text.strip().split("\n")

                # Find the column header line (#CHROM) and data line
                column_header = None
                data_line = None

                for line in lines:
                    if line.startswith("#CHROM"):
                        column_header = line
                    elif not line.startswith("#") and line.strip():
                        data_line = line
                        break

                if data_line:
                    fields = data_line.split("\t")
                    chrom = fields[0]
                    pos = int(fields[1])
                    info_field = fields[7] if len(fields) > 7 else ""

                    # Extract all PROB_ fields from INFO column
                    prob_fields = re.findall(r"PROB_(\w+)=", info_field)

                    # Generate VCF header
                    header_lines = [
                        "##fileformat=VCFv4.2",
                        f"##contig=<ID={chrom},length={pos + 1000}>",
                    ]

                    # Add PROB_ INFO fields (Number=. to allow multiple values)
                    for prob_field in prob_fields:
                        header_lines.append(
                            f"##INFO=<ID=PROB_{prob_field},Number=.,Type=Float>"
                        )

                    # Add standard Varlociraptor FORMAT fields
                    format_fields = [
                        "##FORMAT=<ID=DP,Number=1,Type=Integer>",
                        "##FORMAT=<ID=AF,Number=1,Type=Float>",
                        "##FORMAT=<ID=AFD,Number=.,Type=String>",
                        "##FORMAT=<ID=OBS,Number=1,Type=String>",
                        "##FORMAT=<ID=HINTS,Number=.,Type=String>",
                    ]
                    header_lines.extend(format_fields)

                    # Generate column header if not present
                    if not column_header:
                        if len(fields) < 8:
                            raise ValueError(
                                "VCF record must have at least 8 tab-separated columns"
                            )

                        columns = [
                            "#CHROM",
                            "POS",
                            "ID",
                            "REF",
                            "ALT",
                            "QUAL",
                            "FILTER",
                            "INFO",
                        ]

                        if len(fields) > 8:
                            num_samples = len(fields) - 9
                            sample_names = [f"sample{i + 1}" for i in range(num_samples)]
                            columns.extend(["FORMAT", *sample_names])

                        column_header = "\t".join(columns)

                    # Assemble complete VCF
                    record_text = (
                        "\n".join(header_lines)
                        + "\n"
                        + column_header
                        + "\n"
                        + data_line
                    )

            # Write to temp file
            tmp_fd, tmp_path = tempfile.mkstemp(suffix=".vcf", text=True)
            try:
                with os.fdopen(tmp_fd, "w") as tmp:
                    tmp.write(record_text)

                # Parse with context manager
                with pysam.VariantFile(tmp_path) as vcf:
                    record = next(vcf)
                    sample_names = list(record.samples.keys())

                    st.success(
                        f"Successfully parsed VCF record at {record.chrom}:{record.pos} with {len(sample_names)} sample(s)"
                    )

                    # Display Event Probabilities
                    st.header("Event Probabilities")
                    chart1 = plotting.visualize_event_probabilities(record)
                    st.altair_chart(chart1, use_container_width=True)

                    # Only show sample plots if samples exist
                    if len(sample_names) == 0:
                        st.warning(
                            "No sample columns found. Only Event Probabilities are shown."
                        )
                    else:
                        # Display plots for each sample
                        for idx, sample_name in enumerate(sample_names, 1):
                            st.divider()
                            st.header(f"Sample {idx}: {sample_name}")

                            st.subheader("Allele Frequency Distribution")
                            chart2 = plotting.visualize_allele_frequency_distribution(
                                record, sample_name
                            )
                            st.altair_chart(chart2, use_container_width=True)

                            st.subheader("Observations")
                            chart3 = plotting.visualize_observations(record, sample_name)
                            st.altair_chart(chart3, use_container_width=True)
            finally:
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)

        except Exception as e:
            st.error(f"Error parsing VCF record: {str(e)}")