import streamlit as st
import pysam
import tempfile
import os
from varlociraptor_inspect import plotting


def main_view():
    st.set_page_config(
        page_title="Varlociraptor Inspect",
    )
    st.title("Varlociraptor Inspect")
    st.text("Visual inspection of Varlociraptor VCF records.")

    # Load record from text input
    record_text = st.text_area(
        "Paste your Varlociraptor VCF record here (including header lines starting with #)",
        height=200,
    )

    if record_text:
        try:
            # Check if header is included
            if not record_text.startswith("##fileformat"):
                # No header - add a dummy header
                # Extract contig from first data line
                lines = record_text.strip().split("\n")
                contig = "chr1"
                for line in lines:
                    if not line.startswith("#"):
                        contig = line.split("\t")[0]
                        break

                # Add minimal header
                header = f"""##fileformat=VCFv4.2
##contig=<ID={contig}>
"""
                record_text = header + record_text

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
                # Clean up temp file
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)

        except Exception as e:
            st.error(f"Error parsing VCF record: {str(e)}")
