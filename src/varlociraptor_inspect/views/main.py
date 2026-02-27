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
    record = st.text_area(
        "Paste your Varlociraptor VCF record here (including header lines starting with #)",
        height=200
    )
    
    if record:
        try:
            # Check if we need to add header definitions
            needs_prob_defs = "INFO=<ID=PROB_" not in record
            needs_format_defs = "FORMAT=<ID=AF," not in record
            
            if needs_prob_defs or needs_format_defs:
                # Extract contig from the data line
                lines = record.strip().split('\n')
                contig = "chr1"
                for line in lines:
                    if not line.startswith('#'):
                        contig = line.split('\t')[0]
                        break
                
                # Build header additions
                header_lines = []
                
                if not record.startswith("##fileformat"):
                    header_lines.append("##fileformat=VCFv4.2")
                
                if "##contig=<ID=" not in record:
                    header_lines.append(f"##contig=<ID={contig}>")
                
                if needs_prob_defs:
                    header_lines.extend([
                        "##INFO=<ID=PROB_PRESENT,Number=1,Type=Float,Description=\"Probability present\">",
                        "##INFO=<ID=PROB_ABSENT,Number=1,Type=Float,Description=\"Probability absent\">",
                        "##INFO=<ID=PROB_ARTIFACT,Number=1,Type=Float,Description=\"Probability artifact\">"
                    ])
                
                if needs_format_defs:
                    header_lines.extend([
                        "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
                        "##FORMAT=<ID=AFD,Number=.,Type=String,Description=\"Allele frequency distribution\">",
                        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                        "##FORMAT=<ID=OBS,Number=1,Type=String,Description=\"Observations\">"
                    ])
                
                # Insert headers at the beginning
                if header_lines:
                    record = '\n'.join(header_lines) + '\n' + record
            
            # Write to temporary file
            tmp_fd, tmp_path = tempfile.mkstemp(suffix='.vcf', text=True)
            try:
                with os.fdopen(tmp_fd, 'w') as tmp:
                    tmp.write(record)
                
                # Parse VCF using pysam
                vcf = pysam.VariantFile(tmp_path)
                record = next(vcf)
                sample_names = list(record.samples.keys())
                
                st.success(f"Successfully parsed VCF record at {record.chrom}:{record.pos} with {len(sample_names)} sample(s)")
                
                # Display Event Probabilities
                st.subheader("Event Probabilities")
                chart1 = plotting.visualize_event_probabilities(record)
                st.altair_chart(chart1, use_container_width=True)
                
                # Display plots for each sample
                for idx, sample_name in enumerate(sample_names, 1):
                    st.markdown("---")
                    st.subheader(f"Sample {idx}: {sample_name}")
                    
                    st.markdown(f"**Allele Frequency Distribution - Sample {idx}**")
                    chart2 = plotting.visualize_allele_frequency_distribution(record, sample_name)
                    st.altair_chart(chart2, use_container_width=True)
                    
                    st.markdown(f"**Observations - Sample {idx}**")
                    chart3 = plotting.visualize_observations(record, sample_name)
                    st.altair_chart(chart3, use_container_width=True)
                
                vcf.close()
            finally:
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)
            
        except Exception as e:
            st.error(f"Error parsing VCF record: {str(e)}")