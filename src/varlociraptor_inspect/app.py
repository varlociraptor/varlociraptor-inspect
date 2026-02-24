import sys
print(sys.path)

import streamlit as st
from varlociraptor_inspect import plotting  # noqa

st.set_page_config(
    page_title="Varlociraptor Inspect",
)

st.title("Varlociraptor Inspect")
st.text("Visual inspection of Varlociraptor VCF records.")

# load record, either from text_input, or from URL params
record = st.text_input("Paste your Varlociraptor VCF record here")
# parse (you will have to add a header) and plot
