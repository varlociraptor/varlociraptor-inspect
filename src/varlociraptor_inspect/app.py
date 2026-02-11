import streamlit as st
from . import plotting  # noqa

st.title("Varlociraptor Inspect")
st.text("Visual inspection of Varlociraptor VCF records.")

# load record, either from text_input, or from URL params
record = st.text_input("Paste your Varlociraptor VCF record here")
# parse (you will have to add a header) and plot
