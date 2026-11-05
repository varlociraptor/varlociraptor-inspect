import streamlit as st

try:
    from varlociraptor_inspect.views.main import main_view
except ImportError as error:
    st.error(
        "Some components could not be loaded. Please use the latest version of "
        "your web browser and reload the page."
    )
    st.caption(f"Details: {error}")
    st.stop()

await main_view()  # noqa: F704, PLE1142  # pyrefly: ignore[invalid-syntax] — stlite compiles with top-level await support
