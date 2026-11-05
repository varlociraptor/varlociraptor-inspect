import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from itertools import chain
import hashlib
import pysam
import streamlit as st

from varlociraptor_inspect import plotting
from varlociraptor_inspect.description import build_record_description
from varlociraptor_inspect.plotting import AFDData, OBSData, ProbData


def normalize_whitespace(text: str) -> str:
    """Normalize whitespace in VCF records - replace spaces with tabs in data lines."""
    lines = []
    for line in text.split("\n"):
        if line.startswith("#") or not line.strip():
            lines.append(line)
        else:
            lines.append(re.sub(r"[ \t]+", "\t", line.strip()))
    return "\n".join(lines)


def build_vcf_from_url_params() -> tuple[ProbData, list[AFDData], list[OBSData]] | None:
    """Build dataclass instances directly from URL query parameters."""
    params = st.query_params
    if not params:
        return None

    prob_fields: dict[str, str] = {}
    afd_fields: dict[str, str] = {}
    obs_fields: dict[str, str] = {}

    for key in params:
        if key.startswith("PROB_"):
            prob_fields[key] = params[key]
        elif key.startswith("AFD_"):
            afd_fields[key.removeprefix("AFD_")] = params[key]
        elif key.startswith("OBS_"):
            obs_fields[key.removeprefix("OBS_")] = params[key]

    if not prob_fields:
        return None

    prob_data = ProbData.from_dict(prob_fields)
    sample_names = sorted(set(chain(afd_fields.keys(), obs_fields.keys())))

    afd_data_list: list[AFDData] = []
    for sample in sample_names:
        afd = AFDData.from_string(sample, afd_fields.get(sample, ""))
        if afd is not None:
            afd_data_list.append(afd)

    obs_data_list: list[OBSData] = [
        OBSData.from_string(sample, obs_fields.get(sample, ""))
        for sample in sample_names
    ]

    return prob_data, afd_data_list, obs_data_list


def build_query_string(
    prob_fields: dict[str, str], afd_fields: dict[str, str], obs_fields: dict[str, str]
) -> str:
    """Build a URL query string from raw PROB/AFD/OBS field values."""
    from urllib.parse import urlencode

    params: list[tuple[str, str]] = []
    for event, phred in prob_fields.items():
        params.append((f"PROB_{event}", str(phred)))
    for sample, value in obs_fields.items():
        if value and value != ".":
            params.append((f"OBS_{sample}", str(value)))
    for sample, value in afd_fields.items():
        if value and value != ".":
            params.append((f"AFD_{sample}", str(value)))
    return urlencode(params)


def render_copy_link_button(query_string: str) -> None:
    """Show a shareable link to this record. st.code has a built-in copy-to-clipboard button."""
    try:
        import js  # type: ignore[import-not-found]

        # Under stlite, Python runs in a Web Worker with its own blob: location,
        # unrelated to the real page URL - use the known public deployment URL.
        _ = js  # confirms we're actually running under stlite
        url = f"https://varlociraptor.github.io/varlociraptor-inspect/?{query_string}"
    except ImportError:
        host = st.context.headers.get("Host", "localhost:8501")
        scheme = "http" if host.startswith(("localhost", "127.0.0.1")) else "https"
        url = f"{scheme}://{host}/?{query_string}"
    st.caption("Shareable link to this record:")
    st.code(url, language=None, wrap_lines=True)


def _answer_read_support(obs_by_sample: Mapping[str, OBSData | None]) -> str:
    lines = []
    for sample_name, obs in obs_by_sample.items():
        if obs is None:
            continue
        alt_total = sum(e.count for e in obs.alt_observations)
        lines.append(
            f"**{sample_name}:** {alt_total} fragments support the ALT (variant) allele."
        )
    return "\n\n".join(lines) if lines else "No read observation data available."


def _answer_bias(prob_data: ProbData) -> str:
    artifact = next((e for e in prob_data.entries if e.event == "ARTIFACT"), None)
    if artifact is None:
        return "No PROB_ARTIFACT value available for this record."
    if artifact.probability > 0.5:
        return (
            f"Likely artifact present - PROB_ARTIFACT probability is "
            f"{artifact.probability:.4f} (high)."
        )
    if artifact.probability > 0.05:
        return (
            f"Some evidence of an artifact - PROB_ARTIFACT probability is "
            f"{artifact.probability:.4f} (moderate)."
        )
    return (
        f"No strong evidence of an artifact - PROB_ARTIFACT probability is "
        f"{artifact.probability:.4f} (low)."
    )


def _answer_most_likely_event(prob_data: ProbData) -> str:
    if not prob_data.entries:
        return "No event probabilities available for this record."
    best = max(prob_data.entries, key=lambda e: e.probability)
    return (
        f"**{best.event}** is the most likely event, "
        f"with probability {best.probability:.4f}."
    )


def _answer_top_sample(obs_by_sample: Mapping[str, OBSData | None]) -> str:
    best_sample = None
    best_count = -1
    for sample_name, obs in obs_by_sample.items():
        if obs is None:
            continue
        alt_total = sum(e.count for e in obs.alt_observations)
        if alt_total > best_count:
            best_count = alt_total
            best_sample = sample_name
    if best_sample is None:
        return "No read observation data available."
    return f"**{best_sample}** has the most ALT-supporting fragments ({best_count})."


def _answer_total_observations(obs_by_sample: Mapping[str, OBSData | None]) -> str:
    lines = []
    for sample_name, obs in obs_by_sample.items():
        if obs is None:
            continue
        total = sum(e.count for e in obs.alt_observations) + sum(
            e.count for e in obs.ref_observations
        )
        lines.append(f"**{sample_name}:** {total} total observations.")
    return "\n\n".join(lines) if lines else "No read observation data available."


async def render_webllm_chat(
    description: str,
    key: str,
    prob_data: ProbData,
    afd_by_sample: Mapping[str, AFDData | None],
    obs_by_sample: Mapping[str, OBSData | None],
) -> None:
    """Render quick answers, an in-browser chat, and the underlying data, in tabs."""
    tab_quick, tab_chat, tab_data = st.tabs(
        ["⚡ Quick answers", "💬 Chat", "📄 Data used"]
    )

    with tab_quick:
        st.caption("Instant answers computed directly from the data - no LLM involved.")
        faq: list[tuple[str, str]] = [
            ("ALT-supporting fragments (per sample)?", "reads"),
            ("Any bias?", "bias"),
            ("Most likely event?", "event"),
            ("Which sample has most support?", "top_sample"),
            ("Total observations (per sample)?", "total_obs"),
        ]
        cols = st.columns(len(faq))
        for col, (label, action) in zip(cols, faq, strict=True):
            with col:
                if st.button(label, key=f"quick_{action}_{key}"):
                    st.session_state[f"quick_answer_{key}"] = {
                        "reads": lambda: _answer_read_support(obs_by_sample),
                        "bias": lambda: _answer_bias(prob_data),
                        "event": lambda: _answer_most_likely_event(prob_data),
                        "top_sample": lambda: _answer_top_sample(obs_by_sample),
                        "total_obs": lambda: _answer_total_observations(obs_by_sample),
                    }[action]()
        answer = st.session_state.get(f"quick_answer_{key}")
        if answer:
            st.info(answer)

    with tab_data:
        st.caption("This is the exact text the chat model is given as context.")
        st.text(description)

    with tab_chat:
        try:
            import js  # type: ignore[import-not-found]
            from pyodide.ffi import JsException, to_js  # type: ignore[import-not-found]
        except ImportError:
            st.info(
                "In-browser chat is only available in the "
                "[web deployment](https://varlociraptor.github.io/varlociraptor-inspect/)."
            )
            return

        engine_ready_key = f"webllm_ready_{key}"
        history_key = f"webllm_history_{key}"
        fingerprint_key = f"webllm_fingerprint_{key}"
        description_fingerprint = hashlib.md5(description.encode()).hexdigest()[:8]

        if engine_ready_key not in st.session_state:
            st.session_state[engine_ready_key] = False
        if st.session_state.get(fingerprint_key) != description_fingerprint:
            st.session_state[fingerprint_key] = description_fingerprint
            st.session_state[history_key] = []

        system_prompt = (
            description
            + "\n\nThere is exactly ONE variant described above. The 'Event "
            "probabilities' section lists different hypotheses about that SAME "
            "single variant (e.g. is it somatic, absent, an artifact, loss of "
            "heterozygosity) - these are not separate variants, and there is no "
            "'the LOH variant' or 'the artifact variant' as a distinct thing. When "
            "asked which hypothesis is most likely, compare the probabilities and "
            "name the one with the highest value, still referring to it as the "
            "same single variant under different hypotheses.\n\n"
            "You are a helpful assistant answering questions about the variant "
            "record above. Always cite the exact numbers from the data above, copied "
            "verbatim - never estimate, round, or recompute them yourself. When asked "
            "how many fragments support the variant, use the 'ALT-supporting fragments' count, "
            "never the REF-supporting count or the allele frequency (AF is a fraction "
            "between 0 and 1, not a read count - never confuse the two). Double-check "
            "which sample and which number you are citing before answering.\n\n"
            "When asked about strand bias, orientation bias, read-position bias, or "
            "softclip bias: the PROB_ARTIFACT event captures exactly these biases "
            "combined. A low PROB_ARTIFACT probability means there is little to no "
            "evidence of bias - state this explicitly, and do not conclude there is a "
            "bias just because the ALT-supporting fragments happen to share a strand or "
            "orientation category; a small number of reads naturally cluster by "
            "chance. Only call it a bias if PROB_ARTIFACT is high.\n\n"
            "Never answer in a single word or a short phrase - always explain your "
            "reasoning in 2-3 sentences, referencing the specific numbers that led to "
            "your answer. Do not give generic definitions - explain what these "
            "specific numbers imply about this specific variant. If something is not "
            "present in the data, say so instead of guessing. Never invent or "
            "calculate a ratio, percentage, or fraction that is not explicitly "
            "written in the data above - only copy numbers verbatim.\n\n"
            "Below are two examples of the correct answer format and reasoning "
            "style. These examples use invented placeholder numbers unrelated to "
            "the real record above - never reuse these specific numbers, only "
            "copy the reasoning pattern."
        )

        models = {
            "Llama 3.2 1B (recommended, ~0.9GB)": "Llama-3.2-1B-Instruct-q4f16_1-MLC",
            "Gemma 2 2B (~2GB)": "gemma-2-2b-it-q4f16_1-MLC",
            "Phi-3.5 mini (better quality, ~2.2GB)": "Phi-3.5-mini-instruct-q4f16_1-MLC",
            "Qwen2.5 0.5B (fastest, less accurate)": "Qwen2.5-0.5B-Instruct-q4f16_1-MLC",
        }

        st.caption(
            "Chat with this record using a small language model that runs entirely "
            "in your browser (nothing is sent to a server). Requires a WebGPU-capable "
            "browser (recent Chrome/Edge)."
        )

        if not st.session_state[engine_ready_key]:
            if not hasattr(js.navigator, "gpu") or js.navigator.gpu is None:
                st.warning(
                    "Your browser does not support WebGPU, which is required for "
                    "in-browser chat. Please use a recent version of Chrome or Edge."
                )
                return

            selected_model = st.selectbox(
                "Select model",
                options=list(models.keys()),
                key=f"webllm_model_{key}",
            )
            if st.button("Load model & start chat", key=f"webllm_load_{key}"):
                model_id = models[selected_model]
                with st.spinner(f"Downloading and loading {selected_model}..."):
                    try:
                        webllm = await js.eval(
                            "import('https://esm.run/@mlc-ai/web-llm@0.2.74')"
                        )
                        engine = await webllm.CreateMLCEngine(model_id)
                    except (JsException, OSError) as e:
                        st.error(f"Failed to load the model: {e!s}")
                        return
                    js.globalThis._webllmEngine = engine
                st.session_state[engine_ready_key] = True
                st.rerun()
            return

        for msg in st.session_state[history_key]:
            with st.chat_message(msg["role"]):
                st.write(msg["content"])

        if prompt := st.chat_input(
            "Ask a question about this record...",
            key=f"webllm_input_{key}",
        ):
            st.session_state[history_key].append({"role": "user", "content": prompt})
            with st.chat_message("user"):
                st.write(prompt)

            with st.chat_message("assistant"):
                with st.spinner("Thinking..."):
                    engine = js.globalThis._webllmEngine
                    few_shot_examples = [
                        {
                            "role": "user",
                            "content": "How many fragments support the variant in sample X?",
                        },
                        {
                            "role": "assistant",
                            "content": (
                                "4 fragments support the ALT (variant) allele in sample X. "
                                "This is the 'ALT-supporting fragments' count copied directly "
                                "from the data above, not calculated or estimated."
                            ),
                        },
                        {"role": "user", "content": "Is there a strand bias?"},
                        {
                            "role": "assistant",
                            "content": (
                                "No, there is no meaningful evidence of strand bias. "
                                "PROB_ARTIFACT is very low (close to 0). A low "
                                "PROB_ARTIFACT means little to no evidence of bias - a "
                                "high PROB_ARTIFACT would indicate bias instead."
                            ),
                        },
                        {
                            "role": "user",
                            "content": "Why is observation posterior odds very strong although the observation has an edit distance of 5 to the ALT allele?",
                        },
                        {
                            "role": "assistant",
                            "content": (
                                "Support can still be strong if the respective read has many low quality bases "
                                "that would explain the high edit distance while at the same time having even "
                                "worse likelihoods for all other considered alleles."
                            ),
                        },
                    ]
                    messages = [{"role": "system", "content": system_prompt}]
                    messages.extend(few_shot_examples)
                    messages.extend(st.session_state[history_key])
                    opts = to_js(
                        {"messages": messages, "temperature": 0},
                        dict_converter=js.Object.fromEntries,
                    )
                    try:
                        response = await engine.chat.completions.create(opts)
                        reply = str(response.choices[0].message.content)
                    except (JsException, OSError) as e:
                        reply = f"Error during inference: {e!s}"
                st.write(reply)
                st.session_state[history_key].append(
                    {"role": "assistant", "content": reply}
                )


async def main_view():
    st.set_page_config(page_title="Varlociraptor Inspect")
    st.title("Varlociraptor Inspect")
    st.text("Visual inspection of Varlociraptor VCF records.")
    with st.expander("ℹ️ How to use direct URL links"):
        st.markdown("""
        You can link directly to a visualization by encoding VCF data as URL parameters.

        **Parameters:**
        - `PROB_<event>` — Event probability in PHRED scale (e.g. `PROB_SOMATIC=2.5`)
        - `AFD_<sample>` — Allele frequency distribution (e.g. `AFD_tumor=0.0=0.01,0.5=10.5`)
        - `OBS_<sample>` — Observation string (e.g. `OBS_tumor=31Rv.p.+**..`)

        **Example:**
        http://localhost:8501/?PROB_SOMATIC=2.5&PROB_GERMLINE=4.6&AFD_tumor=0.0=0.01,0.1=10.5&OBS_tumor=31Rv.p.+**..
        This will immediately render the visualizations without needing to paste any VCF data manually.
                """)

    url_data = build_vcf_from_url_params()

    if url_data is not None:
        prob_data, afd_data_list, obs_data_list = url_data
        st.info("Loaded data from URL parameters")

        params = st.query_params
        prob_fields = {
            k.removeprefix("PROB_"): v
            for k, v in params.items()
            if k.startswith("PROB_")
        }
        afd_fields = {
            k.removeprefix("AFD_"): v for k, v in params.items() if k.startswith("AFD_")
        }
        obs_fields = {
            k.removeprefix("OBS_"): v for k, v in params.items() if k.startswith("OBS_")
        }
        record_key = hashlib.md5(str(sorted(params.items())).encode()).hexdigest()[:8]

        render_copy_link_button(build_query_string(prob_fields, afd_fields, obs_fields))

        st.header("Event Probabilities")
        event_chart = plotting.visualize_event_probabilities(prob_data)
        if event_chart is not None:
            st.altair_chart(event_chart, key=f"event_probs_{record_key}")
        else:
            st.warning("No event probability data available.")

        afd_by_sample: dict[str, AFDData | None] = {}
        obs_by_sample: dict[str, OBSData] = {}

        if not afd_data_list and not obs_data_list:
            st.warning(
                "No sample data found in URL parameters. Showing event probabilities only."
            )
        else:
            afd_by_sample = {d.sample_name: d for d in afd_data_list}
            obs_by_sample = {d.sample_name: d for d in obs_data_list}
            all_sample_names = sorted(
                set(chain(afd_by_sample.keys(), obs_by_sample.keys()))
            )

            for idx, sample_name in enumerate(all_sample_names, 1):
                st.divider()
                st.header(f"Sample {idx}: {sample_name}")

                afd = afd_by_sample.get(sample_name)
                st.subheader("Allele Frequency Distribution")
                if afd is not None:
                    st.altair_chart(
                        plotting.visualize_allele_frequency_distribution(afd),
                        key=f"afd_{record_key}_{sample_name}",
                    )
                else:
                    st.warning(
                        f"No allele frequency data available for sample {sample_name}."
                    )

                obs = obs_by_sample.get(sample_name)
                st.subheader("Observations")
                if obs is not None:
                    st.altair_chart(
                        plotting.visualize_observations(obs),
                        key=f"obs_{record_key}_{sample_name}",
                    )
                else:
                    st.warning(
                        f"No observation data available for sample {sample_name}."
                    )

        description = build_record_description(prob_data, afd_by_sample, obs_by_sample)
        st.divider()
        st.header("Chat with this record")
        await render_webllm_chat(
            description,
            key="url",
            prob_data=prob_data,
            afd_by_sample=afd_by_sample,
            obs_by_sample=obs_by_sample,
        )

    else:
        record_text = st.text_area(
            "Paste your Varlociraptor VCF record here (including header lines starting with #)",
            value="",
            height=200,
        )

        if record_text:
            try:
                record_text = normalize_whitespace(record_text)

                if not record_text.startswith("##fileformat"):
                    lines = record_text.strip().split("\n")
                    column_header = None
                    data_line = None

                    for line in lines:
                        if line.startswith("#CHROM"):
                            column_header = line
                            continue
                        if not line.strip():
                            continue
                        first_fields = line.split("\t")
                        if (
                            column_header is None
                            and len(first_fields) >= 8
                            and first_fields[0].strip().lstrip("#") == "CHROM"
                            and first_fields[1].strip() == "POS"
                        ):
                            column_header = line if line.startswith("#") else f"#{line}"
                            continue
                        if not line.startswith("#"):
                            data_line = line
                            break

                    if data_line:
                        fields = data_line.split("\t")

                        if len(fields) < 8:
                            raise ValueError(
                                "VCF record must have at least 8 tab-separated columns"
                            )

                        if len(fields) > 9:
                            format_field_count = len(fields[8].split(":"))
                            for i, sample_value in enumerate(fields[9:], start=1):
                                sample_field_count = len(sample_value.split(":"))
                                if sample_field_count != format_field_count:
                                    raise ValueError(
                                        f"Sample {i} has {sample_field_count} "
                                        f"colon-separated values, but the FORMAT "
                                        f"column declares {format_field_count} "
                                        f"fields ({fields[8]}). Check for a "
                                        f"missing/extra value, or a stray tab or "
                                        f"space inside a sample's OBS string."
                                    )

                        chrom = fields[0]
                        pos = int(fields[1])
                        info_field = fields[7] if len(fields) > 7 else ""
                        prob_fields = re.findall(r"PROB_(\w+)=", info_field)

                        header_lines = [
                            "##fileformat=VCFv4.2",
                            f"##contig=<ID={chrom},length={pos + 1000}>",
                        ]
                        for prob_field in prob_fields:
                            header_lines.append(
                                f"##INFO=<ID=PROB_{prob_field},Number=.,Type=Float>"
                            )
                        header_lines.extend(
                            [
                                "##FORMAT=<ID=DP,Number=1,Type=Integer>",
                                "##FORMAT=<ID=AF,Number=1,Type=Float>",
                                "##FORMAT=<ID=AFD,Number=.,Type=String>",
                                "##FORMAT=<ID=OBS,Number=1,Type=String>",
                                "##FORMAT=<ID=HINTS,Number=.,Type=String>",
                            ]
                        )

                        if not column_header:
                            if len(fields) >= 9:
                                num_samples = len(fields) - 9
                                sample_names = [
                                    f"sample{i + 1}" for i in range(num_samples)
                                ]
                                column_header = (
                                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                                    + "\t".join(sample_names)
                                )
                            else:
                                column_header = (
                                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                                )

                        record_text = (
                            "\n".join(header_lines)
                            + "\n"
                            + column_header
                            + "\n"
                            + data_line
                        )

                tmp_fd, tmp_path = tempfile.mkstemp(suffix=".vcf", text=True)
                try:
                    with os.fdopen(tmp_fd, "w") as tmp:
                        tmp.write(record_text)

                    with pysam.VariantFile(tmp_path) as vcf:
                        try:
                            record = next(vcf)
                        except StopIteration:
                            raise ValueError(
                                "No VCF record found in the pasted text."
                            ) from None
                        except Exception as parse_err:
                            raise ValueError(
                                f"Could not parse this VCF record - check that "
                                f"the INFO/FORMAT fields and sample values are "
                                f"well-formed ({parse_err!s})."
                            ) from parse_err
                        sample_names = list(record.samples.keys())

                        st.success(
                            f"Successfully parsed VCF record at {record.chrom}:{record.pos} "
                            f"with {len(sample_names)} sample(s)"
                        )

                        prob_fields = {}
                        for info_key, info_value in record.info.items():
                            if not info_key.startswith("PROB_"):
                                continue
                            if isinstance(info_value, Sequence) and not isinstance(
                                info_value, str
                            ):
                                info_value = info_value[0] if info_value else None
                            if info_value is not None:
                                prob_fields[info_key.removeprefix("PROB_")] = str(
                                    info_value
                                )

                        afd_fields = {}
                        obs_fields = {}
                        for sname in sample_names:
                            s = record.samples[sname]
                            afd = s.get("AFD")
                            if isinstance(afd, Sequence) and not isinstance(afd, str):
                                # pysam splits comma-separated Number=. fields into a
                                # tuple - AFD's full distribution has multiple points
                                # joined by commas, so rejoin rather than take [0].
                                afd = ",".join(str(v) for v in afd) if afd else None
                            if afd is not None:
                                afd_fields[str(sname)] = str(afd)

                            obs = s.get("OBS")
                            if isinstance(obs, Sequence) and not isinstance(obs, str):
                                obs = ",".join(str(v) for v in obs) if obs else None
                            if obs is not None:
                                obs_fields[str(sname)] = str(obs)

                        render_copy_link_button(
                            build_query_string(prob_fields, afd_fields, obs_fields)
                        )

                        record_key = hashlib.md5(
                            f"{record.chrom}:{record.pos}:{record.ref}:"
                            f"{','.join(str(a) for a in record.alts or ())}".encode()
                        ).hexdigest()[:8]

                        prob_data = ProbData.from_record(record)
                        st.header("Event Probabilities")
                        event_chart = plotting.visualize_event_probabilities(prob_data)
                        if event_chart is not None:
                            st.altair_chart(
                                event_chart,
                                key=f"event_probs_{record_key}",
                            )
                        else:
                            st.warning("No event probability data available.")

                        afd_by_sample: dict[str, AFDData | None] = {}
                        obs_by_sample: dict[str, OBSData] = {}

                        if not sample_names:
                            st.warning(
                                "No sample data found. Showing event probabilities only."
                            )
                        else:
                            for idx, sample_name in enumerate(sample_names, 1):
                                st.divider()
                                st.header(f"Sample {idx}: {sample_name}")

                                afd = AFDData.from_record(record, str(sample_name))
                                st.subheader("Allele Frequency Distribution")
                                if afd is not None:
                                    st.altair_chart(
                                        plotting.visualize_allele_frequency_distribution(
                                            afd
                                        ),
                                        key=f"afd_{record_key}_{sample_name}",
                                    )
                                else:
                                    st.warning(
                                        f"No allele frequency data available for sample {sample_name}."
                                    )

                                obs = OBSData.from_record(record, str(sample_name))
                                st.subheader("Observations")
                                st.altair_chart(
                                    plotting.visualize_observations(obs),
                                    key=f"obs_{record_key}_{sample_name}",
                                )

                            afd_by_sample = {
                                str(s): AFDData.from_record(record, str(s))
                                for s in sample_names
                            }
                            obs_by_sample = {
                                str(s): OBSData.from_record(record, str(s))
                                for s in sample_names
                            }

                        description = build_record_description(
                            prob_data,
                            afd_by_sample,
                            obs_by_sample,
                            variant_info={
                                "chrom": str(record.chrom),
                                "pos": str(record.pos),
                                "ref": str(record.ref),
                                "alt": ",".join(str(a) for a in record.alts or ()),
                            },
                        )
                        st.divider()
                        st.header("Chat with this record")
                        await render_webllm_chat(
                            description,
                            key="paste",
                            prob_data=prob_data,
                            afd_by_sample=afd_by_sample,
                            obs_by_sample=obs_by_sample,
                        )
                finally:
                    if os.path.exists(tmp_path):
                        os.unlink(tmp_path)

            except Exception as e:  # noqa: BLE001
                st.error(f"Error parsing VCF record: {e!s}")
