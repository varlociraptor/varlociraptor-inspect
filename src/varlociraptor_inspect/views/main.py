import os
import re
import tempfile
from collections.abc import Sequence
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

        origin = str(js.location.origin)
        pathname = str(js.location.pathname)
        url = f"{origin}{pathname}?{query_string}"
    except ImportError:
        host = st.context.headers.get("Host", "localhost:8501")
        scheme = "http" if host.startswith(("localhost", "127.0.0.1")) else "https"
        url = f"{scheme}://{host}/?{query_string}"
    st.caption("Shareable link to this record:")
    st.code(url, language=None, wrap_lines=True)


async def render_webllm_chat(description: str, key: str) -> None:
    """Render an in-browser chat using WebLLM with native Streamlit widgets."""
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
        description + "\n\nThere is exactly ONE variant described above. The 'Event "
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
        "how many reads support the variant, use the 'ALT-supporting reads' count, "
        "never the REF-supporting count or the allele frequency (AF is a fraction "
        "between 0 and 1, not a read count - never confuse the two). Double-check "
        "which sample and which number you are citing before answering.\n\n"
        "When asked about strand bias, orientation bias, read-position bias, or "
        "softclip bias: the PROB_ARTIFACT event captures exactly these biases "
        "combined. A low PROB_ARTIFACT probability means there is little to no "
        "evidence of bias - state this explicitly, and do not conclude there is a "
        "bias just because the ALT-supporting reads happen to share a strand or "
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
                        "content": "How many reads support the variant in sample X?",
                    },
                    {
                        "role": "assistant",
                        "content": (
                            "4 reads support the ALT (variant) allele in sample X. "
                            "This is the 'ALT-supporting reads' count copied directly "
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
        await render_webllm_chat(description, key="url")

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
                        elif not line.startswith("#") and line.strip():
                            data_line = line
                            break

                    if data_line:
                        fields = data_line.split("\t")

                        if len(fields) < 8:
                            raise ValueError(
                                "VCF record must have at least 8 tab-separated columns"
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
                            num_samples = len(fields) - 9
                            sample_names = [
                                f"sample{i + 1}" for i in range(num_samples)
                            ]
                            column_header = (
                                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                                + "\t".join(sample_names)
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
                        record = next(vcf)
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
                                afd = afd[0] if afd else None
                            if afd is not None:
                                afd_fields[str(sname)] = str(afd)

                            obs = s.get("OBS")
                            if isinstance(obs, Sequence) and not isinstance(obs, str):
                                obs = obs[0] if obs else None
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
                        await render_webllm_chat(description, key="paste")
                finally:
                    if os.path.exists(tmp_path):
                        os.unlink(tmp_path)

            except Exception as e:  # noqa: BLE001
                st.error(f"Error parsing VCF record: {e!s}")
