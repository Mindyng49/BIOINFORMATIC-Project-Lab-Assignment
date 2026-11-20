from html import escape
from inspect import signature

import pandas as pd
import streamlit as st

import pairwise_alignment as pa
from pairwise_alignment import alignment_blocks, match_line

# Lab 1 is imported from the shared module extracted from lab1-Mindy.py.
try:
    from lab1_core import get_basic_analysis, retrieve_data
    LAB1_AVAILABLE = True
except Exception as exc:  # Biopython missing -> Lab 1 tab shows a notice
    LAB1_AVAILABLE = False
    LAB1_IMPORT_ERROR = str(exc)


# Page setup

st.set_page_config(page_title="Lab 2 - Pairwise Sequence Alignment", layout="wide")


STRETCH = (
    {"width": "stretch"}
    if "width" in signature(st.button).parameters
    else {"use_container_width": True}
)

LAB1 = "Lab 1 - Protein Retrieval"
LAB2 = "Lab 2 - Pairwise Alignment"

CSS = """
<style>
.aln-wrap {display:inline-block; max-width:100%; overflow:auto; max-height:620px; border:1px solid #d7dce5; border-radius:10px; background:#fff;}
table.aln {border-collapse:separate; border-spacing:0; font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size:12px;}
table.aln th {position:sticky; top:0; background:#1f2937; color:#fff; padding:2px 4px; min-width:32px; font-weight:600; z-index:2; line-height:1.15;}
table.aln th.idx {font-size:9px; font-weight:400; color:#cbd5e1;}
table.aln th.rowhead {position:sticky; left:0; top:auto; background:#1f2937; color:#fff; min-width:52px; text-align:right; padding-right:6px; z-index:3; line-height:1.15;}
table.aln th.corner {position:sticky; left:0; top:0; z-index:4; background:#111827;}
table.aln td {padding:1px 3px; text-align:center; min-width:32px; border-right:1px solid rgba(15,23,42,.07); border-bottom:1px solid rgba(15,23,42,.07); line-height:1.1;}
table.aln td .v {display:block; font-weight:600;}
table.aln td .d {display:block; font-size:9px; opacity:.8;}
.alnblock {font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size:13px; line-height:1.75; background:#ffffff; color:#0f172a; padding:10px 14px; border:1px solid #d7dce5; border-radius:10px; margin-bottom:8px; overflow-x:auto;}
.alnblock .row {white-space:pre; margin:0;}
.alnblock .lbl {display:inline-block; width:52px; margin-right:12px; color:#6b7280; user-select:none;}
.alnblock .tail {color:#94a3b8; margin-left:14px; font-size:11px;}
.alnblock .mt {color:#0f766e;}
.alnblock .mm {color:#b91c1c; background:#fee2e2; border-radius:3px;}
.alnblock .gp {color:#64748b; background:#eef2f7; border-radius:3px;}
.alnblock .mk {color:#0f172a;}
.legend span {display:inline-block; padding:2px 8px; border-radius:999px; font-size:11px; margin-right:6px;}
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)


# Session state

PERSISTED = {
    "seq_a": "MKTAYIAKQRQISFVKSHFSRQ",
    "seq_b": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
    "p_mode": "Global (Needleman-Wunsch)",
    "p_match": 2,
    "p_mismatch": -1,
    "p_gap": -2,
}

# Streamlit throws away the state of a widget that is not rendered on a run, so
# leaving Lab 2 for Lab 1 would silently reset the sequences and the scores.
# Therefore the truth lives in plain (non-widget) keys, and each widget reads
# its value from there. A "nonce" is part of the widget key so that a
# programmatic change (e.g. "Use as Sequence 1") recreates the widget with the
# new value instead of fighting the widget's own state.
for _key, _default in PERSISTED.items():
    st.session_state.setdefault(_key, _default)
    st.session_state.setdefault("n_" + _key, 0)

for _key, _value in {"page": LAB1, "lab1_records": {}, "last_result": None}.items():
    st.session_state.setdefault(_key, _value)


def widget_key(key: str) -> str:
    """Versioned widget key; bump the nonce to force a fresh widget."""
    return f"w_{key}_{st.session_state['n_' + key]}"


def pull(key: str) -> None:
    """on_change callback: copy the widget value into the authoritative store."""
    st.session_state[key] = st.session_state[widget_key(key)]


def remember(key, value) -> None:
    """Set a value programmatically and recreate the widget with it."""
    st.session_state[key] = value
    st.session_state["n_" + key] += 1


def send_to_lab2(slot: str, sequence: str) -> None:
    """Callback: push a Lab 1 sequence into Lab 2 and switch tab."""
    remember(slot, sequence)
    st.session_state["page"] = LAB2


def load_example(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> None:
    remember("seq_a", seq_a)
    remember("seq_b", seq_b)
    remember("p_match", match)
    remember("p_mismatch", mismatch)
    remember("p_gap", gap)


def swap_sequences() -> None:
    a, b = st.session_state["seq_a"], st.session_state["seq_b"]
    remember("seq_a", b)
    remember("seq_b", a)



# HTML helpers

def _ramp(value, vmin, vmax):
    """Light -> deep blue colour for a matrix value; returns (bg, fg)."""
    if vmax == vmin:
        t = 0.45
    else:
        t = (value - vmin) / (vmax - vmin)
    t = max(0.0, min(1.0, t))
    r = int(248 - t * (248 - 30))
    g = int(250 - t * (250 - 64))
    b = int(255 - t * (255 - 175))
    fg = "#ffffff" if t > 0.55 else "#0f172a"
    return f"rgb({r},{g},{b})", fg


def matrix_html(res, row_lo, row_hi, col_lo, col_hi, show_path=True):
    """Render the DP matrix as an HTML table with heat colours and the path."""
    visible = [
        res.matrix[i][j]
        for i in range(row_lo, row_hi + 1)
        for j in range(col_lo, col_hi + 1)
    ]
    vmin, vmax = (min(visible), max(visible)) if visible else (0, 0)

    role = {}
    if show_path:
        for k, cell in enumerate(res.path):
            if k == 0:
                role[cell] = "start"
            elif k == len(res.path) - 1:
                role[cell] = "end"
            else:
                role[cell] = "path"

    style_path = {
        "path": ("background:#f59e0b;color:#1f2937;font-weight:700;box-shadow:inset 0 0 0 2px #b45309;"),
        "start": ("background:#16a34a;color:#ffffff;font-weight:700;box-shadow:inset 0 0 0 2px #14532d;"),
        "end": ("background:#dc2626;color:#ffffff;font-weight:700;box-shadow:inset 0 0 0 2px #7f1d1d;"),
    }

    out = ['<div class="aln-wrap"><table class="aln"><thead><tr>']
    out.append('<th class="corner" rowspan="2">i \\ j</th>')
    for j in range(col_lo, col_hi + 1):
        letter = escape(res.seq2[j - 1]) if j > 0 else "&middot;"
        out.append(f'<th>{letter}</th>')
    out.append("</tr><tr>")
    for j in range(col_lo, col_hi + 1):
        out.append(f'<th class="idx">{j}</th>')
    out.append("</tr></thead><tbody>")

    for i in range(row_lo, row_hi + 1):
        letter = escape(res.seq1[i - 1]) if i > 0 else "&middot;"
        out.append(f'<tr><th class="rowhead">{letter} {i}</th>')
        for j in range(col_lo, col_hi + 1):
            value = res.matrix[i][j]
            arrows = "".join(pa.ARROW[d] for d in res.traceback[i][j])
            if not res.is_local and i == 0 and j == 0:
                arrows = ""
            if res.is_local and res.matrix[i][j] == 0:
                arrows = ""
            bg, fg = _ramp(value, vmin, vmax)
            extra = ""
            if (i, j) in role:
                extra = style_path[role[(i, j)]]
            out.append(
                f'<td style="background:{bg};color:{fg};{extra}">'
                f'<span class="v">{value:g}</span>'
                f'<span class="d">{arrows}</span></td>'
            )
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def _seq_line(label, seq, ops):
    spans = []
    for ch, op in zip(seq, ops):
        if ch == "-":
            cls = "gp"
        elif op == "match":
            cls = "mt"
        elif op == "mismatch":
            cls = "mm"
        else:  # gap opposite this residue
            cls = "mk"
        spans.append(f'<span class="{cls}">{escape(ch)}</span>')
    return f'<div class="row"><span class="lbl">{label}</span>{"".join(spans)}</div>'


def alignment_html(res, width=60):
    if not res.alignment_length:
        return None
    parts = []
    for b, (s1, mid, s2) in enumerate(alignment_blocks(res.aligned1, res.aligned2, width)):
        c0 = b * width
        ops = res.operations[c0:c0 + len(s1)]
        marks = []
        for ch in mid:
            cls = "mt" if ch == "|" else ("mm" if ch == "." else "mk")
            marks.append(f'<span class="{cls}">{escape(ch)}</span>')
        parts.append('<div class="alnblock">')
        parts.append(_seq_line("Seq 1", s1, ops))
        parts.append(f'<div class="row"><span class="lbl"></span>{"".join(marks)}'
                     f'<span class="tail">cols {c0 + 1}-{c0 + len(s1)}</span></div>')
        parts.append(_seq_line("Seq 2", s2, ops))
        parts.append("</div>")
    return "".join(parts)


def traceback_rows(res):
    rows = []
    for k in range(len(res.path) - 1, 0, -1):
        i, j = res.path[k]
        pi, pj = res.path[k - 1]
        op = res.operations[k - 1]
        if (pi, pj) == (i - 1, j - 1):
            move, arrow = "diagonal", "\u2196"
        elif (pi, pj) == (i - 1, j):
            move, arrow = "up", "\u2191"
        else:
            move, arrow = "left", "\u2190"

        if op == "match":
            detail = f"{res.seq1[i - 1]} = {res.seq2[j - 1]}  (match)"
        elif op == "mismatch":
            detail = f"{res.seq1[i - 1]} \u2260 {res.seq2[j - 1]}  (mismatch)"
        elif op == "gap_seq2":
            detail = f"{res.seq1[i - 1]} vs -  (gap in seq 2)"
        else:
            detail = f"- vs {res.seq2[j - 1]}  (gap in seq 1)"

        rows.append(
            {
                "Step": len(res.path) - k,
                "Cell (i, j)": f"({i}, {j})",
                "Score": res.matrix[i][j],
                "Move": f"{arrow} {move}",
                "From": f"({pi}, {pj})",
                "Aligned column": detail,
            }
        )
    return rows


def matrix_csv(res) -> str:
    header = ["i \\ j"] + ["-"] + list(res.seq2)
    lines = [",".join(header)]
    for i in range(len(res.matrix)):
        label = res.seq1[i - 1] if i > 0 else "-"
        lines.append(",".join([f"{label}_{i}"] + [str(v) for v in res.matrix[i]]))
    return "\n".join(lines)



# Sidebar navigation

st.sidebar.title("Bioinformatics Lab")
st.sidebar.radio("Navigate", [LAB1, LAB2], key="page")
st.sidebar.divider()
if st.session_state["lab1_records"]:
    st.sidebar.caption("Sequences retrieved this session: "
                       + ", ".join(st.session_state["lab1_records"]))


# LAB 1

if st.session_state["page"] == LAB1:
    st.title("Lab 1 - Mindy Ng")

    if not LAB1_AVAILABLE:
        st.error(f"Biopython is required for Lab 1: {LAB1_IMPORT_ERROR}")
        st.stop()

    with st.form("protein_search"):
        protein_id = st.text_input("Enter Uniprot ID", placeholder="e.g. P04637")
        retrieve = st.form_submit_button("Retrieve")

    if retrieve:
        if protein_id.strip() == "":
            st.warning("Please enter a UniProt ID.")
        else:
            try:
                with st.spinner("Retrieving protein data..."):
                    record = retrieve_data(protein_id.strip())
                    sequence = str(record.seq)
                    length, composition, weight, pi = get_basic_analysis(sequence)

                st.session_state["lab1_records"][protein_id.strip().upper()] = {
                    "name": record.name if record.name != "<unknown name>" else protein_id.strip().upper(),
                    "description": record.description,
                    "sequence": sequence,
                }
                st.success("Protein data retrieved successfully!")

                left, right = st.columns([1.2, 1])

                with left:
                    st.subheader("\U0001F4C4 Retrieved Protein")
                    st.write("**Name:**", record.name)
                    st.write("**Description:**", record.description)
                    st.write("**Sequence:**")
                    formatted_sequence = "\n".join(
                        sequence[i:i + 60] for i in range(0, len(sequence), 60)
                    )
                    st.code(formatted_sequence)

                    c1, c2, _ = st.columns([1, 1, 1.4])
                    c1.button(
                        "\u2192 Use as Sequence 1",
                        key="lab1_to_a",
                        on_click=send_to_lab2,
                        args=("seq_a", sequence),
                        **STRETCH,
                    )
                    c2.button(
                        "\u2192 Use as Sequence 2",
                        key="lab1_to_b",
                        on_click=send_to_lab2,
                        args=("seq_b", sequence),
                        **STRETCH,
                    )

                with right:
                    st.subheader("Basic Protein Analysis")
                    st.write("Sequence Length:", length)
                    st.write("Amino Acid Composition:")
                    composition_text = ", ".join(
                        [f"{amino_acid}: {frequency}"
                         for amino_acid, frequency in composition.items()]
                    )
                    st.write(composition_text)
                    st.write("Molecular Weight:", weight)
                    st.write("Isoelectric Point:", pi)

            except Exception as e:
                st.error(f"Unable to retrieve protein data: {e}")

    # ---- session cache ----------------------------------------------------
    if st.session_state["lab1_records"]:
        with st.expander(f"Retrieved sequences ({len(st.session_state['lab1_records'])})"):
            for acc, rec in st.session_state["lab1_records"].items():
                st.write(f"**{acc}** - {rec['name']} - {len(rec['sequence'])} aa")
                cc1, cc2, _ = st.columns([1, 1, 3])
                cc1.button("\u2192 Seq 1", key=f"cache_a_{acc}",
                           on_click=send_to_lab2, args=("seq_a", rec["sequence"]))
                cc2.button("\u2192 Seq 2", key=f"cache_b_{acc}",
                           on_click=send_to_lab2, args=("seq_b", rec["sequence"]))

    st.info(
        "Retrieve two proteins here, then send each one to Lab 2 to align them. "
        "Use the sidebar to open **Lab 2 - Pairwise Alignment**."
    )


# LAB 2

else:
    st.title("Lab 2 - Pairwise Sequence Alignment")
    st.caption(
        "Global = Needleman-Wunsch \u00b7 Local = Smith-Waterman \u00b7 "
    )

    # ---------------- examples / loaders ---------------------------------
    lc1, lc2, lc3, lc4, lc5 = st.columns([1, 1, 1, 1, 2])
    lc1.button(
        "Protein example",
        on_click=load_example,
        args=("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", 2, -1, -2),
        **STRETCH,
    )
    lc2.button(
        "DNA example",
        on_click=load_example,
        args=("GATTACA", "GCATGCU", 1, -1, -1),
        **STRETCH,
    )
    lc3.button("Swap 1 \u2194 2", on_click=swap_sequences, **STRETCH)
    lc4.button(
        "Clear",
        on_click=load_example,
        args=("", "", 2, -1, -2),
        **STRETCH,
    )
    with lc5:
        if st.session_state["lab1_records"]:
            pick = st.selectbox(
                "From Lab 1",
                list(st.session_state["lab1_records"].keys()),
                label_visibility="collapsed",
            )
            a_col, b_col = st.columns(2)
            a_col.button(
                "Use as Seq 1",
                key=f"l2_a_{pick}",
                on_click=send_to_lab2,
                args=("seq_a", st.session_state["lab1_records"][pick]["sequence"]),
                **STRETCH,
            )
            b_col.button(
                "Use as Seq 2",
                key=f"l2_b_{pick}",
                on_click=send_to_lab2,
                args=("seq_b", st.session_state["lab1_records"][pick]["sequence"]),
                **STRETCH,
            )
        else:
            st.caption("Retrieve a protein in Lab 1 to load it here.")

    # ---------------- inputs ---------------------------------------------
    st.subheader("Input sequences")
    in1, in2 = st.columns(2)
    with in1:
        st.text_area(
            "Sequence 1 - matrix rows (i)",
            value=st.session_state["seq_a"],
            key=widget_key("seq_a"),
            on_change=pull,
            args=("seq_a",),
            height=150,
            placeholder="Paste a sequence or FASTA record...",
        )
        st.caption(f"{len(pa.clean_sequence(st.session_state['seq_a']))} residues")
    with in2:
        st.text_area(
            "Sequence 2 - matrix columns (j)",
            value=st.session_state["seq_b"],
            key=widget_key("seq_b"),
            on_change=pull,
            args=("seq_b",),
            height=150,
            placeholder="Paste a sequence or FASTA record...",
        )
        st.caption(f"{len(pa.clean_sequence(st.session_state['seq_b']))} residues")

    # ---------------- scoring --------------------------------------------
    st.subheader("Scoring & alignment type")
    MODES = ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"]
    s1, s2, s3, s4 = st.columns(4)
    with s1:
        st.radio(
            "Alignment",
            MODES,
            index=MODES.index(st.session_state["p_mode"]),
            key=widget_key("p_mode"),
            on_change=pull,
            args=("p_mode",),
            horizontal=False,
        )
    with s2:
        st.number_input("Match score", value=st.session_state["p_match"],
                        key=widget_key("p_match"), on_change=pull,
                        args=("p_match",), step=1, format="%d")
    with s3:
        st.number_input("Mismatch score", value=st.session_state["p_mismatch"],
                        key=widget_key("p_mismatch"), on_change=pull,
                        args=("p_mismatch",), step=1, format="%d")
    with s4:
        st.number_input("Gap score", value=st.session_state["p_gap"],
                        key=widget_key("p_gap"), on_change=pull,
                        args=("p_gap",), step=1, format="%d")

    if st.session_state["p_gap"] > 0:
        st.warning("A positive gap score is unusual - gaps would be rewarded, "
                   "so the alignment can be dominated by gaps.")

    run = st.button("\u25B6 Run alignment", type="primary")

    if run:
        mode = "local" if st.session_state["p_mode"].startswith("Local") else "global"
        try:
            result = pa.align(
                st.session_state["seq_a"],
                st.session_state["seq_b"],
                match=st.session_state["p_match"],
                mismatch=st.session_state["p_mismatch"],
                gap=st.session_state["p_gap"],
                mode=mode,
            )
            st.session_state["last_result"] = result
        except ValueError as exc:
            st.session_state["last_result"] = None
            st.error(str(exc))

    res = st.session_state["last_result"]

    # ---------------- results --------------------------------------------
    if res is None:
        st.info("Enter two sequences, choose the scoring, then press **Run alignment**.")
    else:
        st.divider()
        kind = "Local (Smith-Waterman)" if res.is_local else "Global (Needleman-Wunsch)"
        st.subheader(f"Result \u00b7 {kind}")

        current = (
            pa.clean_sequence(st.session_state["seq_a"]),
            pa.clean_sequence(st.session_state["seq_b"]),
            float(st.session_state["p_match"]),
            float(st.session_state["p_mismatch"]),
            float(st.session_state["p_gap"]),
            "local" if st.session_state["p_mode"].startswith("Local") else "global",
        )
        used = (
            res.seq1,
            res.seq2,
            float(res.match_score),
            float(res.mismatch_score),
            float(res.gap_score),
            res.mode,
        )
        if current != used:
            st.warning(
                "The inputs or scoring changed after this alignment was computed - "
                "press **Run alignment** to refresh the result."
            )

        st.caption(
            f"Scoring used: match {res.match_score:+g}, mismatch {res.mismatch_score:+g}, "
            f"gap {res.gap_score:+g} \u00b7 seq 1 = {res.n_rows} residues, seq 2 = {res.n_cols} residues "
            f"({(res.n_rows + 1) * (res.n_cols + 1):,} matrix cells)"
        )

        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Alignment score", f"{res.score:g}")
        m2.metric("Identity", f"{res.identity:.1f}%")
        m3.metric("Matches / mismatches", f"{res.n_matches} / {res.n_mismatches}")
        m4.metric("Gaps", res.n_gaps)
        m5.metric("Alignment length", res.alignment_length)

        if res.is_local:
            i0, j0 = res.start_cell
            i1, j1 = res.end_cell
            st.caption(
                f"Local alignment spans seq 1 positions {i0 + 1}-{i1} and "
                f"seq 2 positions {j0 + 1}-{j1}; the maximum value {res.score:g} is at "
                f"cell ({i1}, {j1})."
            )

        st.markdown("#### Aligned sequences")
        if res.alignment_length == 0:
            st.warning(
                "No positive-scoring local alignment exists for these sequences and "
                "scores (best score = 0)."
            )
        else:
            html = alignment_html(res)
            st.markdown(html, unsafe_allow_html=True)
            st.caption(
                "Green = match, red = mismatch, grey = gap "
                "(\u2018-\u2019); the middle row marks identity with \u2018|\u2019."
            )

        # ---------------- matrix -----------------------------------------
        st.markdown("#### Dynamic-programming matrix")
        n, m = res.n_rows, res.n_cols
        total = (n + 1) * (m + 1)
        show_path = st.toggle("Highlight trace-back path", value=True)

        if total <= 6000:
            view = st.radio(
                "Matrix view", ["Full matrix", "Traceback window"],
                horizontal=True, index=0,
            )
        else:
            view = "Traceback window"
            st.info(
                f"The matrix has {total:,} cells, so only a window is drawn. "
                "Download the CSV for the complete matrix."
            )

        MAX_ROWS, MAX_COLS = 60, 60
        if view == "Full matrix":
            r0, r1, c0, c1 = 0, n, 0, m
        else:
            pi = [p[0] for p in res.path] or [0]
            pj = [p[1] for p in res.path] or [0]
            # centre a 40x40 window on the middle of the trace-back
            ci = (min(pi) + max(pi)) // 2
            cj = (min(pj) + max(pj)) // 2
            r0 = max(0, min(ci - 20, n - 39))
            r1 = min(n, r0 + 39)
            c0 = max(0, min(cj - 20, m - 39))
            c1 = min(m, c0 + 39)
            if n <= 150 and m <= 150:
                r0, r1 = st.slider("Rows (i) shown", 0, n, (r0, r1))
                c0, c1 = st.slider("Columns (j) shown", 0, m, (c0, c1))
            r1 = min(r1, r0 + MAX_ROWS)
            c1 = min(c1, c0 + MAX_COLS)
            if r1 < n or c1 < m:
                st.caption(
                    f"Window: rows {r0}-{r1} of {n}, columns {c0}-{c1} of {m}. "
                    "Use the CSV download to see every cell."
                )

        st.markdown(
            matrix_html(res, r0, r1, c0, c1, show_path=show_path),
            unsafe_allow_html=True,
        )
        st.markdown(
            '<div class="legend">'
            '<span style="background:#f59e0b;color:#1f2937">trace-back path</span>'
            '<span style="background:#16a34a;color:#fff">path start</span>'
            '<span style="background:#dc2626;color:#fff">path end / max</span>'
            '<span style="background:#dbeafe;color:#0f172a">\u2196 diagonal</span>'
            '<span style="background:#dbeafe;color:#0f172a">\u2191 up (gap in seq 2)</span>'
            '<span style="background:#dbeafe;color:#0f172a">\u2190 left (gap in seq 1)</span>'
            "</div>",
            unsafe_allow_html=True,
        )
        st.caption(
            f"Showing rows {r0}-{r1} and columns {c0}-{c1}. "
            "Arrows in each cell show the move(s) that produced the optimal value."
        )

        # ---------------- traceback --------------------------------------
        st.markdown("#### Trace-back path")
        rows = traceback_rows(res)
        if rows:
            df = pd.DataFrame(rows)
            st.dataframe(df, hide_index=True, **STRETCH,
                         height=min(420, 45 + 35 * min(len(df), 10)))
            st.caption(
                f"{len(rows)} steps traced backwards from cell "
                f"{res.end_cell} to {res.start_cell}."
            )
            if res.is_local:
                st.caption("Tracing stops at the first cell whose value is 0.")
            else:
                st.caption("Tracing stops at cell (0, 0).")
        else:
            st.write("No trace-back steps (empty alignment).")