"""NetPharm Pro - network pharmacology workflow (Streamlit UI). Logic lives in netpharm_core.py."""
import html
from datetime import date

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import netpharm_core as C

st.set_page_config(page_title="NetPharm Pro - Network Pharmacology Web Tool", page_icon="🧬", layout="wide")
ss = st.session_state

st.markdown("""<style>
.gene-chip{display:inline-block;border:1px solid #475569;color:#38bdf8;border-radius:14px;padding:1px 9px;
margin:2px;font-size:.78rem;font-family:ui-monospace,monospace}
</style>""", unsafe_allow_html=True)

CONF = {"Low (0.15)": 0.15, "Medium (0.40) - STRING default": 0.40, "High (0.70)": 0.70,
        "Highest (0.90)": 0.90, "Custom": None}
ORG = {"Homo sapiens": 9606, "Mus musculus": 10090, "Rattus norvegicus": 10116}
LAYOUTS = ["Degree-sorted circular", "Hubs in centre (shells)", "Spring (force-directed)", "Kamada-Kawai"]
UPLOAD_TYPES = ["xlsx", "xls", "csv", "tsv", "txt"]

# ───────────────────────── sidebar: every threshold is a setting ─────────────────────────
S = {}
with st.sidebar:
    st.markdown("## 🧬 NetPharm Pro")
    st.caption("Compound → targets → disease overlap → STRING hubs → enrichment → pathways")
    with st.expander("⚙️ Analysis settings", expanded=True):
        S["swiss"] = st.slider("SwissTargetPrediction probability ≥", 0.0, 1.0, 0.0, 0.01,
                               help="0 keeps every predicted target.")
        S["ot"] = st.slider("Open Targets association score ≥", 0.0, 1.0, 0.10, 0.01)
        S["gc"] = st.number_input("GeneCards relevance score ≥", 0.0, 1000.0, 0.0, 0.5)
        S["combine"] = st.radio("Combine several sources / search terms", ["Union", "Intersection"], horizontal=True)
        choice = st.selectbox("STRING confidence", list(CONF), index=1)
        S["conf"] = CONF[choice] if CONF[choice] is not None else st.slider("Custom confidence", 0.0, 1.0, 0.40, 0.01)
        org = st.selectbox("Organism", list(ORG))
        S["org"], S["taxid"] = org, ORG[org]
        if S["taxid"] != 9606:
            st.caption("Enrichment libraries and disease sources are human-only.")
        S["extra"] = st.number_input("Max extra interactors to add", 0, 50, 0)
        S["hide"] = st.checkbox("Hide disconnected nodes", True)
        S["ntype"] = st.selectbox("STRING network type", ["functional", "physical"])
        S["hub_rule"] = st.radio("Hub rule", ["Top N by degree", "Degree above mean"])
        S["hub_n"] = st.number_input("N (for Top N; ties included)", 1, 100, 10)
        S["fdr"] = st.slider("Enrichment FDR ≤", 0.001, 1.0, 0.05, 0.001)
        S["enr_src"] = st.radio("Enrichment source", ["In-app (Enrichr)", "ShinyGO upload"])
        S["n_path"] = st.number_input("Pathways per source to select", 1, 30, 5)
    with st.expander("🔌 Connection check"):
        if st.button("Test web services"):
            for k, v in C.check_connections().items():
                st.write(f"**{k}**: {v}")
    st.caption(f"Session date: {date.today().isoformat()}")

st.markdown("# 🧬 NetPharm Pro")
st.caption("Network pharmacology workflow for any compound and any disease. Nothing is invented: if a web service "
           "fails you will see the error and can upload your own file instead.")


def dl(label, data, name, mime, key):
    st.download_button(label, data, file_name=name, mime=mime, key=key)


def chips(genes):
    st.markdown("".join(f'<span class="gene-chip">{html.escape(g)}</span>' for g in genes), unsafe_allow_html=True)


def get_shared():
    if ss.get("shared_override"):
        return ss["shared_override"]
    if "comp_df" in ss and "dise_df" in ss:
        return sorted(set(ss["comp_df"]["Gene"]) & set(ss["dise_df"]["Gene"]))
    return []


# ───────────────────────── reusable gene-source panel ─────────────────────────
def source_panel(kind):
    k = kind[:4]
    opts = (["Saved library", "SwissTargetPrediction file", "Other gene file", "Paste list"] if kind == "compound"
            else ["Saved library", "Open Targets (auto)", "UniProt (auto)", "GeneCards file", "Other gene file", "Paste list"])
    chosen = st.multiselect("Gene sources (combine as many as you like)", opts, default=[opts[1]], key=f"{k}_src")
    lib, files, terms, pasted = [], {}, [], ""
    if "Saved library" in chosen:
        names = C.list_library(kind)
        lib = st.multiselect("Saved lists", names, key=f"{k}_lib") if names else []
        if not names:
            st.info(f"No saved {kind} lists yet. Build one below and use 'Save to library'.")
    if "SwissTargetPrediction file" in chosen:
        files["SwissTargetPrediction"] = st.file_uploader("SwissTargetPrediction export", type=UPLOAD_TYPES, key=f"{k}_sw")
    if kind == "disease" and any(o in chosen for o in ("Open Targets (auto)", "UniProt (auto)")):
        terms = [t.strip() for t in st.text_area("Disease search terms (one per line)", key=f"{k}_terms",
                                                 placeholder="e.g. generalized anxiety disorder").splitlines() if t.strip()]
    if "GeneCards file" in chosen:
        files["GeneCards"] = st.file_uploader("GeneCards export", type=UPLOAD_TYPES, key=f"{k}_gc")
    if "Other gene file" in chosen:
        files["Uploaded file"] = st.file_uploader("Gene file (column named Gene / Gene Symbol)", type=UPLOAD_TYPES, key=f"{k}_of")
    if "Paste list" in chosen:
        pasted = st.text_area("Paste gene symbols", key=f"{k}_paste")
    std = st.checkbox("Standardize symbols to official names (MyGene.info)", True, key=f"{k}_std")

    if st.button(f"Build {kind} gene list", type="primary", key=f"{k}_go"):
        frames, errs = {}, []

        def add(name, fn):
            try:
                frames[name] = fn()
            except Exception as e:  # show every failure, never hide it
                errs.append(f"{name}: {e}")

        with st.spinner("Collecting genes..."):
            for n in lib:
                add(f"Library: {n}", lambda n=n: C.load_library(kind, n)[["Gene", "Score"]])
            for label, f in files.items():
                if f is None:
                    continue

                def parse(f=f, label=label):
                    df = C.read_upload(f.name, f.getvalue())
                    if label == "SwissTargetPrediction":
                        return C.apply_cutoff(C.parse_gene_table(df, "swiss"), S["swiss"])
                    if label == "GeneCards":
                        return C.apply_cutoff(C.parse_gene_table(df), S["gc"])
                    return C.parse_gene_table(df)
                add(label, parse)
            for t in terms:
                if "Open Targets (auto)" in chosen:
                    add(f"Open Targets: {t}", lambda t=t: C.apply_cutoff(C.opentargets_genes(t), S["ot"]))
                if "UniProt (auto)" in chosen:
                    add(f"UniProt: {t}", lambda t=t: C.uniprot_genes(t))
            if pasted.strip():
                add("Pasted list", lambda: pd.DataFrame({"Gene": C.clean_genes([pasted]), "Score": float("nan")}))
            changes, unmapped = pd.DataFrame(), []
            if frames and std:
                try:
                    frames, changes, unmapped = C.standardize_frames(frames, S["taxid"])
                except C.FetchError as e:
                    errs.append(f"Symbol standardization skipped ({e}). Symbols were only cleaned, not mapped.")
        if not frames:
            st.error("No genes collected. " + " | ".join(errs) if errs else "Choose a source and provide input first.")
            return
        ss[f"{k}_df"] = C.merge_sources(frames, S["combine"])
        ss[f"{k}_frames"] = {n: set(f["Gene"]) for n, f in frames.items()}
        ss[f"{k}_notes"] = (errs, changes, unmapped)
        ss.pop("ppi", None)
        ss.pop("shared_override", None)

    if f"{k}_df" in ss:
        df, (errs, changes, unmapped) = ss[f"{k}_df"], ss[f"{k}_notes"]
        for e in errs:
            st.warning(e)
        c1, c2 = st.columns(2)
        c1.metric("Unique genes", len(df))
        c2.metric("Sources used", len(ss[f"{k}_frames"]))
        st.dataframe(pd.DataFrame({"Source": list(ss[f"{k}_frames"]), "Genes": [len(v) for v in ss[f"{k}_frames"].values()]}),
                     hide_index=True)
        if len(changes):
            with st.expander(f"{len(changes)} symbols were renamed - please review"):
                st.dataframe(changes, hide_index=True)
        if unmapped:
            with st.expander(f"{len(unmapped)} entries were not recognised as human genes and were dropped"):
                st.write(", ".join(unmapped))
        st.dataframe(df, hide_index=True)
        dl("⬇ Gene list (CSV)", df.to_csv(index=False), f"{kind}_genes.csv", "text/csv", f"{k}_dl")
        name = st.text_input("Library name", value=ss.get("cname", "") if kind == "compound" else
                             (terms[0] if terms else ""), key=f"{k}_libname")
        cA, cB = st.columns(2)
        if cA.button("💾 Save to library", key=f"{k}_save") and name.strip():
            try:
                st.success(f"Saved: {C.save_library(kind, name, df)}. On Streamlit Cloud this is temporary - "
                           f"download the file and commit it to data/{kind}s/ in GitHub to keep it.")
            except OSError as e:
                st.error(f"Could not write to disk ({e}). Use the download and commit it to GitHub instead.")
        with cB:
            dl("⬇ Library file (commit to GitHub)", C.library_csv(df).to_csv(index=False),
               f"{C.safe_name(name or kind)}.csv", "text/csv", f"{k}_libdl")


tabs = st.tabs(["1 Compound", "2 Disease", "3 Overlap", "4 STRING PPI", "5 Enrichment",
                "6 Pathways", "7 Network editor", "8 Report"])

# ───────────────────────── 1 compound ─────────────────────────
with tabs[0]:
    st.subheader("Compound and predicted targets")
    c1, c2 = st.columns([3, 1])
    cname = c1.text_input("Compound / drug name", key="cname", placeholder="e.g. quercetin")
    c2.write("")
    if c2.button("Fetch from PubChem") and cname.strip():
        try:
            ss["pc"] = C.pubchem_info(cname)
        except Exception as e:
            ss.pop("pc", None)
            st.error(f"PubChem lookup failed: {e}")
    if "pc" in ss:
        p = ss["pc"]
        viol, checks = C.lipinski(p)
        m = st.columns(4)
        m[0].metric("Formula", p.get("MolecularFormula", "-"))
        m[1].metric("MW", p.get("MolecularWeight", "-"))
        m[2].metric("XLogP", p.get("XLogP", "-"))
        m[3].metric("Lipinski violations", viol)
        st.caption("Rule of five: " + "; ".join(f"{k} {'✓' if v else '✗'}" for k, v in checks.items()))
        st.code(p["SMILES_final"], language=None)
        st.link_button("Open SwissTargetPrediction (paste the SMILES, download the result)",
                       "http://www.swisstargetprediction.ch")
    source_panel("compound")

# ───────────────────────── 2 disease ─────────────────────────
with tabs[1]:
    st.subheader("Disease genes")
    st.caption("Tip: search the condition and its symptoms as separate terms, then choose Union or Intersection in settings.")
    source_panel("disease")

# ───────────────────────── 3 overlap ─────────────────────────
with tabs[2]:
    st.subheader("Common genes")
    if "comp_df" not in ss or "dise_df" not in ss:
        st.info("Build both the compound list (tab 1) and the disease list (tab 2) first.")
    else:
        A, B = set(ss["comp_df"]["Gene"]), set(ss["dise_df"]["Gene"])
        shared = get_shared()
        m = st.columns(3)
        m[0].metric("Compound targets", len(A))
        m[1].metric("Disease genes", len(B))
        m[2].metric("Common genes", len(shared))
        if shared:
            try:
                png = C.venn_bytes({"Compound targets": A, "Disease genes": B})
                st.image(png, width=480)
                dl("⬇ Venn PNG (300 dpi)", png, "venn.png", "image/png", "v_png")
                dl("⬇ Venn PDF", C.venn_bytes({"Compound targets": A, "Disease genes": B}, "pdf"), "venn.pdf",
                   "application/pdf", "v_pdf")
            except Exception as e:
                st.warning(f"Venn diagram unavailable: {e}")
            chips(shared)
            tbl = (pd.DataFrame({"Gene": shared})
                   .merge(ss["comp_df"].rename(columns={"Score": "Compound score", "Sources": "Compound sources"}), on="Gene")
                   .merge(ss["dise_df"].rename(columns={"Score": "Disease score", "Sources": "Disease sources"}), on="Gene"))
            st.dataframe(tbl, hide_index=True)
            dl("⬇ Common genes (CSV)", tbl.to_csv(index=False), "common_genes.csv", "text/csv", "sh_csv")
            ss["shared_tbl"] = tbl
        else:
            st.error("No common genes. Loosen a cutoff or add a source.")
        for side, key in (("compound", "comp_frames"), ("disease", "dise_frames")):
            fr = ss.get(key, {})
            if len(fr) in (2, 3):
                with st.expander(f"Venn of the {side} sources"):
                    try:
                        st.image(C.venn_bytes(fr), width=480)
                    except Exception as e:
                        st.warning(str(e))
        with st.expander("Use my own common-gene list instead (e.g. from Google Sheets)"):
            own = st.text_area("Paste genes", key="own_txt")
            if st.button("Use this list") and own.strip():
                ss["shared_override"] = C.clean_genes([own])
                st.rerun()
            if ss.get("shared_override") and st.button("Go back to the computed overlap"):
                ss.pop("shared_override")
                st.rerun()

# ───────────────────────── 4 STRING ─────────────────────────
with tabs[3]:
    st.subheader("STRING protein-protein interaction network")
    shared = get_shared()
    txt = st.text_area("Genes sent to STRING (editable)", "\n".join(shared), height=110, key="str_genes")
    genes = C.clean_genes([txt])
    mode = st.radio("Network source", ["Query STRING API", "Upload STRING short tabular output"], horizontal=True)
    up = st.file_uploader("STRING export (.tsv / .txt)", type=UPLOAD_TYPES, key="str_up") if mode.startswith("Upload") else None
    if st.button("Run STRING analysis", type="primary"):
        if len(genes) < 2:
            st.error("Need at least 2 genes.")
        else:
            try:
                with st.spinner("Building network..."):
                    if mode.startswith("Query"):
                        edges = C.string_network(genes, S["conf"], S["taxid"], S["extra"], S["ntype"])
                        pe = C.string_ppi_enrichment(genes, S["conf"], S["taxid"], S["extra"], S["ntype"])
                        ver = C.string_version()
                    else:
                        if up is None:
                            raise ValueError("Upload the STRING file first.")
                        edges = C.edges_from_table(C.read_upload(up.name, up.getvalue()), S["conf"])
                        pe, ver = None, "from uploaded file"
                    G = C.build_graph(genes, edges, S["hide"])
                ss["ppi"] = {"G": G, "edges": edges, "pe": pe, "ver": ver, "conf": S["conf"], "mode": mode}
            except Exception as e:
                st.error(f"STRING step failed: {e}")
    if "ppi" in ss:
        P = ss["ppi"]
        G = P["G"]
        M = C.node_metrics(G)
        hubs = C.select_hubs(M, S["hub_rule"], S["hub_n"])
        stt = C.network_stats(G)
        m = st.columns(5)
        for i, (k, v) in enumerate(stt.items()):
            m[i].metric(k, v)
        m[4].metric("PPI enrichment p", f"{P['pe']['p_value']:.2e}" if P["pe"] else "n/a")
        st.caption(f"STRING version: {P['ver']}; confidence ≥ {P['conf']}; source: {P['mode']}")
        if len(M):
            st.markdown(f"**Core hubs ({S['hub_rule']}, {len(hubs)} genes)**")
            chips(hubs)
            st.bar_chart(M.head(20).set_index("Gene")["Degree"])
            st.dataframe(M.assign(Hub=M["Gene"].isin(hubs)), hide_index=True)
            png = C.draw_network_png(G, C.layout_positions(G, LAYOUTS[0]), M, hubs=hubs)
            st.image(png, width=640)
            c = st.columns(3)
            with c[0]:
                dl("⬇ Network image (300 dpi PNG)", png, "string_network.png", "image/png", "p_png")
            with c[1]:
                dl("⬇ Interactions (TSV)", P["edges"].to_csv(sep="\t", index=False), "string_edges.tsv", "text/tab-separated-values", "p_tsv")
            with c[2]:
                dl("⬇ Degree table (CSV)", M.to_csv(index=False), "degree_table.csv", "text/csv", "p_deg")
        else:
            st.warning("No interactions at this confidence. Lower it in the settings.")

# ───────────────────────── 5 enrichment ─────────────────────────
with tabs[4]:
    st.subheader("KEGG and GO enrichment")
    genes = get_shared()
    if S["enr_src"].startswith("In-app"):
        libs = st.multiselect("Enrichr libraries", ["KEGG_2021_Human", "GO_Biological_Process_2023", "KEGG_2019_Human",
                              "GO_Biological_Process_2021", "Reactome_2022", "WikiPathway_2023_Human"],
                              default=["KEGG_2021_Human", "GO_Biological_Process_2023"])
        st.caption("Enrichr reports Benjamini-Hochberg adjusted p-values; they are shown in the FDR column.")
        if st.button("Run enrichment", type="primary"):
            try:
                with st.spinner("Contacting Enrichr..."):
                    ss["enr"] = C.enrichr(genes, libs)
            except Exception as e:
                st.error(f"Enrichment failed: {e}. Switch to 'ShinyGO upload' in the settings.")
    else:
        st.caption("Upload the full results table downloaded from ShinyGO (one file for KEGG, one for GO BP).")
        fk = st.file_uploader("ShinyGO KEGG table", type=UPLOAD_TYPES, key="sg_k")
        fg = st.file_uploader("ShinyGO GO Biological Process table", type=UPLOAD_TYPES, key="sg_g")
        if st.button("Read ShinyGO files", type="primary"):
            try:
                parts = [C.parse_shinygo(C.read_upload(f.name, f.getvalue()), lab)
                         for f, lab in ((fk, "KEGG (ShinyGO)"), (fg, "GO BP (ShinyGO)")) if f is not None]
                if not parts:
                    raise ValueError("Upload at least one file.")
                ss["enr"] = pd.concat(parts, ignore_index=True)
            except Exception as e:
                st.error(str(e))
    if "enr" in ss:
        E = ss["enr"][ss["enr"]["FDR"] <= S["fdr"]].sort_values("FDR")
        st.metric(f"Terms with FDR ≤ {S['fdr']}", len(E))
        for src, d in E.groupby("Source", sort=False):
            st.markdown(f"**{src}** ({len(d)} terms)")
            st.dataframe(d.drop(columns="Source"), hide_index=True)
            png = C.enrichment_bar_png(d, src)
            st.image(png, width=640)
            dl(f"⬇ {src} chart (PNG)", png, f"{C.safe_name(src)}.png", "image/png", f"e_{src}")
        ss["enr_f"] = E

# ───────────────────────── 6 pathways ─────────────────────────
with tabs[5]:
    st.subheader("Pathway-gene table")
    if "enr_f" not in ss or ss["enr_f"].empty:
        st.info("Run enrichment (tab 5) first.")
    else:
        picked = []
        for src, d in ss["enr_f"].groupby("Source", sort=False):
            terms = d["Term"].tolist()
            sel = st.multiselect(f"{src}", terms, default=terms[:int(S["n_path"])], key=f"pw_{src}")
            picked.append(d[d["Term"].isin(sel)])
        sel_df = pd.concat(picked) if picked else pd.DataFrame()
        if len(sel_df):
            long, wide = C.pathway_tables(sel_df)
            ss["pw_long"], ss["pw_wide"], ss["pw_sel"] = long, wide, sel_df
            st.markdown("**Pathway → genes (each pathway is a column)**")
            st.dataframe(wide, hide_index=True)
            c = st.columns(2)
            with c[0]:
                dl("⬇ Wide table (CSV)", wide.to_csv(index=False), "pathway_genes_wide.csv", "text/csv", "pw_w")
            with c[1]:
                dl("⬇ Long table (CSV)", long.to_csv(index=False), "pathway_genes_long.csv", "text/csv", "pw_l")
            png = C.pathway_graph_png(ss.get("cname") or "Compound", sel_df)
            st.image(png, width=720)
            dl("⬇ Compound-pathway-gene network (PNG)", png, "compound_pathway_gene.png", "image/png", "pw_png")

# ───────────────────────── 7 network editor ─────────────────────────
with tabs[6]:
    st.subheader("Interactive network editor")
    if "ppi" not in ss:
        st.info("Run the STRING analysis (tab 4) first.")
    else:
        G0 = ss["ppi"]["G"]
        r = st.columns(4)
        layout = r[0].selectbox("Layout", LAYOUTS)
        color_by = r[1].selectbox("Colour by", ["Degree", "Betweenness", "Closeness"])
        size_by = r[2].selectbox("Size by", ["Degree", "Betweenness", "Closeness"])
        cmap = r[3].selectbox("Colour map", ["YlOrRd", "viridis", "plasma", "coolwarm", "Blues"])
        f = st.columns(2)
        min_edge = f[0].slider("Minimum edge confidence", float(ss["ppi"]["conf"]), 1.0, float(ss["ppi"]["conf"]), 0.01)
        min_deg = f[1].number_input("Minimum node degree", 0, 50, 0)
        H = G0.copy()
        H.remove_edges_from([(a, b) for a, b, d in G0.edges(data=True) if d["score"] < min_edge])
        H.remove_nodes_from([n for n in list(H) if H.degree(n) < min_deg])
        if H.number_of_nodes() == 0:
            st.warning("Filters removed every node.")
        else:
            M = C.node_metrics(H)
            hubs = C.select_hubs(M, S["hub_rule"], S["hub_n"])
            pos = C.layout_positions(H, layout)
            try:
                components.html(C.pyvis_html(H, pos, M, color_by, size_by, cmap, hubs), height=720)
                st.caption("Drag nodes to rearrange, scroll to zoom. Dark outline = core hub. Hover for metrics.")
            except Exception as e:
                st.warning(f"Interactive view unavailable ({e}); showing a static image instead.")
                st.image(C.draw_network_png(H, pos, M, color_by, size_by, cmap, hubs), width=700)
            c = st.columns(4)
            with c[0]:
                dl("⬇ PNG (300 dpi)", C.draw_network_png(H, pos, M, color_by, size_by, cmap, hubs), "network.png", "image/png", "n_png")
            with c[1]:
                dl("⬇ GraphML (Cytoscape)", C.graphml_bytes(H, M), "network.graphml", "application/xml", "n_gml")
            with c[2]:
                dl("⬇ Node table (CSV)", M.to_csv(index=False), "nodes.csv", "text/csv", "n_nodes")
            with c[3]:
                ed = pd.DataFrame([(a, b, d["score"]) for a, b, d in H.edges(data=True)], columns=["source", "target", "score"])
                dl("⬇ Edge table (CSV)", ed.to_csv(index=False), "edges.csv", "text/csv", "n_edges")
            st.caption("In Cytoscape: File → Import → Network from File (GraphML), or import edges.csv then nodes.csv as a table.")

# ───────────────────────── 8 report ─────────────────────────
with tabs[7]:
    st.subheader("Report and methods log")
    P = ss.get("ppi")
    settings = pd.DataFrame([
        ("Date", date.today().isoformat()), ("Compound", ss.get("cname", "")), ("Organism", S["org"]),
        ("SwissTargetPrediction probability >=", S["swiss"]), ("Open Targets score >=", S["ot"]),
        ("GeneCards relevance >=", S["gc"]), ("Source combination", S["combine"]),
        ("STRING confidence", P["conf"] if P else S["conf"]), ("STRING version", P["ver"] if P else ""),
        ("STRING network type", S["ntype"]), ("Extra interactors", S["extra"]), ("Hide disconnected", S["hide"]),
        ("Hub rule", f"{S['hub_rule']} (N={S['hub_n']})"), ("Enrichment source", S["enr_src"]),
        ("Enrichment FDR <=", S["fdr"]), ("Pathways per source", S["n_path"]),
    ], columns=["Parameter", "Value"])
    sheets = {"Settings": settings}
    if "pc" in ss:
        sheets["Compound info"] = pd.DataFrame(ss["pc"].items(), columns=["Property", "Value"])
    for key, name in (("comp_df", "Compound genes"), ("dise_df", "Disease genes"), ("shared_tbl", "Common genes")):
        if key in ss:
            sheets[name] = ss[key]
    if P:
        M = C.node_metrics(P["G"])
        hubs = C.select_hubs(M, S["hub_rule"], S["hub_n"])
        sheets["STRING interactions"] = P["edges"]
        sheets["Degree and hubs"] = M.assign(Hub=M["Gene"].isin(hubs))
    for key, name in (("enr_f", "Enrichment"), ("pw_long", "Pathway-gene long"), ("pw_wide", "Pathway-gene wide")):
        if key in ss:
            sheets[name] = ss[key]
    st.dataframe(settings, hide_index=True)
    st.write("Sheets in the workbook:", ", ".join(sheets))
    dl("⬇ Full Excel report", C.build_report(sheets), f"netpharm_report_{date.today().isoformat()}.xlsx",
       "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "rep")
    if P and "comp_df" in ss and "dise_df" in ss:
        st.markdown("**Draft methods sentence (edit before use)**")
        st.text_area("Methods draft", (
            f"Targets of {ss.get('cname', 'the compound')} were collected from "
            f"{'; '.join(ss['comp_frames'])} and disease-associated genes from {'; '.join(ss['dise_frames'])}. "
            f"Overlapping genes ({len(get_shared())}) were submitted to STRING (v{P['ver']}, {S['org']}, confidence ≥ "
            f"{P['conf']}). Hubs were defined by degree ({S['hub_rule']}). Enrichment used {S['enr_src']} with FDR ≤ {S['fdr']}."),
            height=130, label_visibility="collapsed")
st.divider()
st.caption("For research and education. Predicted targets and associations need experimental validation.")
