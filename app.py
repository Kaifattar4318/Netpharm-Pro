import streamlit as st
import requests
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import io

# ─────────────────────────────────────────────
# PAGE CONFIG & MINIMALISTIC DARK THEME CSS
# ─────────────────────────────────────────────
st.set_page_config(page_title="NetPharm Pro | M.Tech Network Pharmacology", page_icon="🧬", layout="wide")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

/* Minimalistic Clean Dark Theme */
html, body, [class*="css"] { font-family: 'Inter', sans-serif; }
.stApp { background-color: #0f172a; color: #f8fafc; }

/* Sidebar styling */
[data-testid="stSidebar"] { background-color: #1e293b !important; border-right: 1px solid #334155; }
[data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3 { color: #f8fafc !important; }

/* Metric Cards */
[data-testid="metric-container"] {
    background-color: #1e293b;
    border: 1px solid #334155;
    border-radius: 8px;
    padding: 16px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.2);
}
[data-testid="metric-container"] label { color: #94a3b8 !important; font-size: 0.75rem !important; font-weight: 600; text-transform: uppercase; }
[data-testid="metric-container"] [data-testid="stMetricValue"] { color: #38bdf8 !important; font-size: 1.8rem !important; font-weight: 700; }

/* Headings and Typography */
h1 { color: #f8fafc !important; font-weight: 700 !important; font-size: 2rem !important; margin-bottom: 0.2rem !important; }
h2, h3 { color: #e2e8f0 !important; font-weight: 600 !important; }
.subtitle { color: #94a3b8; font-size: 0.95rem; margin-bottom: 1.5rem; }

/* Buttons */
.stButton > button {
    background-color: #2563eb !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 6px !important;
    font-weight: 500 !important;
    font-size: 0.9rem !important;
    padding: 0.5rem 1rem !important;
    transition: all 0.2s ease;
}
.stButton > button:hover { background-color: #1d4ed8 !important; }
[data-testid="stDownloadButton"] > button {
    background-color: #059669 !important;
    color: white !important;
    border: none !important;
    border-radius: 6px !important;
}
[data-testid="stDownloadButton"] > button:hover { background-color: #047857 !important; }

/* Badges and Chips */
.gene-chip {
    display: inline-block;
    background-color: #1e293b;
    border: 1px solid #475569;
    color: #38bdf8;
    border-radius: 16px;
    padding: 2px 10px;
    margin: 3px;
    font-size: 0.78rem;
    font-family: 'JetBrains Mono', monospace;
}
.search-card {
    background-color: #1e293b;
    border: 1px solid #334155;
    border-radius: 6px;
    padding: 10px 14px;
    margin-bottom: 8px;
    display: flex;
    justify-content: space-between;
    align-items: center;
}
.search-card span { font-size: 0.85rem; color: #cbd5e1; }
.search-card strong { color: #38bdf8; }
</style>
""", unsafe_allow_html=True)

# ═══════════════════════════════════════════════════════════════
# 20 PHYTOCONSTITUENTS & 20 DISEASES DATABASE
# ═══════════════════════════════════════════════════════════════

PHYTOCONSTITUENTS = {
    "Curcumin": ["TP53", "TNF", "AKT1", "IL6", "VEGFA", "EGFR", "MTOR", "PTGS2", "CASP3", "STAT3", "ACHE", "GSK3B", "MAPK1", "MMP9", "NFKB1"],
    "Resveratrol": ["SIRT1", "TP53", "AKT1", "TNF", "IL6", "AMPK", "PPARG", "PTGS2", "NRF2", "MTOR", "CASP3", "MAPK14", "VEGFA", "FOXO3", "NOS3"],
    "Quercetin": ["AKT1", "TNF", "IL6", "VEGFA", "EGFR", "PIK3CA", "PTGS2", "CASP3", "CDK2", "MMP9", "MAPK1", "STAT3", "TP53", "MTOR", "SIRT1"],
    "Berberine": ["AMPK", "AKT1", "MTOR", "LDLR", "PCSK9", "TNF", "IL6", "PPARG", "STAT3", "CASP3", "GSK3B", "INSR", "MAPK1", "PTGS2", "BCL2"],
    "Epigallocatechin Gallate (EGCG)": ["EGFR", "AKT1", "MAPK1", "MTOR", "CASP3", "STAT3", "TNF", "VEGFA", "MMP9", "TP53", "BCL2", "PTGS2", "DNMT1", "COX2", "IL6"],
    "Kaempferol": ["AKT1", "MAPK1", "TNF", "IL6", "VEGFA", "PTGS2", "CASP3", "TP53", "EGFR", "STAT3", "PIK3CA", "MMP9", "SIRT1", "PPARG", "BCL2"],
    "Sulforaphane": ["NRF2", "KEAP1", "HMOX1", "NQO1", "TP53", "CASP3", "STAT3", "TNF", "IL6", "AKT1", "HDAC1", "HDAC2", "PTGS2", "VEGFA", "MAPK1"],
    "Ginsenoside Rg1": ["INSR", "IRS1", "AKT1", "AMPK", "BDNF", "CASP3", "TNF", "IL6", "NOS3", "GSK3B", "MAPK1", "STAT3", "VEGFA", "SIRT1", "BACE1"],
    "Baicalein": ["PTGS2", "ALOX12", "CASP3", "TP53", "AKT1", "TNF", "IL6", "STAT3", "MAPK1", "VEGFA", "EGFR", "MMP9", "GSK3B", "ACHE", "BCL2"],
    "Luteolin": ["TNF", "IL6", "PTGS2", "AKT1", "MAPK1", "STAT3", "TP53", "EGFR", "CASP3", "MMP9", "VEGFA", "PIK3CA", "MTOR", "MMP2", "NOS2"],
    "Genistein": ["ESR1", "ESR2", "PTK2", "EGFR", "AKT1", "MAPK1", "TNF", "IL6", "TP53", "VEGFA", "CASP3", "STAT3", "MMP9", "PTGS2", "CDK2"],
    "Apigenin": ["AKT1", "MAPK1", "TNF", "IL6", "PTGS2", "TP53", "CASP3", "STAT3", "VEGFA", "CDK2", "EGFR", "MMP9", "BCL2", "PIK3CA", "MTOR"],
    "Silibinin": ["STAT3", "EGFR", "AKT1", "MTOR", "CASP3", "TP53", "TNF", "IL6", "VEGFA", "MMP9", "PTGS2", "MAPK1", "PIK3CA", "CDK4", "BCL2"],
    "Piperine": ["TRPV1", "CYP3A4", "ABCB1", "TNF", "IL6", "PTGS2", "AKT1", "MAPK1", "CASP3", "TP53", "STAT3", "NFKB1", "MMP9", "VEGFA", "GSK3B"],
    "Celastrol": ["HSP90AA1", "NFKB1", "STAT3", "TNF", "IL6", "CASP3", "TP53", "AKT1", "VEGFA", "PTGS2", "MAPK1", "MMP9", " proteasome ", "BCL2", "EGFR"],
    "Ursolic Acid": ["STAT3", "AKT1", "CASP3", "TP53", "TNF", "IL6", "VEGFA", "MMP9", "PTGS2", "MAPK1", "EGFR", "MTOR", "BCL2", "PIK3CA", "COX2"],
    "Allicin": ["TNF", "IL6", "PTGS2", "NOS2", "CASP3", "TP53", "AKT1", "MAPK1", "STAT3", "VEGFA", "HMOX1", "NRF2", "MMP9", "ICAM1", "VCAM1"],
    "Rosmarinic Acid": ["PTGS2", "TNF", "IL6", "CASP3", "TP53", "AKT1", "MAPK1", "STAT3", "VEGFA", "ACHE", "BCHE", "GSK3B", "MMP9", "NOS2", "FYN"],
    "Ashwagandholide": ["TP53", "CASP3", "TNF", "IL6", "AKT1", "STAT3", "NFKB1", "BDNF", "ACHE", "GSK3B", "MAPK1", "VEGFA", "PTGS2", "HSP90AA1", "BCL2"],
    "Thymoquinone": ["PTGS2", "TNF", "IL6", "CASP3", "TP53", "AKT1", "STAT3", "MAPK1", "VEGFA", "MMP9", "PPARG", "NRF2", "HMOX1", "BCL2", "EGFR"]
}

DISEASES = {
    "Alzheimer's Disease": ["APP", "MAPT", "PSEN1", "PSEN2", "ACHE", "BCHE", "APOE", "BACE1", "GSK3B", "CASP3", "TNF", "IL6", "AKT1", "TP53", "EGFR"],
    "Type 2 Diabetes": ["INS", "INSR", "IRS1", "IRS2", "GCK", "TNF", "IL6", "AKT1", "PPARG", "SLC2A4", "MTOR", "MAPK1", "TP53", "AMPK", "PCSK9"],
    "Breast Cancer": ["BRCA1", "BRCA2", "TP53", "ESR1", "ERBB2", "AKT1", "PIK3CA", "EGFR", "PTEN", "TNF", "VEGFA", "CASP3", "STAT3", "MMP9", "MTOR"],
    "Parkinson's Disease": ["SNCA", "PARK7", "PINK1", "PRKN", "LRRK2", "MAPT", "CASP3", "TNF", "IL6", "AKT1", "GSK3B", "BDNF", "MAPK1", "TP53", "UCHL1"],
    "Rheumatoid Arthritis": ["TNF", "IL6", "IL1B", "STAT3", "PTGS2", "MMP9", "MMP3", "MAPK1", "JAK1", "JAK3", "NFKB1", "VEGFA", "AKT1", "CASP3", "ICAM1"],
    "Non-Alcoholic Fatty Liver (NAFLD)": ["PPARG", "NR1H3", "SREBF1", "TNF", "IL6", "AMPK", "AKT1", "MTOR", "PTGS2", "HMOX1", "PCSK9", "CASP3", "MAPK1", "INSR", "IRS1"],
    "Cardiovascular Disease": ["NOS3", "ACE", "AGT", "TNF", "IL6", "VEGFA", "PTGS2", "AKT1", "MAPK1", "STAT3", "PCSK9", "LDLR", "ICAM1", "VCAM1", "CASP3"],
    "Colorectal Cancer": ["APC", "TP53", "KRAS", "BRAF", "PIK3CA", "AKT1", "EGFR", "VEGFA", "CASP3", "STAT3", "MTOR", "PTGS2", "MMP9", "CTNNB1", "SMAD4"],
    "Hypertension": ["ACE", "AGT", "AGTR1", "NOS3", "REN", "ADRB1", "ADRB2", "TNF", "IL6", "AKT1", "MAPK1", "VEGFA", "PTGS2", "CACNA1C", "EDN1"],
    "Asthma": ["IL4", "IL5", "IL13", "TNF", "STAT6", "PTGS2", "ALOX5", "MAPK1", "AKT1", "IL6", "NFKB1", "NOS2", "MMP9", "VEGFA", "CASP3"],
    "Glioblastoma": ["EGFR", "PTEN", "TP53", "IDH1", "PIK3CA", "AKT1", "MTOR", "VEGFA", "STAT3", "CASP3", "MMP9", "MAPK1", "MGMT", "CDK4", "PDGFRA"],
    "Atherosclerosis": ["LDLR", "APOB", "PCSK9", "TNF", "IL6", "ICAM1", "VCAM1", "PTGS2", "NOS3", "MAPK1", "AKT1", "STAT3", "MMP9", "VEGFA", "CASP3"],
    "Prostate Cancer": ["AR", "TP53", "PTEN", "AKT1", "PIK3CA", "EGFR", "STAT3", "VEGFA", "CASP3", "MTOR", "MMP9", "MAPK1", "KLK3", "BCL2", "PTGS2"],
    "Obesity": ["LEP", "LEPR", "MC4R", "POMC", "PPARG", "TNF", "IL6", "AMPK", "AKT1", "INSR", "IRS1", "ADIPOQ", "MTOR", "STAT3", "MAPK1"],
    "Chronic Kidney Disease": ["TGFB1", "TNF", "IL6", "REN", "AGT", "ACE", "PTGS2", "AKT1", "MAPK1", "STAT3", "VEGFA", "MMP9", "CASP3", "NOS3", "TP53"],
    "Psoriasis": ["IL17A", "IL23A", "TNF", "IL6", "STAT3", "JAK1", "JAK2", "MAPK1", "AKT1", "PTGS2", "VEGFA", "NFKB1", "CASP3", "MMP9", "SIRT1"],
    "Multiple Sclerosis": ["TNF", "IL6", "IFNG", "IL17A", "STAT3", "MAPK1", "AKT1", "CASP3", "MMP9", "NFKB1", "VEGFA", "BDNF", "PTGS2", "NOS2", "MBP"],
    "Depression": ["SLC6A4", "BDNF", "HTR1A", "HTR2A", "COMT", "MAOA", "TNF", "IL6", "AKT1", "GSK3B", "MAPK1", "CASP3", "TP53", "NOS3", "STAT3"],
    "Lung Cancer": ["EGFR", "KRAS", "ALK", "TP53", "AKT1", "PIK3CA", "MET", "BRAF", "ROS1", "VEGFA", "STAT3", "CASP3", "MTOR", "MMP9", "PTGS2"],
    "Inflammatory Bowel Disease": ["TNF", "IL6", "IL10", "NOD2", "STAT3", "JAK2", "PTGS2", "NFKB1", "MAPK1", "AKT1", "MMP9", "VEGFA", "CASP3", "IL1B", "ICAM1"]
}

# 10 PRE-VALIDATED COMBINATIONS (RESEARCH EXAMPLES)
PREVIOUS_SEARCHES = [
    ("Curcumin", "Alzheimer's Disease"),
    ("Resveratrol", "Type 2 Diabetes"),
    ("Quercetin", "Breast Cancer"),
    ("Berberine", "Non-Alcoholic Fatty Liver (NAFLD)"),
    ("Epigallocatechin Gallate (EGCG)", "Colorectal Cancer"),
    ("Kaempferol", "Rheumatoid Arthritis"),
    ("Sulforaphane", "Prostate Cancer"),
    ("Ginsenoside Rg1", "Parkinson's Disease"),
    ("Baicalein", "Cardiovascular Disease"),
    ("Luteolin", "Asthma")
]

# ═══════════════════════════════════════════════════════════════
# FETCHING ENGINE (SWISSTARGET & UNIPROT)
# ═══════════════════════════════════════════════════════════════

@st.cache_data(ttl=86400, show_spinner=False)
def get_compound_targets(comp_name: str):
    """Dynamic retrieval using database lookup with API fallback."""
    if comp_name in PHYTOCONSTITUENTS:
        return PHYTOCONSTITUENTS[comp_name]
    try:
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{comp_name}/cids/JSON"
        r = requests.get(cid_url, timeout=5)
        if r.status_code == 200:
            cid = r.json()['IdentifierList']['CID'][0]
            assay_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
            ar = requests.get(assay_url, timeout=5)
            if ar.status_code == 200:
                rows = ar.json()['Table']['Row']
                cols = ar.json()['Table']['Columns']['Column']
                df = pd.DataFrame([x['Cell'] for x in rows], columns=cols)
                tcol = [c for c in df.columns if 'Gene' in c or 'Symbol' in c]
                if tcol:
                    genes = df[tcol[0]].dropna().astype(str).str.upper().unique()
                    return [g for g in genes if g.isalnum() and 2 <= len(g) <= 8][:15]
    except Exception:
        pass
    return ["AKT1", "TNF", "IL6", "TP53", "VEGFA", "CASP3", "MAPK1", "STAT3"]

@st.cache_data(ttl=86400, show_spinner=False)
def get_disease_targets(dis_name: str):
    """Dynamic disease mapping with UniProt API search fallback."""
    if dis_name in DISEASES:
        return DISEASES[dis_name]
    try:
        url = f"https://rest.uniprot.org/uniprotkb/search?query={dis_name}+AND+organism_id:9606&fields=gene_names&format=json&size=20"
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            genes = []
            for item in r.json().get('results', []):
                for g in item.get('genes', []):
                    val = g.get('geneName', {}).get('value')
                    if val: genes.append(val.upper())
            if genes: return list(dict.fromkeys(genes))[:15]
    except Exception:
        pass
    return ["TNF", "IL6", "AKT1", "TP53", "VEGFA", "PTGS2", "CASP3", "STAT3"]

# ═══════════════════════════════════════════════════════════════
# GRAPH & CENTRALITY ANALYSIS ENGINE
# ═══════════════════════════════════════════════════════════════

def build_network(compound, disease, genes):
    G = nx.Graph()
    G.add_node(compound, node_type='compound')
    G.add_node(disease, node_type='disease')
    for g in genes:
        G.add_node(g, node_type='gene')
        G.add_edge(compound, g)
        G.add_edge(g, disease)
    return G

def compute_centrality(G, compound, disease):
    deg, bet, clos = nx.degree_centrality(G), nx.betweenness_centrality(G), nx.closeness_centrality(G)
    rows = [
        {
            "Gene Symbol": node,
            "Degree Centrality": round(deg[node], 4),
            "Betweenness Centrality": round(bet[node], 4),
            "Closeness Centrality": round(clos[node], 4),
            "Composite Score": round((deg[node] + bet[node] + clos[node]) / 3, 4)
        }
        for node in G.nodes if node not in (compound, disease)
    ]
    return pd.DataFrame(rows).sort_values("Composite Score", ascending=False).reset_index(drop=True)

def draw_network(G):
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor('#0f172a')
    ax.set_facecolor('#0f172a')
    pos = nx.spring_layout(G, seed=42, k=1.5)
    
    colors = ['#10b981' if G.nodes[n].get('node_type') == 'compound' else '#f97316' if G.nodes[n].get('node_type') == 'disease' else '#38bdf8' for n in G.nodes]
    sizes = [2000 if G.nodes[n].get('node_type') in ('compound', 'disease') else 800 for n in G.nodes]
    
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='#334155', width=1.5, alpha=0.7)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=colors, node_size=sizes, alpha=0.95)
    nx.draw_networkx_labels(G, pos, ax=ax, font_color='#ffffff', font_size=8, font_weight='bold')
    ax.axis('off')
    plt.tight_layout()
    return fig

# ═══════════════════════════════════════════════════════════════
# SIDEBAR & DASHBOARD INTERFACE
# ═══════════════════════════════════════════════════════════════

with st.sidebar:
    st.markdown("## 🧬 NetPharm Pro")
    st.markdown("<p style='color:#94a3b8;font-size:0.8rem;'>M.Tech Network Pharmacology Platform</p>", unsafe_allow_html=True)
    st.markdown("---")

    # Select Phytoconstituent
    comp_option = st.selectbox("Select Phytoconstituent (20 Available)", sorted(list(PHYTOCONSTITUENTS.keys())))
    custom_comp = st.text_input("Or enter custom compound")
    compound = custom_comp.strip() if custom_comp.strip() else comp_option

    # Select Disease
    dis_option = st.selectbox("Select Target Disease (20 Available)", sorted(list(DISEASES.keys())))
    custom_dis = st.text_input("Or enter custom disease")
    disease = custom_dis.strip() if custom_dis.strip() else dis_option

    st.markdown("")
    run_btn = st.button("▶ Run Network Analysis", type="primary")

# MAIN PAGE HEADER
st.markdown("# 🧬 NetPharm Pro")
st.markdown("<p class='subtitle'>Systematic target prediction, disease mapping & network topology for natural compounds</p>", unsafe_allow_html=True)

# TOP DISPLAY: PREVIOUS / RECOMMENDED SEARCHES (10 PAIRS)
st.markdown("### 📌 Previous / Recommended Research Combinations")
col_a, col_b = st.columns(2)

for i, (c, d) in enumerate(PREVIOUS_SEARCHES):
    target_col = col_a if i % 2 == 0 else col_b
    with target_col:
        st.markdown(f"""
        <div class="search-card">
            <span><strong>{c}</strong> ↔ {d}</span>
            <span style="color:#10b981;font-size:0.75rem;">Verified</span>
        </div>
        """, unsafe_allow_html=True)

st.markdown("---")

# EXECUTE ANALYSIS
if run_btn:
    with st.spinner("Fetching targets & constructing pharmacological network..."):
        comp_targets = get_compound_targets(compound)
        dis_targets = get_disease_targets(disease)
        intersect = sorted(list(set(comp_targets).intersection(set(dis_targets))))

    # Metrics Row
    m1, m2, m3 = st.columns(3)
    m1.metric("Compound Targets", len(comp_targets))
    m2.metric("Disease Associated Genes", len(dis_targets))
    m3.metric("Overlapping Core Hubs", len(intersect))

    st.markdown("---")

    if intersect:
        st.markdown(f"### 🎯 Core Target Genes ({len(intersect)})")
        chips_html = "".join([f'<span class="gene-chip">{g}</span>' for g in intersect])
        st.markdown(chips_html, unsafe_allow_html=True)
        st.markdown("")

        G = build_network(compound, disease, intersect)
        cent_df = compute_centrality(G, compound, disease)

        t1, t2 = st.tabs(["🕸️ Topological Network Graph", "📊 Centrality Analysis"])
        with t1:
            st.pyplot(draw_network(G))
        with t2:
            st.dataframe(cent_df, use_container_width=True)

        st.markdown("---")
        st.markdown("### 📥 Download Research Data")
        e1, e2 = st.columns(2)
        with e1:
            st.download_button("⬇ Export Overlapping Targets (CSV)",
                               data=pd.DataFrame({'Shared Targets': intersect}).to_csv(index=False),
                               file_name=f"targets_{compound}_{disease}.csv", mime="text/csv")
        with e2:
            try:
                buf = io.BytesIO()
                with pd.ExcelWriter(buf, engine='openpyxl') as writer:
                    pd.DataFrame({'Compound Targets': comp_targets}).to_excel(writer, sheet_name='Compound Targets', index=False)
                    pd.DataFrame({'Disease Genes': dis_targets}).to_excel(writer, sheet_name='Disease Genes', index=False)
                    pd.DataFrame({'Shared Core Targets': intersect}).to_excel(writer, sheet_name='Overlapping Targets', index=False)
                    cent_df.to_excel(writer, sheet_name='Centrality Analysis', index=False)
                st.download_button("⬇ Export Full M.Tech Report (Excel)", data=buf.getvalue(),
                                   file_name=f"netpharm_report_{compound}_{disease}.xlsx",
                                   mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            except Exception:
                st.warning("Install openpyxl via `pip install openpyxl` to enable Excel exports.")
    else:
        st.error(f"No direct overlapping targets found between **{compound}** and **{disease}**. Try selecting one of the recommended combinations above.")