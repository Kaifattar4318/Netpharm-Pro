"""NetPharm Pro - core logic (no Streamlit imports, so every function can be unit-tested).

Design rules
* Never invent data. If a web service fails we raise FetchError; the UI then offers an upload instead.
* Every threshold is a function argument, so the app can log exactly what was used.
"""
from __future__ import annotations

import glob
import io
import os
import re
from datetime import date
from urllib.parse import quote

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIRS = {
    "compound": os.path.join(BASE_DIR, "data", "compounds"),
    "disease": os.path.join(BASE_DIR, "data", "diseases"),
}
HEADERS = {"User-Agent": "NetPharmPro/1.0 (academic research tool)"}
TIMEOUT = 30
STRING_API = "https://string-db.org/api"
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
ENRICHR = "https://maayanlab.cloud/Enrichr"


class FetchError(Exception):
    """A web service could not be reached or returned nothing usable."""


def _req(method, url, **kw):
    kw.setdefault("timeout", TIMEOUT)
    kw["headers"] = {**HEADERS, **kw.get("headers", {})}
    try:
        r = requests.request(method, url, **kw)
        r.raise_for_status()
        return r
    except requests.RequestException as e:
        raise FetchError(f"{url.split('/')[2]} - {e.__class__.__name__}: {str(e)[:150]}") from e


# ───────────────────────── gene-symbol handling ─────────────────────────
_GENE_RE = re.compile(r"^[A-Z0-9][A-Z0-9\-\.@_]{0,24}$")
_STOP = {"GENE", "GENES", "SYMBOL", "NAME", "NA", "NAN", "NONE"}


def clean_genes(items, first_only=False):
    """Split, upper-case, de-duplicate (order kept) and drop anything that is not gene-shaped."""
    out = {}
    for x in items:
        if x is None or (isinstance(x, float) and pd.isna(x)):
            continue
        toks = [t for t in re.split(r"[;,\s/|]+", str(x).strip()) if t]
        if first_only:
            toks = toks[:1]
        for t in toks:
            t = t.upper()
            if t not in _STOP and _GENE_RE.match(t):
                out[t] = None
    return list(out)


def standardize_symbols(genes, taxid=9606):
    """Map aliases to official symbols with MyGene.info. Returns {entered: official or None}."""
    genes = list(dict.fromkeys(genes))
    mapping = {}
    for i in range(0, len(genes), 500):
        chunk = genes[i:i + 500]
        r = _req("POST", "https://mygene.info/v3/query",
                 data={"q": ",".join(chunk), "scopes": "symbol,alias", "fields": "symbol",
                       "species": str(taxid)})
        hits = {}
        for h in r.json():
            if h.get("notfound") or "symbol" not in h:
                continue
            q, sym = str(h.get("query", "")).upper(), str(h["symbol"]).upper()
            if q not in hits or sym == q:  # an exact symbol match always wins
                hits[q] = sym
        for g in chunk:
            mapping[g] = hits.get(g)
    return mapping


def standardize_frames(frames, taxid=9606):
    allg = sorted({g for f in frames.values() for g in f["Gene"]})
    mp = standardize_symbols(allg, taxid)
    new = {}
    for k, f in frames.items():
        d = f.copy()
        d["Gene"] = d["Gene"].map(mp)
        d = d.dropna(subset=["Gene"])
        new[k] = d.groupby("Gene", as_index=False, sort=False)["Score"].max()
    changes = pd.DataFrame([(a, b) for a, b in mp.items() if b and a != b],
                           columns=["Entered", "Official symbol"])
    unmapped = [a for a, b in mp.items() if not b]
    return new, changes, unmapped


# ───────────────────────── file parsing ─────────────────────────
def read_upload(name, data: bytes):
    ext = name.lower().rsplit(".", 1)[-1]
    if ext in ("xlsx", "xlsm", "xls"):
        return pd.read_excel(io.BytesIO(data))
    text = data.decode("utf-8-sig", errors="replace")
    if ext == "txt" and "\t" not in text:
        return pd.DataFrame({"Gene": clean_genes([text])})
    return pd.read_csv(io.StringIO(text), sep=None, engine="python")


def _find_col(df, *keys):
    low = {c: str(c).strip().lower().lstrip("#") for c in df.columns}
    for k in keys:  # exact match first, then substring
        for c, lc in low.items():
            if lc == k:
                return c
    for k in keys:
        for c, lc in low.items():
            if k in lc:
                return c
    return None


def parse_gene_table(df, kind="generic"):
    """Return DataFrame[Gene, Score]. kind='swiss' reads SwissTargetPrediction exports."""
    if kind == "swiss":
        gcol = _find_col(df, "common name", "gene symbol", "gene", "symbol")
        scol = _find_col(df, "probability")
    else:
        gcol = _find_col(df, "gene symbol", "gene", "symbol", "gene names", "gene name",
                         "approved symbol", "gene_primary")
        scol = _find_col(df, "relevance score", "relevance", "association score", "score")
    if gcol is None:
        if df.shape[1] == 1:
            gcol = df.columns[0]
        else:
            raise ValueError("Could not find a gene column. Rename it 'Gene' or 'Gene Symbol'.")
    scores = pd.to_numeric(df[scol], errors="coerce") if scol is not None else pd.Series(np.nan, index=df.index)
    first_only = "name" in str(gcol).lower()
    rows = [(g, s) for val, s in zip(df[gcol], scores) for g in clean_genes([val], first_only)[:1 if first_only else None]]
    if not rows:
        raise ValueError("No valid gene symbols found in the file.")
    out = pd.DataFrame(rows, columns=["Gene", "Score"])
    return out.groupby("Gene", as_index=False, sort=False)["Score"].max()


def apply_cutoff(df, cutoff):
    if cutoff is None or cutoff <= 0 or df["Score"].isna().all():
        return df
    return df[df["Score"] >= cutoff].reset_index(drop=True)


def merge_sources(frames, mode="Union"):
    cols = ["Gene", "Score", "Sources"]
    if not frames:
        return pd.DataFrame(columns=cols)
    parts = []
    for k, f in frames.items():
        d = f[["Gene", "Score"]].copy()
        d["Source"] = k
        parts.append(d)
    allr = pd.concat(parts, ignore_index=True)
    g = allr.groupby("Gene", sort=False).agg(
        Score=("Score", "max"),
        Sources=("Source", lambda s: "; ".join(dict.fromkeys(s))),
        n=("Source", "nunique")).reset_index()
    if mode == "Intersection":
        g = g[g["n"] == len(frames)]
    g = g.drop(columns="n").sort_values("Score", ascending=False, na_position="last", kind="stable")
    return g.reset_index(drop=True)[cols]


# ───────────────────────── offline library ─────────────────────────
def safe_name(s):
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", str(s).strip()).strip("_").lower() or "unnamed"


def list_library(kind):
    return sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(LIB_DIRS[kind], "*.csv")))


def load_library(kind, name):
    df = pd.read_csv(os.path.join(LIB_DIRS[kind], f"{name}.csv"))
    if "Score" not in df.columns:
        df["Score"] = np.nan
    df["Gene"] = df["Gene"].astype(str).str.upper().str.strip()
    return df


def library_csv(df, source_note=""):
    out = df[["Gene", "Score"]].copy()
    out["Source"] = df["Sources"] if "Sources" in df.columns else source_note
    out["Date"] = date.today().isoformat()
    return out


def save_library(kind, name, df, source_note=""):
    os.makedirs(LIB_DIRS[kind], exist_ok=True)
    path = os.path.join(LIB_DIRS[kind], f"{safe_name(name)}.csv")
    library_csv(df, source_note).to_csv(path, index=False)
    return path


# ───────────────────────── web sources ─────────────────────────
def pubchem_info(name):
    base = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(name.strip(), safe='')}/property/%s/JSON"
    common = "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount"
    try:
        r = _req("GET", base % f"SMILES,ConnectivitySMILES,{common}")
    except FetchError:
        r = _req("GET", base % f"CanonicalSMILES,IsomericSMILES,{common}")
    p = r.json()["PropertyTable"]["Properties"][0]
    p["SMILES_final"] = p.get("SMILES") or p.get("IsomericSMILES") or p.get("ConnectivitySMILES") or p.get("CanonicalSMILES")
    return p


def lipinski(p):
    """Return (violations, details) using Lipinski's rule of five."""
    checks = {
        "MW <= 500": float(p.get("MolecularWeight", 0)) <= 500,
        "H-bond donors <= 5": p.get("HBondDonorCount", 0) <= 5,
        "H-bond acceptors <= 10": p.get("HBondAcceptorCount", 0) <= 10,
        "XLogP <= 5": float(p.get("XLogP", 0) or 0) <= 5,
    }
    return sum(not v for v in checks.values()), checks


def opentargets_genes(term, size=500):
    q1 = ('query($q:String!){search(queryString:$q, entityNames:["disease"], page:{index:0,size:1})'
          '{hits{id name}}}')
    hits = _req("POST", OT_URL, json={"query": q1, "variables": {"q": term}}).json()["data"]["search"]["hits"]
    if not hits:
        raise FetchError(f"Open Targets found no disease matching '{term}'")
    q2 = ('query($id:String!,$n:Int!){disease(efoId:$id){name associatedTargets(page:{index:0,size:$n})'
          '{rows{target{approvedSymbol} score}}}}')
    d = _req("POST", OT_URL, json={"query": q2, "variables": {"id": hits[0]["id"], "n": size}}).json()["data"]["disease"]
    rows = [(r["target"]["approvedSymbol"].upper(), r["score"]) for r in d["associatedTargets"]["rows"]]
    return pd.DataFrame(rows, columns=["Gene", "Score"])


def uniprot_genes(term, size=500):
    params = {"query": f"({term}) AND (reviewed:true) AND (organism_id:9606)",
              "fields": "gene_primary", "format": "tsv", "size": size}
    df = pd.read_csv(io.StringIO(_req("GET", "https://rest.uniprot.org/uniprotkb/search", params=params).text), sep="\t")
    genes = clean_genes(df.iloc[:, 0]) if len(df) else []
    if not genes:
        raise FetchError(f"UniProt returned no reviewed human genes for '{term}'")
    return pd.DataFrame({"Gene": genes, "Score": np.nan})


# ───────────────────────── STRING / graph ─────────────────────────
def _sp(genes, score, taxid, extra, ntype):
    return {"identifiers": "\r".join(genes), "species": taxid, "required_score": int(round(score * 1000)),
            "add_white_nodes": int(extra), "network_type": ntype, "caller_identity": "netpharm_pro"}


def string_network(genes, score, taxid=9606, extra=0, ntype="functional"):
    r = _req("POST", f"{STRING_API}/tsv/network", data=_sp(genes, score, taxid, extra, ntype))
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    if df.empty or "preferredName_A" not in df.columns:
        return pd.DataFrame(columns=["Node A", "Node B", "Score"])
    e = df[["preferredName_A", "preferredName_B", "score"]].copy()
    e.columns = ["Node A", "Node B", "Score"]
    return e.drop_duplicates().reset_index(drop=True)


def string_ppi_enrichment(genes, score, taxid=9606, extra=0, ntype="functional"):
    try:
        return _req("POST", f"{STRING_API}/json/ppi_enrichment", data=_sp(genes, score, taxid, extra, ntype)).json()[0]
    except (FetchError, IndexError, ValueError):
        return None


def string_version():
    try:
        return _req("GET", f"{STRING_API}/json/version").json()[0]["string_version"]
    except (FetchError, IndexError, KeyError, ValueError):
        return "unknown"


def edges_from_table(df, min_score=0.0):
    """Read STRING's 'short tabular text output' (or the API tsv) into [Node A, Node B, Score]."""
    a = _find_col(df, "node1", "preferredname_a")
    b = _find_col(df, "node2", "preferredname_b")
    s = _find_col(df, "combined_score", "score")
    if a is None or b is None or s is None:
        raise ValueError("Not a STRING tabular export (need node1, node2 and combined_score columns).")
    e = pd.DataFrame({"Node A": df[a].astype(str), "Node B": df[b].astype(str),
                      "Score": pd.to_numeric(df[s], errors="coerce")}).dropna()
    return e[e["Score"] >= min_score].drop_duplicates().reset_index(drop=True)


def build_graph(genes, edges, hide_isolated=True):
    G = nx.Graph()
    G.add_nodes_from(genes)
    for a, b, s in edges[["Node A", "Node B", "Score"]].itertuples(index=False):
        a, b = str(a).upper(), str(b).upper()
        if a == b:
            continue
        s = float(s)
        if G.has_edge(a, b):
            G[a][b]["score"] = max(G[a][b]["score"], s)
        else:
            G.add_edge(a, b, score=s)
    if hide_isolated:
        G.remove_nodes_from(list(nx.isolates(G)))
    return G


def node_metrics(G):
    cols = ["Rank", "Gene", "Degree", "Betweenness", "Closeness"]
    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=cols)
    deg, bet, clo = dict(G.degree()), nx.betweenness_centrality(G), nx.closeness_centrality(G)
    df = pd.DataFrame({"Gene": list(G.nodes), "Degree": [deg[n] for n in G.nodes],
                       "Betweenness": [round(bet[n], 4) for n in G.nodes],
                       "Closeness": [round(clo[n], 4) for n in G.nodes]})
    df = df.sort_values(["Degree", "Betweenness"], ascending=False, kind="stable").reset_index(drop=True)
    df.insert(0, "Rank", range(1, len(df) + 1))
    return df


def network_stats(G):
    n, e = G.number_of_nodes(), G.number_of_edges()
    return {"Nodes": n, "Edges": e, "Average degree": round(2 * e / n, 2) if n else 0.0,
            "Isolated nodes": sum(1 for _ in nx.isolates(G))}


def select_hubs(df, rule="Top N by degree", n=10):
    if df.empty:
        return []
    if rule.startswith("Top"):
        n = max(1, min(int(n), len(df)))
        return df.loc[df["Degree"] >= df["Degree"].iloc[n - 1], "Gene"].tolist()  # ties included
    return df.loc[df["Degree"] > df["Degree"].mean(), "Gene"].tolist()


# ───────────────────────── enrichment ─────────────────────────
def enrichr(genes, libraries):
    r = _req("POST", f"{ENRICHR}/addList", files={"list": (None, "\n".join(genes)), "description": (None, "NetPharmPro")})
    uid = r.json()["userListId"]
    rows = []
    for lib in libraries:
        res = _req("GET", f"{ENRICHR}/enrich", params={"userListId": uid, "backgroundType": lib}).json().get(lib, [])
        for t in res:  # [rank, term, p, z, combined, [genes], adj_p, ...]
            rows.append({"Source": lib, "Term": t[1], "P-value": t[2], "FDR": t[6],
                         "Overlap": len(t[5]), "Genes": ", ".join(t[5])})
    return pd.DataFrame(rows, columns=["Source", "Term", "P-value", "FDR", "Overlap", "Genes"])


def parse_shinygo(df, label):
    p, f = _find_col(df, "pathway", "term"), _find_col(df, "enrichment fdr", "fdr")
    g = next((c for c in df.columns if str(c).strip().lower() == "genes"), None)
    if p is None or f is None or g is None:
        raise ValueError("ShinyGO file needs 'Pathway', 'Enrichment FDR' and 'Genes' columns "
                         "(use the full results download, not the chart).")
    genes = [", ".join(clean_genes([x])) for x in df[g]]
    return pd.DataFrame({"Source": label, "Term": df[p].astype(str),
                         "P-value": np.nan, "FDR": pd.to_numeric(df[f], errors="coerce"),
                         "Overlap": [len(x.split(", ")) if x else 0 for x in genes], "Genes": genes})


def pathway_tables(sel):
    long = pd.DataFrame([(r.Term, g) for r in sel.itertuples() for g in clean_genes([r.Genes])],
                        columns=["Pathway", "Gene"])
    wide = pd.DataFrame({r.Term: pd.Series(clean_genes([r.Genes])) for r in sel.itertuples()})
    return long, wide


# ───────────────────────── drawing / exports ─────────────────────────
def _png(fig, dpi=300, fmt="png"):
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return buf.getvalue()


def _scale(vals, lo, hi):
    v = np.asarray(vals, dtype=float)
    if v.size == 0 or v.max() == v.min():
        return np.full(v.shape, (lo + hi) / 2)
    return lo + (v - v.min()) / (v.max() - v.min()) * (hi - lo)


def layout_positions(G, name):
    if G.number_of_nodes() == 0:
        return {}
    order = sorted(G.nodes, key=lambda n: (-G.degree(n), n))
    if name == "Degree-sorted circular":
        return nx.circular_layout(order)
    if name == "Hubs in centre (shells)":
        k = max(1, len(order) // 5)
        shells = [order[:k], order[k:k * 3], order[k * 3:]]
        return nx.shell_layout(G, nlist=[s for s in shells if s])
    if name == "Kamada-Kawai" and G.number_of_edges() > 0:
        return nx.kamada_kawai_layout(G)
    return nx.spring_layout(G, seed=42, k=1.2 / max(1, np.sqrt(G.number_of_nodes())))


def draw_network_png(G, pos, metrics, color_by="Degree", size_by="Degree", cmap="YlOrRd", hubs=(), dpi=300):
    m = metrics.set_index("Gene")
    fig, ax = plt.subplots(figsize=(9, 7))
    fig.patch.set_facecolor("white")
    nodes = list(G.nodes)
    cv = [m.loc[n, color_by] for n in nodes]
    sz = _scale([m.loc[n, size_by] for n in nodes], 250, 1600)
    w = [0.4 + 3 * G[a][b].get("score", 0.5) for a, b in G.edges]
    nx.draw_networkx_edges(G, pos, ax=ax, width=w, edge_color="#94a3b8", alpha=0.6)
    sc = nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=nodes, node_color=cv, cmap=cmap, node_size=sz,
                                edgecolors=["#111827" if n in hubs else "white" for n in nodes],
                                linewidths=[2.5 if n in hubs else 0.8 for n in nodes])
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7, font_weight="bold")
    cb = fig.colorbar(sc, ax=ax, shrink=0.6)
    cb.set_label(color_by)
    ax.axis("off")
    return _png(fig, dpi)


def pyvis_html(G, pos, metrics, color_by="Degree", size_by="Degree", cmap="YlOrRd", hubs=(), height="700px"):
    from matplotlib import colormaps
    from matplotlib.colors import Normalize, to_hex
    from pyvis.network import Network

    m = metrics.set_index("Gene")
    net = Network(height=height, width="100%", bgcolor="#ffffff", font_color="#111827", cdn_resources="in_line")
    nodes = list(G.nodes)
    cv = np.array([m.loc[n, color_by] for n in nodes], dtype=float)
    norm = Normalize(cv.min() if len(cv) else 0, cv.max() if len(cv) and cv.max() > cv.min() else (cv.min() + 1 if len(cv) else 1))
    sz = _scale([m.loc[n, size_by] for n in nodes], 12, 45)
    for n, c, s in zip(nodes, cv, sz):
        x, y = pos[n]
        tip = (f"{n}\nDegree {int(m.loc[n, 'Degree'])}\nBetweenness {m.loc[n, 'Betweenness']}\n"
               f"Closeness {m.loc[n, 'Closeness']}")
        net.add_node(n, label=n, x=float(x) * 450, y=float(y) * 450, size=float(s), title=tip, physics=False,
                     color={"background": to_hex(colormaps[cmap](norm(c))), "border": "#111827" if n in hubs else "#ffffff"},
                     borderWidth=4 if n in hubs else 1)
    for a, b, d in G.edges(data=True):
        net.add_edge(a, b, width=0.5 + 4 * d.get("score", 0.5), color="#94a3b8", title=f"STRING score {d.get('score', 0):.3f}")
    net.toggle_physics(False)
    return net.generate_html()


def graphml_bytes(G, metrics):
    H = G.copy()
    m = metrics.set_index("Gene")
    for n in H.nodes:
        H.nodes[n]["degree"] = int(m.loc[n, "Degree"])
        H.nodes[n]["betweenness"] = float(m.loc[n, "Betweenness"])
        H.nodes[n]["closeness"] = float(m.loc[n, "Closeness"])
    buf = io.BytesIO()
    nx.write_graphml(H, buf)
    return buf.getvalue()


def venn_bytes(sets, fmt="png", dpi=300):
    from matplotlib_venn import venn2, venn3

    names, vals = list(sets), [set(v) for v in sets.values()]
    if len(vals) not in (2, 3) or any(not v for v in vals):
        raise ValueError("A Venn diagram needs 2 or 3 non-empty sets.")
    fig, ax = plt.subplots(figsize=(6, 5))
    fig.patch.set_facecolor("white")
    cols = ["#2563eb", "#f97316", "#10b981"][:len(vals)]
    (venn2 if len(vals) == 2 else venn3)(vals, set_labels=names, set_colors=cols, alpha=0.6, ax=ax)
    return _png(fig, dpi, fmt)


def enrichment_bar_png(df, title, top=10, dpi=300):
    d = df.sort_values("FDR").head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, max(3, 0.45 * len(d) + 1)))
    fig.patch.set_facecolor("white")
    ax.barh([t[:60] for t in d["Term"]], -np.log10(d["FDR"].clip(lower=1e-300)), color="#0f766e")
    ax.set_xlabel("-log10(FDR)")
    ax.set_title(title)
    return _png(fig, dpi)


def pathway_graph_png(compound, sel, dpi=300):
    G = nx.Graph()
    G.add_node(compound, t="c")
    for r in sel.itertuples():
        G.add_edge(compound, r.Term)
        G.nodes[r.Term]["t"] = "p"
        for g in clean_genes([r.Genes]):
            G.add_edge(r.Term, g)
            G.nodes[g].setdefault("t", "g")
    pos = nx.spring_layout(G, seed=7, k=0.9)
    col = {"c": "#10b981", "p": "#a78bfa", "g": "#38bdf8"}
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.patch.set_facecolor("white")
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#cbd5e1")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=[col[G.nodes[n]["t"]] for n in G],
                           node_size=[1500 if G.nodes[n]["t"] != "g" else 450 for n in G])
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=6)
    ax.axis("off")
    return _png(fig, dpi)


def build_report(sheets):
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as w:
        for name, df in sheets.items():
            if df is not None and len(df.columns):
                df.to_excel(w, sheet_name=re.sub(r"[\[\]\*\?/\\:]", "", name)[:31], index=False)
    return buf.getvalue()


def check_connections():
    tests = {
        "PubChem": ("GET", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/cids/JSON", {}),
        "STRING": ("GET", f"{STRING_API}/json/version", {}),
        "Enrichr": ("GET", f"{ENRICHR}/datasetStatistics", {}),
        "Open Targets": ("POST", OT_URL, {"json": {"query": "{ meta { name } }"}}),
        "UniProt": ("GET", "https://rest.uniprot.org/uniprotkb/search", {"params": {"query": "TP53", "size": 1, "format": "tsv"}}),
        "MyGene.info": ("GET", "https://mygene.info/v3/metadata", {}),
    }
    out = {}
    for n, (m, u, kw) in tests.items():
        try:
            _req(m, u, timeout=10, **kw)
            out[n] = "✅ reachable"
        except FetchError as e:
            out[n] = f"❌ {e}"
    return out
