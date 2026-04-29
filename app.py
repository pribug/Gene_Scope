import streamlit as st
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np
import plotly.express as px
import re
import gc

sns.set_style("whitegrid")
st.set_page_config(page_title="GeneScope", layout="wide")

st.markdown("""
<style>
[data-testid="stMetricLabel"] { font-size: 12px !important; }
[data-testid="stMetricValue"] { font-size: 18px !important; }
</style>
""", unsafe_allow_html=True)

st.title("GeneScope: Alzheimer's Gene Expression Dashboard")
st.caption("Explore differential gene expression across brain regions and conditions.")


# ============================================================
#  LOAD + OPTIMISE IN ONE PASS
#  Changes vs original:
#  A. float64 → float32 on Expression  (50 % saving on biggest col)
#  B. Low-cardinality string cols → category dtype (80-90 % saving)
#  C. Intermediates freed with gc.collect() before returning
#  D. global_data = data.copy() removed — was a full second copy
# ============================================================
@st.cache_data
def load_data():
    expression = pd.read_csv("data/processed/final_expression.csv", encoding="latin1")
    metadata   = pd.read_csv("data/processed/metadata.csv",        encoding="latin1")

    for col in metadata.columns:
        if metadata[col].dtype == "object":
            metadata[col] = metadata[col].str.replace("\xa0", "").str.strip()

    metadata["Condition"] = (
        metadata["Condition"]
        .str.lower().str.strip().str.replace("\u2019", "'")
        .replace({
            "alzheimer's disease": "AD",
            "alzheimers disease":  "AD",
            "ad":                  "AD",
            "control":             "Control",
        })
    )
    metadata["Sex"] = (
        metadata["Sex"]
        .str.lower().str.strip()
        .replace({"male": "Male", "female": "Female"})
    )

    def is_fake_date(g):
        return bool(re.match(r"^\d{1,2}-[A-Za-z]{3}$", str(g)))

    expression = expression[~expression["Gene"].apply(is_fake_date)]
    expression["Gene"] = expression["Gene"].str.upper()

    expression_long = expression.melt(
        id_vars=["Gene", "Gene Title"],
        var_name="Sample",
        value_name="Expression",
    )
    merged = expression_long.merge(metadata, on="Sample")

    # FIX A: halve memory on the Expression column
    merged["Expression"] = merged["Expression"].astype("float32")

    # FIX B: category dtype for repeated string columns
    for col in ["Gene", "Gene Title", "Sample", "Condition",
                "Sex", "Brain_region", "Cell_type"]:
        if col in merged.columns:
            merged[col] = merged[col].astype("category")

    # FIX C: free intermediates before returning
    del expression, expression_long, metadata
    gc.collect()

    return merged


data = load_data()
# FIX D: no global_data = data.copy() — use `data` directly everywhere


# ============================================================
#  HEAVY PRE-COMPUTATIONS — cached, run once at startup
# ============================================================
@st.cache_data
def compute_volcano(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for gene_name, gdf in df.groupby("Gene", observed=True):
        ctrl = gdf[gdf["Condition"] == "Control"]["Expression"].astype(float)
        ad   = gdf[gdf["Condition"] == "AD"     ]["Expression"].astype(float)
        if len(ctrl) >= 2 and len(ad) >= 2:
            ctrl_mean = float(ctrl.mean())
            ad_mean   = float(ad.mean())
            if ctrl_mean != 0:
                log2fc  = np.log2(ad_mean / ctrl_mean)
                _, pval = ttest_ind(ctrl, ad)
                pval    = max(pval, 1e-300)   # guard against log10(0)
                records.append({
                    "Gene":    gene_name,
                    "log2FC":  log2fc,
                    "-log10p": -np.log10(pval),
                })
    return pd.DataFrame(records)


@st.cache_data
def compute_top_genes(df: pd.DataFrame) -> pd.DataFrame:
    top = (
        df.groupby(["Gene", "Condition"], observed=True)["Expression"]
        .mean()
        .astype("float32")
        .unstack()
    )
    if "AD" in top.columns and "Control" in top.columns:
        top = top.dropna()
        top["Difference"]  = top["AD"] - top["Control"]
        top["Fold Change"] = top["AD"] / top["Control"]
    return top


with st.spinner("Pre-computing statistics… (first load only)"):
    volcano_base = compute_volcano(data)
    top_genes_df = compute_top_genes(data)


# ============================================================
#  SUMMARY METRICS
# ============================================================
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total Genes",   data["Gene"].nunique())
col2.metric("Total Samples", data["Sample"].nunique())
col3.metric("Conditions",    "Control vs AD")
col4.metric("Brain Regions", data["Brain_region"].nunique())

st.divider()

# ============================================================
#  GENE SEARCH
# ============================================================
st.header("Gene Search")

search    = st.text_input("Search gene")
all_genes = sorted(data["Gene"].astype(str).unique())
genes     = [g for g in all_genes if search.upper() in g] if search else all_genes
selected  = st.multiselect("Select up to 3 genes", genes)
gene      = selected[0] if len(selected) == 1 else None

st.divider()

# ============================================================
#  SIDEBAR FILTERS
# ============================================================
st.sidebar.header("Filters")

region    = st.sidebar.multiselect("Brain Region", sorted(data["Brain_region"].dropna().astype(str).unique()))
cell_type = st.sidebar.multiselect("Cell Type",    sorted(data["Cell_type"].dropna().astype(str).unique()))
condition = st.sidebar.multiselect("Condition",    sorted(data["Condition"].dropna().astype(str).unique()))
sex       = st.sidebar.multiselect("Sex",          sorted(data["Sex"].dropna().astype(str).unique()))

# FIX E: boolean mask approach — one slice instead of 4 chained DataFrame copies
mask = pd.Series(True, index=data.index)
if region:    mask &= data["Brain_region"].isin(region)
if cell_type: mask &= data["Cell_type"].isin(cell_type)
if condition: mask &= data["Condition"].isin(condition)
if sex:       mask &= data["Sex"].isin(sex)

filtered = data[mask]   # single slice, no intermediate copies


# ============================================================
#  MULTI-GENE COMPARISON
# ============================================================
if len(selected) > 1:
    st.subheader("Multi-Gene Comparison")

    multi = filtered[filtered["Gene"].isin(selected)]

    fig, ax = plt.subplots()
    sns.boxplot(x="Gene", y="Expression", hue="Condition", data=multi, ax=ax)
    st.pyplot(fig)
    plt.close(fig)

    mean = (
        multi.groupby(["Gene", "Condition"], observed=True)["Expression"]
        .mean().unstack().round(2)
    )
    st.dataframe(mean)

    if "AD" in mean.columns and "Control" in mean.columns:
        fc = mean.copy()
        fc["Fold Change"] = fc["AD"] / fc["Control"]
        st.subheader("Fold Change Comparison")
        st.dataframe(fc.round(2))


# ============================================================
#  SINGLE GENE ANALYSIS
# ============================================================
if gene:
    gene_data = filtered[filtered["Gene"] == gene]

    if not gene_data.empty:
        description  = gene_data.iloc[0]["Gene Title"]
        control      = gene_data[gene_data["Condition"] == "Control"]["Expression"].astype(float)
        disease      = gene_data[gene_data["Condition"] == "AD"     ]["Expression"].astype(float)

        control_mean = control.mean() if len(control) else None
        disease_mean = disease.mean() if len(disease) else None
        p_value      = None
        fold_change  = None

        if len(control) >= 2 and len(disease) >= 2:
            _, p_value = ttest_ind(control, disease)

        if control_mean and control_mean != 0 and disease_mean is not None:
            fold_change = disease_mean / control_mean

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Gene Overview")
            st.markdown(f"## {gene}")
            st.markdown(str(description))

            a, b, c, d = st.columns(4)
            a.metric("Control Mean", round(control_mean, 2) if control_mean else "NA")
            b.metric("AD Mean",      round(disease_mean, 2) if disease_mean else "NA")
            c.metric("P-value",      f"{p_value:.5f}"       if p_value is not None else "NA")
            d.metric("Fold Change",  round(fold_change, 2)  if fold_change else "NA")

            st.markdown("---")
            st.markdown("**Data Quality Check:**")
            st.write(f"Control samples: {len(control)}")
            st.write(f"AD samples: {len(disease)}")

            st.markdown("---")
            st.markdown("**Context:**")
            st.write(f"Samples: {len(gene_data)}")
            st.write(f"Regions: {gene_data['Brain_region'].nunique()}")
            st.write(f"Cell Types: {gene_data['Cell_type'].nunique()}")

            st.subheader("Interpretation")
            if len(control) == 0 and len(disease) == 0:
                st.warning("No data after filtering")
            elif len(control) == 0:
                st.warning("Only AD samples available")
            elif len(disease) == 0:
                st.warning("Only Control samples available")
            elif len(control) < 2 or len(disease) < 2:
                st.warning("Not enough samples for statistical test")
            elif p_value < 0.05:
                direction = "Upregulated" if disease_mean > control_mean else "Downregulated"
                st.success(f"{direction} in AD (FC: {round(fold_change, 2)})")
            else:
                st.info("No significant difference detected")

        with col2:
            st.subheader("Expression Analysis")

            fig1, ax1 = plt.subplots()
            sns.boxplot(x="Condition", y="Expression", data=gene_data, ax=ax1)
            st.pyplot(fig1)
            plt.close(fig1)

            fig2, ax2 = plt.subplots()
            sns.stripplot(x="Condition", y="Expression", data=gene_data, ax=ax2)
            st.pyplot(fig2)
            plt.close(fig2)


# ============================================================
#  GLOBAL INSIGHTS
# ============================================================
st.divider()
st.header("Global Insights")

rank_mode = st.radio("Rank by:", ["Difference", "Fold Change"])

if not top_genes_df.empty and rank_mode in top_genes_df.columns:
    top_display = top_genes_df.sort_values(rank_mode, ascending=False).head(10)
    st.subheader(f"Top 10 Genes by {rank_mode}")
    st.dataframe(top_display.round(2))

# Volcano — volcano_base pre-computed at startup; only add "Selected" here
st.subheader("Volcano Plot")

if not volcano_base.empty:
    volcano = volcano_base.copy()
    volcano["Selected"] = volcano["Gene"].isin(selected).map({True: "Yes", False: "No"})

    fig = px.scatter(
        volcano,
        x="log2FC",
        y="-log10p",
        color="Selected",
        color_discrete_map={"Yes": "#e63946", "No": "#adb5bd"},
        hover_data=["Gene"],
        title="Volcano Plot",
    )
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray",
                  annotation_text="p = 0.05")
    fig.add_vline(x= 1, line_dash="dash", line_color="gray")
    fig.add_vline(x=-1, line_dash="dash", line_color="gray")

    # FIX F: use_container_width deprecated in Streamlit 1.57 → use width=
    st.plotly_chart(fig, width="stretch")
else:
    st.info("Not enough data to render the volcano plot.")
