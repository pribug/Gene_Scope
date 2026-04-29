import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")   # to be set before any other plt import
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np
import plotly.express as px
import re

sns.set_style("whitegrid")
st.set_page_config(page_title="GeneScope", layout="wide")

# ------------------ STYLE ------------------
st.markdown("""
<style>
[data-testid="stMetricLabel"] { font-size: 12px !important; }
[data-testid="stMetricValue"] { font-size: 18px !important; }
</style>
""", unsafe_allow_html=True)

# ------------------ TITLE ------------------
st.title("GeneScope: Alzheimer's Gene Expression Dashboard")
st.caption("Explore differential gene expression across brain regions and conditions.")

# ==================== DATA LOADING ====================
# FIX 1: @st.cache_data ensures the CSV is read only ONCE per session,
#         not on every widget interaction.
@st.cache_data
def load_data():
    expression = pd.read_csv("data/processed/final_expression.csv", encoding="latin1")
    metadata   = pd.read_csv("data/processed/metadata.csv",        encoding="latin1")

    # --- clean metadata ---
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

    # --- fix Excel-auto-formatted gene names ---
    def is_fake_date(g):
        return bool(re.match(r"^\d{1,2}-[A-Za-z]{3}$", str(g)))

    expression = expression[~expression["Gene"].apply(is_fake_date)]
    expression["Gene"] = expression["Gene"].str.upper()

    # --- melt to long format ---
    expression_long = expression.melt(
        id_vars=["Gene", "Gene Title"],
        var_name="Sample",
        value_name="Expression",
    )
    return expression_long.merge(metadata, on="Sample")


data = load_data()
global_data = data.copy()

# ==================== PRE-COMPUTATIONS ====================

@st.cache_data
def compute_volcano(df: pd.DataFrame) -> pd.DataFrame:
    """Vectorised t-test over all genes. Returns one row per gene."""
    records = []
    for gene_name, gdf in df.groupby("Gene"):
        ctrl = gdf[gdf["Condition"] == "Control"]["Expression"].astype(float)
        ad   = gdf[gdf["Condition"] == "AD"     ]["Expression"].astype(float)
        if len(ctrl) >= 2 and len(ad) >= 2:
            ctrl_mean = ctrl.mean()
            ad_mean   = ad.mean()
            if ctrl_mean != 0:
                fc      = ad_mean / ctrl_mean
                log2fc  = np.log2(fc)
                _, pval = ttest_ind(ctrl, ad)
                pval    = max(pval, 1e-300)
                records.append({
                    "Gene":     gene_name,
                    "log2FC":   log2fc,
                    "-log10p":  -np.log10(pval),
                })
    return pd.DataFrame(records)


@st.cache_data
def compute_top_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Mean expression per gene×condition. Cached to avoid per-interaction recompute."""
    top = df.groupby(["Gene", "Condition"])["Expression"].mean().unstack()
    if "AD" in top.columns and "Control" in top.columns:
        top = top.dropna()
        top["Difference"]   = top["AD"] - top["Control"]
        top["Fold Change"]  = top["AD"] / top["Control"]
    return top


# Run heavy computations once and cached after first run
with st.spinner("Pre-computing statistics for all genes… (first load only)"):
    volcano_base = compute_volcano(global_data)
    top_genes_df = compute_top_genes(global_data)

# ==================== SUMMARY METRICS ====================
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total Genes",    len(data["Gene"].unique()))
col2.metric("Total Samples",  len(data["Sample"].unique()))
col3.metric("Conditions",     "Control vs AD")
col4.metric("Brain Regions",  data["Brain_region"].nunique())

st.divider()

# ==================== GENE SEARCH ====================
st.header("Gene Search")

search = st.text_input("Search gene")
genes  = sorted(data["Gene"].unique())
if search:
    genes = [g for g in genes if search.upper() in g]

selected = st.multiselect("Select up to 3 genes", genes)
gene     = selected[0] if len(selected) == 1 else None

st.divider()

# ==================== SIDEBAR FILTERS ====================
st.sidebar.header("Filters")

region    = st.sidebar.multiselect("Brain Region", sorted(data["Brain_region"].dropna().unique()))
cell_type = st.sidebar.multiselect("Cell Type",    sorted(data["Cell_type"].dropna().unique()))
condition = st.sidebar.multiselect("Condition",    sorted(data["Condition"].dropna().unique()))
sex       = st.sidebar.multiselect("Sex",          sorted(data["Sex"].dropna().unique()))

filtered = data.copy()
if region:    filtered = filtered[filtered["Brain_region"].isin(region)]
if cell_type: filtered = filtered[filtered["Cell_type"].isin(cell_type)]
if condition: filtered = filtered[filtered["Condition"].isin(condition)]
if sex:       filtered = filtered[filtered["Sex"].isin(sex)]

# ==================== MULTI-GENE COMPARISON ====================
if len(selected) > 1:
    st.subheader("Multi-Gene Comparison")

    multi = filtered[filtered["Gene"].isin(selected)]

    fig, ax = plt.subplots()
    sns.boxplot(x="Gene", y="Expression", hue="Condition", data=multi, ax=ax)
    st.pyplot(fig)
    plt.close(fig) 

    mean = multi.groupby(["Gene", "Condition"])["Expression"].mean().unstack().round(2)
    st.dataframe(mean)

    if "AD" in mean.columns and "Control" in mean.columns:
        fc = mean.copy()
        fc["Fold Change"] = fc["AD"] / fc["Control"]
        st.subheader("Fold Change Comparison")
        st.dataframe(fc.round(2))

# ==================== SINGLE GENE ANALYSIS ====================
if gene:
    gene_data = filtered[filtered["Gene"] == gene]

    if not gene_data.empty:
        description  = gene_data.iloc[0]["Gene Title"]
        control      = gene_data[gene_data["Condition"] == "Control"]["Expression"].astype(float)
        disease      = gene_data[gene_data["Condition"] == "AD"     ]["Expression"].astype(float)

        control_mean = control.mean() if len(control) else None
        disease_mean = disease.mean() if len(disease) else None

        if len(control) >= 2 and len(disease) >= 2:
            _, p_value = ttest_ind(control, disease)
        else:
            p_value = None

        if control_mean and control_mean != 0 and disease_mean is not None:
            fold_change = disease_mean / control_mean
        else:
            fold_change = None

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Gene Overview")
            st.markdown(f"## {gene}")
            st.markdown(description)

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

            # --- Interpretation ---
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

# ==================== GLOBAL INSIGHTS ====================
st.divider()
st.header("Global Insights")

# --- Top Genes ---
rank_mode = st.radio("Rank by:", ["Difference", "Fold Change"])

if not top_genes_df.empty and rank_mode in top_genes_df.columns:
    top_display = (
        top_genes_df
        .sort_values(rank_mode, ascending=False)
        .head(10)
    )
    st.subheader(f"Top 10 Genes by {rank_mode}")
    st.dataframe(top_display.round(2))

# --- Volcano Plot ---

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
        title="Volcano Plot (hover to explore genes)",
    )
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray",
                  annotation_text="p=0.05")
    fig.add_vline(x=1,  line_dash="dash", line_color="gray")
    fig.add_vline(x=-1, line_dash="dash", line_color="gray")

    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("Not enough data to render the volcano plot.")
