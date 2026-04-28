import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np
import plotly.express as px
import re

sns.set_style("whitegrid")
st.set_page_config(page_title="GeneScope", layout="wide")

# ------------------ FONT ------------------
st.markdown("""
<style>
[data-testid="stMetricLabel"] { font-size: 12px !important; }
[data-testid="stMetricValue"] { font-size: 18px !important; }
</style>
""", unsafe_allow_html=True)

# ------------------ TITLE ------------------
st.title("GeneScope: Alzheimer's Gene Expression Dashboard")
st.caption("Explore differential gene expression across brain regions and conditions.")

# ------------------ LOAD DATA ------------------
@st.cache_data
def load_data():
    expression = pd.read_csv("data/processed/final_expression.csv", encoding="latin1")
    metadata = pd.read_csv("data/processed/metadata.csv", encoding="latin1")

    for col in metadata.columns:
        if metadata[col].dtype == "object":
            metadata[col] = metadata[col].str.replace("\xa0", "").str.strip()

    metadata["Condition"] = metadata["Condition"].str.lower().str.strip().str.replace("’", "'")
    metadata["Condition"] = metadata["Condition"].replace({
        "alzheimer's disease": "AD",
        "alzheimers disease": "AD",
        "ad": "AD",
        "control": "Control"
    })

    metadata["Sex"] = metadata["Sex"].str.lower().str.strip()
    metadata["Sex"] = metadata["Sex"].replace({
        "male": "Male",
        "female": "Female"
    })

    def is_fake_date(g):
        return bool(re.match(r"^\d{1,2}-[A-Za-z]{3}$", str(g)))

    expression = expression[~expression["Gene"].apply(is_fake_date)]
    expression["Gene"] = expression["Gene"].str.upper()

    expression_long = expression.melt(
        id_vars=["Gene", "Gene Title"],
        var_name="Sample",
        value_name="Expression"
    )

    return expression_long.merge(metadata, on="Sample")


data = load_data()
global_data = data.copy()

# ------------------ SUMMARY ------------------
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total Genes", len(data["Gene"].unique()))
col2.metric("Total Samples", len(data["Sample"].unique()))
col3.metric("Conditions", "Control vs AD")
col4.metric("Brain Regions", data["Brain_region"].nunique())

st.divider()

# ================== SEARCH ==================
st.header("Gene Search")

search = st.text_input("Search gene")

genes = sorted(data["Gene"].unique())
if search:
    genes = [g for g in genes if search.upper() in g]

selected = st.multiselect("Select up to 3 genes", genes)
gene = selected[0] if len(selected) == 1 else None

st.divider()

# ------------------ FILTERS ------------------
st.sidebar.header("Filters")

region = st.sidebar.multiselect("Brain Region", sorted(data["Brain_region"].dropna().unique()))
cell_type = st.sidebar.multiselect("Cell Type", sorted(data["Cell_type"].dropna().unique()))
condition = st.sidebar.multiselect("Condition", sorted(data["Condition"].dropna().unique()))
sex = st.sidebar.multiselect("Sex", sorted(data["Sex"].dropna().unique()))

filtered = data.copy()

if region:
    filtered = filtered[filtered["Brain_region"].isin(region)]
if cell_type:
    filtered = filtered[filtered["Cell_type"].isin(cell_type)]
if condition:
    filtered = filtered[filtered["Condition"].isin(condition)]
if sex:
    filtered = filtered[filtered["Sex"].isin(sex)]

# ================== RESULTS ==================

# MULTI GENE
if len(selected) > 1:
    st.subheader("Multi-Gene Comparison")

    multi = filtered[filtered["Gene"].isin(selected)]

    fig, ax = plt.subplots()
    sns.boxplot(x="Gene", y="Expression", hue="Condition", data=multi, ax=ax)
    st.pyplot(fig)

    mean = multi.groupby(["Gene", "Condition"])["Expression"].mean().unstack().round(2)
    st.dataframe(mean)

    if "AD" in mean.columns and "Control" in mean.columns:
        fc = mean.copy()
        fc["Fold Change"] = fc["AD"] / fc["Control"]
        st.subheader("Fold Change Comparison")
        st.dataframe(fc.round(2))

# SINGLE GENE
if gene:

    gene_data = filtered[filtered["Gene"] == gene]

    if not gene_data.empty:

        description = gene_data.iloc[0]["Gene Title"]

        control = gene_data[gene_data["Condition"] == "Control"]["Expression"].astype(float)
        disease = gene_data[gene_data["Condition"] == "AD"]["Expression"].astype(float)

        control_mean = control.mean() if len(control) else None
        disease_mean = disease.mean() if len(disease) else None

        if len(control) >= 2 and len(disease) >= 2:
            _, p_value = ttest_ind(control, disease)
        else:
            p_value = None

        if control_mean is not None and control_mean != 0 and disease_mean is not None:
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
            b.metric("AD Mean", round(disease_mean, 2) if disease_mean else "NA")
            c.metric("P-value", f"{p_value:.5f}" if p_value is not None else "NA")
            d.metric("Fold Change", round(fold_change, 2) if fold_change else "NA")

            st.markdown("---")

            st.markdown("**Data Quality Check:**")
            st.write(f"Control samples: {len(control)}")
            st.write(f"AD samples: {len(disease)}")

            st.markdown("---")

            st.markdown("**Context:**")
            st.write(f"Samples: {len(gene_data)}")
            st.write(f"Regions: {gene_data['Brain_region'].nunique()}")
            st.write(f"Cell Types: {gene_data['Cell_type'].nunique()}")

            # Interpretation
            st.subheader("Interpretation")

            if len(control) == 0 and len(disease) == 0:
                st.warning("No data after filtering")
            elif len(control) == 0:
                st.warning("Only AD samples available")
            elif len(disease) == 0:
                st.warning("Only Control samples available")
            elif len(control) < 2 or len(disease) < 2:
                st.warning("Not enough samples")
            elif p_value < 0.05:
                if disease_mean > control_mean:
                    st.success(f"Upregulated in AD (FC: {round(fold_change,2)})")
                else:
                    st.success(f"Downregulated in AD (FC: {round(fold_change,2)})")
            else:
                st.info("No significant difference")

        with col2:

            st.subheader("Expression Analysis")

            fig1, ax1 = plt.subplots()
            sns.boxplot(x="Condition", y="Expression", data=gene_data, ax=ax1)
            st.pyplot(fig1)

            fig2, ax2 = plt.subplots()
            sns.stripplot(x="Condition", y="Expression", data=gene_data, ax=ax2)
            st.pyplot(fig2)

# ================== GLOBAL INSIGHTS ==================
st.divider()
st.header("Global Insights")

# -------- Rank toggle --------
rank_mode = st.radio("Rank by:", ["Difference", "Fold Change"])

top = global_data.groupby(["Gene", "Condition"])["Expression"].mean().unstack()

if "AD" in top.columns and "Control" in top.columns:

    top = top.dropna()

    if rank_mode == "Difference":
        top["Score"] = top["AD"] - top["Control"]
    else:
        top["Score"] = top["AD"] / top["Control"]

    top = top.sort_values("Score", ascending=False).head(10)
    st.subheader(f"Top Genes ({rank_mode})")
    st.dataframe(top.round(2))

# -------- Volcano Plot --------
st.subheader("Volcano Plot")

volcano_list = []

for gene_name, df in global_data.groupby("Gene"):

    control = df[df["Condition"] == "Control"]["Expression"].astype(float)
    disease = df[df["Condition"] == "AD"]["Expression"].astype(float)

    if len(control) >= 2 and len(disease) >= 2:

        control_mean = control.mean()
        disease_mean = disease.mean()

        if control_mean != 0:

            fc = disease_mean / control_mean
            log2fc = np.log2(fc)

            _, pval = ttest_ind(control, disease)

            volcano_list.append({
                "Gene": gene_name,
                "log2FC": log2fc,
                "-log10p": -np.log10(pval),
                "Selected": "Yes" if gene_name in selected else "No"
            })

volcano = pd.DataFrame(volcano_list)

fig = px.scatter(
    volcano,
    x="log2FC",
    y="-log10p",
    color="Selected",
    hover_data=["Gene"],
    title="Volcano Plot (Click to explore genes)"
)

st.plotly_chart(fig, use_container_width=True)