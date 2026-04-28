# Gene\_Scope



Alzheimer's Disease Gene Expression Dashboard



GeneScope is an interactive web-based dashboard for exploring differential gene expression patterns in Alzheimer's Disease (AD) across multiple brain regions, cell types, and biological conditions.



The application provides an accessible interface for analyzing transcriptomic data, enabling both single-gene and multi-gene exploration with statistical evaluation and visualization.



\---



\## Application



(Add your Streamlit link here after deployment)  

https://your-app-name.streamlit.app



\---



\## Dataset Information



\- Source: Gene Expression Omnibus (GEO)  

\- Accession: GSE5281  

\- Organism: Homo sapiens  

\- Platform: Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)  

\- Experiment Type: Expression profiling by microarray  



\### Study Overview



This dataset investigates gene expression changes associated with Alzheimer's Disease and normal aging across six anatomically distinct brain regions:



\- Entorhinal Cortex  

\- Hippocampus  

\- Medial Temporal Gyrus  

\- Posterior Cingulate  

\- Superior Frontal Gyrus  

\- Primary Visual Cortex  



Samples were obtained from multiple Alzheimer's Disease Centers and processed using laser capture microdissection (LCM) to reduce tissue heterogeneity. Expression profiling was performed on approximately 55,000 transcripts.



\### Citation



Liang WS, Dunckley T, Beach TG, Grover A et al.  

Gene expression profiles in anatomically and functionally distinct regions of the normal aged human brain.  

Physiol Genomics. 2007;28(3):311–22  

PMID: 17077275  



\---



\## Features



\### Single Gene Analysis

\- Expression comparison between Control and AD

\- Mean expression values

\- Statistical testing using independent t-test

\- Fold change calculation (AD / Control)

\- Automated interpretation (upregulated, downregulated, or not significant)



\### Multi-Gene Comparison

\- Compare up to three genes simultaneously

\- Condition-wise expression distribution

\- Mean expression comparison table

\- Fold change comparison



\### Volcano Plot

\- Genome-wide differential expression visualization

\- Log2 fold change vs −log10(p-value)

\- Highlighting of statistically significant genes

\- Interactive hover information



\### Top Differential Genes Panel

\- Ranking based on:

&#x20; - Mean difference (AD − Control), or

&#x20; - Fold change (AD / Control)

\- Displays top regulated genes across the dataset



\### Dynamic Filtering

Filter the dataset by:

\- Brain Region  

\- Cell Type  

\- Condition (Control vs AD)  

\- Sex  



\### Data Export

\- Download gene-specific datasets for external analysis



\---



\## Statistical Framework



\- Statistical testing performed using independent two-sample t-test

\- Significance threshold: p < 0.05

\- Fold change interpretation:

&#x20; - >1 indicates higher expression in AD

&#x20; - <1 indicates lower expression in AD



\---



\## Technology Stack



\- Language: Python  

\- Framework: Streamlit  

\- Data Processing: Pandas, NumPy  

\- Statistical Analysis: SciPy  

\- Visualization: Matplotlib, Seaborn, Plotly  



\---



\## Project Structure



GeneScope/

├── app.py

├── requirements.txt

├── README.md

└── data/

&#x20;   └── processed/

&#x20;                 ├── final\_expression.csv

&#x20;                 └── metadata.csv





\---



\## Installation and Usage



\### 1. Clone the repository


https://github.com/pribug/Gene_Scope/tree/main


\### 2. Install dependencies



pip install -r requirements.txt





\### 3. Run the application



streamlit run app.py





\---



\## Deployment



This application can be deployed using Streamlit Cloud by connecting the GitHub repository and selecting `app.py` as the entry point. Deployment generates a public URL for sharing and access.



\---



\## Project Objective



To provide an interactive platform for exploring transcriptomic changes associated with Alzheimer's Disease, combining statistical analysis with visualization to support data-driven interpretation.



\---



\## Author



Priyam Tripathi  



\---



\## Notes



\- Initial loading time may be affected by dataset size and preprocessing steps  

\- All analyses are dynamically updated based on selected filters  

\- This tool is intended for exploratory analysis and research use  



\---

