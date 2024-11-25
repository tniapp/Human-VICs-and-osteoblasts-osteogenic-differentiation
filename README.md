<h1>Human-VICs-and-osteoblasts-osteogenic-differentiation</h1>
This code replicates the bioinformatic analysis conducted on proteomic and transcriptomic data obtained from human valve interstitial cells (VIC) and osteoblasts. The source of this analysis stems from the article titled "Similar, but not the same: multi-omics comparison of human valve interstitial cells and osteoblast osteogenic differentiation expanded with an estimation of data-dependent and data-independent PASEF proteomics".
In this study, we employed the timsTOF Pro platform to explore the proteomic profiles of valve interstitial cells (VICs) and osteoblasts during osteogenic differentiation, utilizing three data acquisition/analysis techniques: Data-Dependent Acquisition (DDA-PASEF) and Data-Independent Acquisition (DIA-PASEF) with a classic library based and machine learning-based “library-free” search (DIA-ML). RNA-seq complemented comparative proteome coverage analysis to provide a comprehensive biological reference.

<h2>Repository Structure</h2>

<pre>
Human-VICs-and-osteoblasts-osteogenic-differentiation-main/
├── LICENSE
├── README.md
├── DDA/
│   ├── dda_final.R
│   ├── dda_imp.csv
│   ├── dda_ost_contVSvic_cont.csv
│   ├── dda_ost_diffVSvic_diff.csv
│   └── dda_vic_contVSvic_diff.csv
├── DIA-ML/
│   ├── diaml_final.R
│   ├── diaml_imp.csv
│   ├── diaml_ost_contVSvic_cont.csv
│   ├── diaml_ost_diffVSvic_diff.csv
│   └── diaml_vic_contVSvic_diff.csv
├── DIA/
│   ├── dia_final.R
│   ├── dia_imp.csv
│   ├── dia_ost_contVSvic_cont.csv
│   ├── dia_ost_diffVSvic_diff.csv
│   └── dia_vic_contVSvic_diff.csv
├── data_prep/
│   ├── combined_peptide_DDA.tsv
│   ├── diaml_imp.csv
│   ├── report..gg_matrix_diaml.tsv
│   ├── data_preparation.r
│   ├── diaml_toNA.csv
│   ├── report..pr_matrix_DiA_ph.tsv
│   ├── dda_imp.csv
│   ├── fact_NA.xlsx
│   ├── report..pr_matrix_diaml.tsv
│   ├── dda_toNA.csv
│   ├── prot_dda.xlsx
│   ├── transcriptomics_data.csv
│   ├── dia_imp.csv
│   ├── rename.xlsx
│   ├── dia_toNA.csv
│   └── report..gg_matrix_DiA_ph.tsv
├── omics_comb/
│   ├── Enrichment_contr.R
│   ├── Enrichment_dif.R
│   ├── dia_ost_contVSost_diff.csv
│   ├── dia_ost_contVSvic_cont.csv
│   ├── dia_ost_diffVSvic_diff.csv
│   ├── dia_vic_contVSvic_diff.csv
│   └── deg_combination/
│       ├── c_vs_d/
│       │   └── VIC_OB_dif_deg.r
│       └── vic_vs_ob/
│           └── VIC_OB_c_comb_deg.r
├── qualitative_analysis/
│   └── DIA_qualit_analysis.r
└── transcriptome/
    ├── vic_ost_counts.csv
    ├── vic_ost_metadata.xlsx
    └── vic_ost_rna_analysis.R
</pre>

<h2>Repository Content Description</h2>

<h3>Folder <code>DDA</code></h3>
<ul>
    <li><strong>dda_final.R</strong> — R script for processing DDA data.</li>
    <li><strong>dda_imp.csv</strong> — File with imputed data results.</li>
    <li><strong>dda_ost_contVSvic_cont.csv</strong>, <strong>dda_ost_diffVSvic_diff.csv</strong>, <strong>dda_vic_contVSvic_diff.csv</strong> — Comparative analysis results for various DDA data conditions.</li>
</ul>

<h3>Folder <code>DIA-ML</code></h3>
<ul>
    <li><strong>diaml_final.R</strong> — R script for DIA-ML data analysis.</li>
    <li><strong>diaml_imp.csv</strong> — Imputed data for DIA-ML.</li>
    <li><strong>diaml_ost_contVSvic_cont.csv</strong>, <strong>diaml_ost_diffVSvic_diff.csv</strong>, <strong>diaml_vic_contVSvic_diff.csv</strong> — Comparative analysis results for DIA-ML data.</li>
</ul>

<h3>Folder <code>DIA</code></h3>
<ul>
    <li><strong>dia_final.R</strong> — Main R script for DIA data processing.</li>
    <li><strong>dia_imp.csv</strong> — Imputed data.</li>
    <li><strong>dia_ost_contVSvic_cont.csv</strong>, <strong>dia_ost_diffVSvic_diff.csv</strong>, <strong>dia_vic_contVSvic_diff.csv</strong> — Analysis results for various DIA data conditions.</li>
</ul>

<h3>Folder <code>data_prep</code></h3>
<ul>
    <li>Raw proteome data.</li>
</ul>

<h3>Folder <code>omics_comb</code></h3>
<ul>
    <li><strong>Enrichment_contr.R</strong> and <strong>Enrichment_dif.R</strong> — R scripts for enrichment analysis of control and differentiated cells.</li>
    <li><strong>deg_combination/</strong> — Contains two subdirectories for enrichment analysis of diffexpressed genes (DEG) in various conditions (control against differentiated cells; VIC against OB).</li>
</ul>

<h3>Folder <code>qualitative_analysis</code></h3>
<ul>
    <li><strong>DIA_qualit_analysis.r</strong> — R script for qualitative analysis of DIA data.</li>
</ul>

<h3>Folder <code>transcriptome</code></h3>
<ul>
    <li><strong>vic_ost_counts.csv</strong> — Count matrix for transcriptomic analysis.</li>
    <li><strong>vic_ost_metadata.xlsx</strong> — Metadata for VIC and osteoblast (OB) samples.</li>
    <li><strong>vic_ost_rna_analysis.R</strong> — Script for transcriptome analysis.</li>
</ul>

<h2>Sample Information</h2>

<p>The analyzed samples include:</p>
<ul>
    <li>Valve Interstitial Cells (VIC) and osteoblast (OB) cell samples.</li>
    <li>Divided into two groups: control and differentiated samples for each cell type.</li>
    <li>Raw reads are deposited in the NCBI SRA database with BioProject identifier PRJNA947173.</li>
</ul>
