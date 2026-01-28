
🧬 GSE210616 — Spatial Transcriptomics of TNBC (NCBI GEO)

Title: Spatial transcriptomics of triple negative breast cancer
Accession: GSE210616 (NCBI GEO)
Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616

Summary:
Spatial transcriptomics profiling of Triple Negative Breast Cancer (TNBC) tumors using the 10x Genomics Visium platform. A total of 43 spatial sections from 22 patients — representing 28 tissue sections from 14 primary TNBC tumors — were profiled to capture gene expression with spatial localization. This dataset is designed to support analysis of spatial patterns in tumor heterogeneity and microenvironment among TNBC samples.  ￼

Key Details:
	•	Organism: Homo sapiens
	•	Technology: 10x Genomics Visium spatial transcriptomics (Illumina NovaSeq 6000)
	•	Samples: 43 sections across 14 patients
	•	Data Included: Processed spatial count matrices, high-resolution H&E images, spatial coordinates (spot barcodes), Loupe Browser .cloupe files, and scale factor JSONs.
	•	Clinical / Biological Context: Focused on racial disparities and TNBC architecture, with rich metadata for tumor regions and spatial gene expression patterns.  ￼

Suggested citation:
Bassiouni R, Carpten J, Idowu M, Craig D. Spatial transcriptomics of triple negative breast cancer. GEO Accession GSE210616 (2022).  ￼

⸻

📦 Zenodo: 14247036 — Breast Cancer Visium ST (Mixed Subtypes)

Link: https://zenodo.org/records/14247036

Summary:
This dataset compiles Visium spatial transcriptomics data across various breast cancer subtypes, including TNBC and other subtypes (luminal, HER2+). It was generated to analyze the tumor microenvironment and predict drug responses via functional enrichment and subpopulation detection. Includes processed Seurat objects, signature collections, and deconvolution references.  ￼

Key Details:
	•	Includes:
	•	ST Visium objects with SCTransform normalized counts and spot deconvolution results (predicted clonal composition)
	•	Processed single-cell RNA-seq reference objects
	•	Drug sensitivity and functional enrichment score objects
	•	Biological Focus: Tumor microenvironment interactions, intratumor heterogeneity, and potential drug response signatures.
	•	Platforms: 10x Visium

⸻

📊 Zenodo: 4739739 — Single-Cell + Spatial Breast Cancer Atlas

Link: https://zenodo.org/records/4739739

Summary:
Spatial transcriptomics and pathology-linked metadata from the Wu et al. study: a single-cell and spatially resolved atlas of human breast cancers. The dataset includes expression matrices, spatial coordinate data, H&E tissue images (annotated and raw), and clinical annotation for six primary breast tumors.  ￼

Key Details:
	•	Includes:
	•	Raw and filtered spatial count matrices
	•	Spatial image data
	•	Metadata with clinical subtype and pathological spot annotation
	•	Technology: 10x Genomics Visium
	•	Use Cases: spatially resolved profiling, spot-level clinical annotation, integration of spatial expression with pathology.

⸻

📍 Zenodo: 14204217 — ST TNBC (Comprehensive ST Output)

Link: https://zenodo.org/records/14204217

Summary:
Data for spatial transcriptomics focused specifically on TNBC, including arrays/subarrays, clinical metadata, spot-level classification/regression results, clustering outputs, deconvolved expression per annotation type, and high-resolution images. It provides R objects and annotated raw counts for exploratory analysis and classification of spatial spots.  ￼

Key Details:
	•	Content:
	•	rawCountsMatrices — raw count matrices
	•	Clinical data and slide metadata
	•	Clustering results and spatial megaclusters
	•	Deconvolution and annotation files
	•	Image files (H&E and associated annotation layers)

Use Cases: spatial spot classification, cluster analysis, subtype comparisons in TNBC.

⸻

📂 Zenodo: 3957257 — HER2-Positive Breast Tumor Spatial Deconvolution (not TNBC)

Link: https://zenodo.org/records/3957257

Summary:
Although not specifically TNBC, this dataset includes spatial transcriptomics (Visium) on HER2-positive breast tumors. It contains processed count matrices, H&E images, spot selection coordinates, and annotated meta data — valuable for comparative studies across breast cancer subtypes.  ￼

Key Details:
	•	Includes: processed count matrices, images, spot mapping files
	•	Focus: spatial deconvolution in HER2+ tumors — useful for benchmarking and cross-subtype analysis.

⸻

🔬 Dataset — Mendeley Data: gb83sywsjc/1 (Breast Ecosystem Atlas)

Link: https://data.mendeley.com/datasets/gb83sywsjc/1

Summary:
A breast tumor immune ecosystem dataset containing mass cytometry and immunofluorescence data from 144 tumor and 50 non-tumor samples. Not spatial transcriptomics per se, but valuable for integrative comparisons with spatial gene expression patterns.  ￼

Key Details:
	•	Type: Mass cytometry with immune markers, immunofluorescence
	•	Context: Characterizing tumor immune ecosystem and heterogeneity across breast cancers

⸻

🔍 cellxgene / Human Breast Cell Atlas (Multimodal Resource)

Links:
	•	https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637
	•	https://navinlabcode.github.io/HumanBreastCellAtlas.github.io/index.html#about

Summary:
A collection of spatial and single-cell expression data from the Human Breast Cell Atlas project. It includes multimodal datasets (scRNA-seq, ST, histology) with cell type annotations — useful for integrated analyses and referencing normal vs tumor tissue expression patterns.

Key Details:
	•	Platforms: Single-cell and spatial assays
	•	Utility: Atlas for cell type mapping and reference comparisons

⸻

🧠 HEST (Hugging Face) — Spatial Transcriptomics Image/Gene Dataset

Dataset: MahmoodLab/hest on HuggingFace
Link: https://huggingface.co/datasets/MahmoodLab/hest

Summary:
HEST (Histology and Expression ST dataset) aggregates >1,000 spatial transcriptomics profiles with linked histology images, enabling multi-modal learning (image ↔ expression). While not TNBC-specific, it contains cancer samples and is valuable for training models connecting morphology and gene expression.

Key Details:
	•	Content: 1,229 spatial profiles with image + gene expression
	•	Applications: multimodal representation learning, benchmarking, and model development
