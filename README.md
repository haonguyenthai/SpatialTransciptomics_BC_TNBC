# SpatialTransciptomics_BC_TNBC
Visium data for TNBC 
🧾 How to Use This Repo
	1.	Clone the repo:
  git clone https://github.com/<your-username>/tnbc-spatial-transcriptomics-data.git
	2.	Explore the datasets/ folder to find dataset descriptions, links, and metadata.
	3.	Contribute by adding:
	•	more detailed summaries
	•	sample metadata tables
	•	scripts for automated downloading & preprocessing

Nice, let’s do it properly 👍
Below is a clean, GitHub-ready Markdown table summarizing all TNBC-relevant spatial transcriptomics datasets you listed. You can paste this directly into:
	•	README.md (main summary table), or
	•	datasets/README.md (catalog page)

⸻

📊 TNBC Spatial Transcriptomics Dataset Summary

Dataset	Source	TNBC	Technology	Samples / Scale	Modalities	Key Contents	Link
GSE210616	GEO	✅ Yes	10x Visium	43 sections, 22 patients (14 TNBC tumors)	ST + H&E	Raw & processed count matrices, spatial coordinates, H&E images, Loupe files, metadata	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616
Zenodo 14247036	Zenodo	⚠️ Partial	10x Visium	Multiple breast tumors	ST + scRNA-seq	Seurat objects, deconvolution results, drug response signatures, TME analysis	https://zenodo.org/records/14247036
Zenodo 14204217	Zenodo	✅ Yes	10x Visium	Multiple TNBC samples	ST + H&E	Raw count matrices, clustering results, deconvolution, clinical metadata, images	https://zenodo.org/records/14204217
Zenodo 4739739	Zenodo	⚠️ Mixed	10x Visium	6 breast tumors	ST + H&E	Spatial count matrices, annotated histology, clinical labels	https://zenodo.org/records/4739739
Zenodo 3957257	Zenodo	❌ No (HER2+)	10x Visium	HER2+ tumors	ST + H&E	Processed counts, spot coordinates, pathology images	https://zenodo.org/records/3957257
Human Breast Cell Atlas	cellxgene / Navin Lab	⚠️ Mixed	scRNA-seq + ST	Atlas-scale	scRNA-seq + ST	Cell-type annotations, spatial & single-cell references	https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637
Mendeley gb83sywsjc	Mendeley Data	⚠️ Mixed	Mass cytometry	194 samples	CyTOF + IF	Immune ecosystem profiling, tumor vs non-tumor	https://data.mendeley.com/datasets/gb83sywsjc/1
HEST	HuggingFace	⚠️ Mixed	Multiple ST platforms	1,229 slides	ST + H&E	Large-scale histology–expression pairs for multimodal learning	https://huggingface.co/datasets/MahmoodLab/hest


