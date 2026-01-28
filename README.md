# TNBC Spatial Transcriptomics Data Resources

A curated collection of publicly available spatial transcriptomics datasets for Triple-Negative Breast Cancer (TNBC) research.

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](CONTRIBUTING.md)

## Table of Contents

- [Overview](#overview)
- [TNBC-Specific Datasets](#tnbc-specific-datasets)
- [Mixed Breast Cancer Datasets (Including TNBC)](#mixed-breast-cancer-datasets-including-tnbc)
- [Reference Datasets](#reference-datasets)
- [Large-Scale Multi-Organ Collections](#large-scale-multi-organ-collections)
- [Data Access](#data-access)
- [Citation](#citation)
- [Contributing](#contributing)

## Overview

Triple-Negative Breast Cancer (TNBC) accounts for 10-15% of all breast cancers and is characterized by the lack of estrogen receptor (ER), progesterone receptor (PR), and HER2 expression. TNBC is highly heterogeneous and aggressive, with limited targeted treatment options. Spatial transcriptomics enables the study of gene expression patterns while preserving tissue architecture, providing crucial insights into the tumor microenvironment and cellular interactions.

This repository catalogs publicly available spatial transcriptomics datasets focused on or including TNBC samples, providing researchers with easy access to these valuable resources.


- Expression/morphology pair annotations
- Nuclei segmentation masks

---

## Data Access

### Comprehensive Dataset Table

| Dataset | Source | TNBC | Technology | Samples / Scale | Modalities | Key Contents | Size | Link |
|---------|--------|------|------------|----------------|------------|--------------|------|------|
| GSE210616 | GEO |  Yes | 10x Visium | 43 sections, 22 patients | ST + H&E | Raw & processed matrices, spatial coords, H&E images, Loupe files | 35.1 GB | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616) |
| Zenodo 14204217 | Zenodo | Yes | 10x Visium | Multiple TNBC samples | ST + H&E + IHC | Raw counts, clustering, deconvolution, clinical metadata, 18 annotation types | 58.0 GB | [Zenodo](https://zenodo.org/records/14204217) |
| Zenodo 14247036 | Zenodo |  Partial | 10x Visium | 9 tumors (4 TNBC) | ST + scRNA-seq | Seurat objects, deconvolution, drug response signatures (>1,200 drugs) | 7.0 GB | [Zenodo](https://zenodo.org/records/14247036) |
| Zenodo 4739739 | Zenodo |  Mixed | 10x Visium | 6 tumors (4 TNBC) | ST + H&E | Spatial matrices, annotated histology, pathologist annotations | 920 MB | [Zenodo](https://zenodo.org/records/4739739) |
| Zenodo 3957257 | Zenodo | HER2+ | 10x Visium | HER2+ tumors | ST + H&E | Processed counts, pathology images, spatial deconvolution | 629.6 MB | [Zenodo](https://zenodo.org/records/3957257) |
| Human Breast Cell Atlas | CellxGene |  Mixed | scRNA-seq + ST | 126 women, 714K cells | scRNA-seq + ST + CODEX + smFISH | Cell-type annotations, spatial maps, reference atlas | Variable | [CellxGene](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637) |
| Mendeley gb83sywsjc | Mendeley |  Mixed | Mass Cytometry | 194 samples (144 tumor) | CyTOF + IF | Immune ecosystem profiling, 73 proteins, 26M cells | Variable | [Mendeley](https://data.mendeley.com/datasets/gb83sywsjc/1) |
| HEST-1k | Hugging Face |  Mixed | Multiple ST | 1,255 slides, 367 cancer samples | ST + WSI | Histology-expression pairs, 76M nuclei, multi-organ | >100B | [HuggingFace](https://huggingface.co/datasets/MahmoodLab/hest) |

**Legend:**
- **Yes** - Dataset exclusively or primarily contains TNBC samples
- **Partial/Mixed** - Dataset contains some TNBC samples along with other subtypes
- **No** - Dataset does not contain TNBC, but useful as reference (HER2+, normal tissue, etc.)

### Data Types Available

- ✅ Raw count matrices
- ✅ Processed/normalized counts
- ✅ H&E whole slide images
- ✅ Pathological annotations
- ✅ Clinical metadata
- ✅ Spot coordinates/positions
- ✅ Scalefactors for alignment
- ✅ Analysis code/pipelines
- ✅ Deconvolution results
- ✅ Clustering information

---

## TNBC-Specific Datasets

### 1. USC TNBC Cohort (GSE210616)

**Study:** Spatial transcriptomics of triple negative breast cancer  
**Published:** 2022  
**Platform:** 10x Genomics Visium  
**Samples:** 43 tissue sections from 22 TNBC patients
- 15 African-American patients
- 7 Caucasian patients
- 2 sections per patient (except patient 19)

**Key Features:**
- Focus on racial disparities in TNBC
- Primary tumor samples
- Pathologically annotated

**Access:**
- **GEO:** [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616)
- **Publication:** [PMID: 36283023](https://pubmed.ncbi.nlm.nih.gov/36283023/)
- **Data Type:** Raw count matrices, H&E images (JPG, PNG), metadata
- **Size:** 35.1 GB (RAW.tar)

**Note:** Raw sequencing data not submitted due to patient privacy concerns. Processed data available as supplementary files.

---

### 2. Belgian TNBC Spatial Atlas (Zenodo: 14204217)

**Study:** Spatial transcriptomics reveals substantial heterogeneity in triple-negative breast cancer with potential clinical implications  
**Published:** Nature Communications, 2024  
**Platform:** Custom spatial transcriptomics  
**Institution:** Université Libre de Bruxelles  
**Samples:** Multiple TNBC patients with comprehensive annotations

**Key Features:**
- Extensive pathological annotations (18 categories)
- Detailed spatial clustering analysis
- Megaclustering across patients
- Spot-level classification/regression
- TLS (Tertiary Lymphoid Structures) identification
- Integration with bulk RNA-seq and immunotherapy datasets
- CD3/CD20 IHC images for validation

**Access:**
- **Zenodo:** [10.5281/zenodo.14204217](https://zenodo.org/records/14204217)
- **GitHub:** [BCTL-Bordet/ST](https://github.com/BCTL-Bordet/ST)
- **Data Type:** Count matrices, H&E images, annotations, clusterings, deconvolution results, projections
- **Size:** 58.0 GB

**Data Organization:**
```
├── byArray/           # Data by array/subarray with spot coordinates
├── classification/    # Spot-level classification and regression models
├── Clinical/          # Clinical data and sample metadata
├── clustering/        # Intra-patient clustering and megaclustering
├── deconvolution/     # Per-annotation gene expression profiles
├── external datasets/ # Validation datasets (SCAN-B, METABRIC, I-SPY2)
├── Images/            # Original H&E and annotated images (NDPI, JPG, PNG)
├── patches/           # Patch size distribution data
├── rawCountsMatrices/ # Raw count matrices (TSV)
├── Robjects/          # R objects for analysis
└── projections/       # Spatial visualizations and mappings
```

**Annotation Categories:**
- Tumor, In situ carcinoma, Necrosis
- Stroma (Low TIL, High TIL, Acellular)
- Fat tissue, Vessels, Lymphoid nodules
- Lactiferous ducts, Nerves, Heterologous elements
- Artefacts, Holes (whitespace)

---

## Mixed Breast Cancer Datasets (Including TNBC)

### 3. Spanish CNIO Drug Response Study (Zenodo: 14247036)

**Study:** Spatial Transcriptomics in Breast Cancer Reveals Tumour Microenvironment-Driven Drug Responses and Clonal Therapeutic Heterogeneity  
**Published:** 2024  
**Platform:** 10x Genomics Visium  
**Institution:** Spanish National Cancer Research Centre (CNIO)  
**Samples:** 9 invasive breast adenocarcinomas
- 2 Luminal (CID4290, CID4535)
- 4 TNBC (CID44971, CID4465, 1142243F, 1160920F)
- 3 HER2+ (738811QB, 1168993F, V19L29)

**Key Features:**
- Drug sensitivity predictions for >1,200 drugs using Beyondcell
- Functional pathway enrichment scores
- Spot deconvolution with SCEVAN clonal composition
- SCTransform-normalized counts
- Matched single-cell RNA-seq reference data

**Access:**
- **Zenodo:** [10.5281/zenodo.14247036](https://zenodo.org/records/14247036)
- **Code Repositories:**
  - Drug signatures: [cnio-bu/SSc-breast](https://github.com/cnio-bu/SSc-breast)
  - Preprocessing: [cnio-bu/ST-preprocess](https://github.com/cnio-bu/ST-preprocess)
  - Analysis: [cnio-bu/breast-bcspatial](https://github.com/cnio-bu/breast-bcspatial)
- **Data Type:** Seurat objects, Beyondcell objects, gene signatures
- **Size:** 7.0 GB

**Data Structure:**
```
├── signatures/
│   ├── SSc breast/          # Drug response signatures (>1,200 drugs)
│   └── Functional signatures/ # Pathway enrichment signatures
├── visium/                  # Processed ST Seurat objects
├── single-cell/             # scRNA-seq reference objects
└── beyondcell/
    ├── sensitivity/         # Drug sensitivity predictions
    └── functional/          # Pathway enrichment scores
```

---

### 4. Wu et al. Breast Cancer Atlas (Zenodo: 4739739)

**Study:** A single-cell and spatially resolved atlas of human breast cancers  
**Published:** Nature Genetics, 2021  
**Platform:** 10x Genomics Visium  
**Institution:** Garvan Institute of Medical Research  
**Samples:** 6 primary breast cancers
- 2 ER+ (CID4535, CID4290)
- 4 TNBC (CID44971, CID4465, 1142243F, 1160920F)

**Key Features:**
- Pathologist annotations for each spot
- Clinical subtype information
- Matched with extensive scRNA-seq data (SCP1039)
- Stereoscope deconvolution
- Identification of 9 tumor ecotypes
- Spatial organization of stromal-immune niches

**Access:**
- **Zenodo:** [10.5281/zenodo.4739739](https://zenodo.org/records/4739739)
- **Single Cell Portal:** [SCP1039](https://singlecell.broadinstitute.org/single_cell/study/SCP1039/)
- **Publication:** [Wu et al., Nat Genet 2021](https://www.nature.com/articles/s41588-021-00911-1)
- **Data Type:** SpaceRanger outputs, filtered matrices, metadata, H&E images
- **Size:** 920.8 MB

**Data Contents:**
```
├── raw_count_matrices.tar.gz       # SpaceRanger raw outputs
├── filtered_count_matrices.tar.gz  # Filtered count matrices
├── spatial.tar.gz                  # Images, scalefactors, positions
├── metadata.tar.gz                 # Clinical and pathological annotations
└── images.pdf                      # H&E and annotation visualization
```

**Key Findings:**
- Nine ecotypes associated with cellular heterogeneity and prognosis
- PD-L1/PD-L2+ macrophage populations linked to clinical outcomes
- Spatial organization of tumor-immune interactions

---

### 5. Andersson et al. HER2+ Dataset (Zenodo: 3957257)

**Study:** Spatial Deconvolution of HER2-positive Breast Tumors Reveals Novel Intercellular Relationships  
**Published:** 2020  
**Platform:** 10x Visium  
**Institution:** Royal Institute of Technology (KTH) and Science for Life Laboratory  
**Samples:** HER2+ breast tumors (NOT TNBC)

**Relevance to TNBC Research:**
- **Comparative reference** for HER2+ vs TNBC spatial differences
- Methodology applicable to TNBC analysis
- Spatial deconvolution approach

**Key Features:**
- Processed count matrices
- H&E images (plain and annotated)
- Pathologist annotations
- Spot selection files for visualization

**Access:**
- **Zenodo:** [10.5281/zenodo.3957257](https://zenodo.org/records/3957257)
- **Publication:** [DOI: 10.1101/2020.07.14.200600](https://doi.org/10.1101/2020.07.14.200600)
- **Data Type:** Count matrices, H&E images, annotations, spot coordinates
- **Size:** 629.6 MB

**Data Contents:**
```
├── count-matrices.zip      # Processed counts [spots × genes]
├── images.zip             # H&E and annotated images
├── spot-selection.zip     # Array to pixel coordinate mapping
└── meta.zip              # Spot labels and annotations
```

**Annotation Categories:**
- Breast glands
- Connective tissue
- Immune infiltrates
- Tumor regions

**Note:** This dataset is HER2+, not TNBC. Included as a valuable reference for comparative studies and methodological approaches applicable to TNBC.

---

### 6. Wagner et al. Mass Cytometry Atlas (Mendeley: gb83sywsjc)

**Study:** A single-cell atlas of the tumor and immune ecosystem of human breast cancer  
**Published:** Cell, 2019  
**Platform:** Mass Cytometry (CyTOF)  
**Institution:** University of Zurich  
**Samples:** 194 samples (144 tumor, 50 non-tumor)
- Includes ER+, HER2+, and TNBC samples

**Relevance to TNBC Research:**
- **Immune ecosystem profiling** at single-cell protein level
- Complements spatial transcriptomics with protein data
- Identifies immunosuppressive features in TNBC

**Key Features:**
- 73 proteins measured per cell
- 26 million cells analyzed
- Tumor and immune cell-centric panels
- PD-L1+ macrophage characterization
- Exhausted T cell populations

**Access:**
- **Mendeley Data:** [10.17632/gb83sywsjc.1](https://data.mendeley.com/datasets/gb83sywsjc/1)
- **Publication:** [Wagner et al., Cell 2019](https://doi.org/10.1016/j.cell.2019.01.003)
- **License:** CC BY 4.0
- **Data Type:** CyTOF data, flow cytometry analysis

**Key Findings (relevant to TNBC):**
- High frequencies of PD-L1+ tumor-associated macrophages in high-grade ER+ and ER- tumors
- Exhausted T cell populations in TNBC
- Ecosystem-based classification for precision medicine
- Tumor cell phenotypic heterogeneity

**Technologies:**
- Mass cytometry (CyTOF) for protein quantification
- Imaging mass cytometry (IMC) for spatial protein distribution

**Note:** While not spatial transcriptomics, this dataset provides essential protein-level immune profiling that complements ST data. Particularly valuable for validating immune populations identified in TNBC spatial studies.

---

## Reference Datasets

### 7. Human Breast Cell Atlas (CZI/CellxGene)

**Study:** A spatially resolved single cell genomic atlas of the adult human breast  
**Published:** Nature, 2024  
**Platform:** Multiple (Visium, CODEX, Resolve smFISH, MERSCOPE)  
**Institution:** Multi-institutional (Navin Lab, MD Anderson; Kessenbrock Lab, UC Irvine; Lawson Lab, UC Irvine)  
**Samples:** 126 women (normal breast tissue, not cancer)

**Key Features:**
- **Normal breast tissue reference** (not tumor)
- 714,331 cells profiled by scRNA-seq
- 117,346 nuclei profiled by snRNA-seq
- 12 major cell types, 58 biological cell states
- Diverse biological states (parity, menopause, age, ethnicity)
- Multi-modal spatial profiling

**Relevance to TNBC Research:**
- **Reference for normal tissue comparison**
- Cell type identification and marker discovery
- Spatial organization patterns in healthy breast

**Access:**
- **Website:** [HumanBreastCellAtlas.github.io](https://navinlabcode.github.io/HumanBreastCellAtlas.github.io/)
- **CellxGene:** [Collection 48259aa8](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637)
- **GitHub:** [navinlabcode/HumanBreastCellAtlas](https://github.com/navinlabcode/HumanBreastCellAtlas)
- **Protocols:** Tissue dissociation and analysis code available

**Technologies Used:**
- 10X Genomics Chromium 3' (scRNA-seq, snRNA-seq)
- 10X Genomics Visium (Spatial Transcriptomics)
- Akoya CODEX (34-antibody panel)
- Resolve Biosciences smFISH (100-gene panel)
- Vizgen MERSCOPE (140-300 gene panels)

---

## Large-Scale Multi-Organ Collections

### 8. HEST-1k (Hugging Face)

**Study:** Histology-based Expression Spatial Transcriptomics  
**Published:** 2024  
**Platform:** Multiple spatial transcriptomics platforms  
**Institution:** Mahmood Lab (Harvard Medical School)  
**Samples:** 1,255 spatial transcriptomic profiles

**Coverage:**
- **26 organs** (2 species: Homo Sapiens, Mus Musculus)
- **367 cancer samples** from 25 cancer types
- **Breast cancer samples included** (TNBC among them)
- 1.5M expression/morphology pairs
- 76M nuclei identified

**Key Features:**
- Whole Slide Images (WSI) with pixel size < 1.15 µm/px
- Aligned H&E images with spatial transcriptomics
- Comprehensive metadata
- Multi-organ cancer atlas
- Foundation model training resource

**Access:**
- **Hugging Face:** [MahmoodLab/hest](https://huggingface.co/datasets/MahmoodLab/hest)
- **Paper:** [arXiv:2406.16192](https://arxiv.org/abs/2406.16192)
- **License:** CC BY-NC-SA 4.0
- **Status:** Gated dataset (request access)
- **Size:** >100B parameters

**Data Organization:**
- ST profiles from 131 public and internal cohorts
- Unified format with aligned WSI


### Download Instructions

#### GEO Datasets (GSE210616)
```bash
# Using wget
wget -r -np -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210616/

# Or use GEOquery in R
library(GEOquery)
gse <- getGEO("GSE210616")
```

#### Zenodo Datasets
```bash
# Direct download using wget
wget https://zenodo.org/records/14204217/files-archive  # Belgian TNBC
wget https://zenodo.org/records/14247036/files-archive  # CNIO Drug Response
wget https://zenodo.org/records/4739739/files-archive   # Wu et al.
wget https://zenodo.org/records/3957257/files-archive   # HER2+ Reference

# Or use zenodo_get Python package
pip install zenodo_get
zenodo_get 10.5281/zenodo.14204217
```

#### Mendeley Dataset (Mass Cytometry)
```bash
# Mendeley requires browser download
# 1. Visit: https://data.mendeley.com/datasets/gb83sywsjc/1
# 2. Click "Download All"
# 3. Extract downloaded files
```

#### HEST Dataset (Hugging Face)
```python
from datasets import load_dataset

# Request access first at: https://huggingface.co/datasets/MahmoodLab/hest
dataset = load_dataset("MahmoodLab/hest")
```

#### CellxGene/HBCA
Visit the [CellxGene portal](https://cellxgene.cziscience.com/) and download through the web interface or use the API.

---

## Citation

If you use these datasets in your research, please cite the original publications:

### GSE210616
```bibtex
@article{bassiouni2022spatial,
  title={Spatial transcriptomics of triple negative breast cancer},
  author={Bassiouni, R and Carpten, J and Idowu, M and Craig, D},
  journal={GEO Accession GSE210616},
  year={2022},
  publisher={NCBI}
}
```

### Belgian TNBC Atlas (Zenodo 14204217)
```bibtex
@article{venet2024spatial,
  title={Spatial transcriptomics reveals substantial heterogeneity in triple-negative breast cancer with potential clinical implications},
  author={Venet, David and others},
  journal={Nature Communications},
  year={2024},
  doi={10.5281/zenodo.14204217}
}
```

### CNIO Drug Response Study (Zenodo 14247036)
```bibtex
@dataset{jimenez2024spatial,
  title={Spatial Transcriptomics in Breast Cancer Reveals Tumour Microenvironment-Driven Drug Responses and Clonal Therapeutic Heterogeneity},
  author={Jiménez-Santos, María José and García-Martín, Santiago and others},
  year={2024},
  doi={10.5281/zenodo.14247036}
}
```

### Wu et al. Nature Genetics (Zenodo 4739739)
```bibtex
@article{wu2021single,
  title={A single-cell and spatially resolved atlas of human breast cancers},
  author={Wu, Sunny Z and Al-Eryani, Ghamdan and Roden, Daniel L and others},
  journal={Nature Genetics},
  volume={53},
  number={9},
  pages={1334--1347},
  year={2021},
  doi={10.1038/s41588-021-00911-1}
}
```

### Andersson et al. HER2+ (Zenodo 3957257)
```bibtex
@article{andersson2021spatial,
  title={Spatial deconvolution of HER2-positive breast tumors reveals novel intercellular relationships},
  author={Andersson, Alma and Larsson, Ludvig and Stenbeck, Linnea and others},
  journal={Genome Biology},
  volume={22},
  number={1},
  pages={1--27},
  year={2021},
  doi={10.1186/s13059-021-02271-1}
}
```

### Wagner et al. Mass Cytometry (Mendeley gb83sywsjc)
```bibtex
@article{wagner2019single,
  title={A single-cell atlas of the tumor and immune ecosystem of human breast cancer},
  author={Wagner, Johanna and Rapsomaniki, Maria A and Chevrier, Stéphane and others},
  journal={Cell},
  volume={177},
  number={5},
  pages={1330--1345},
  year={2019},
  doi={10.1016/j.cell.2019.01.003}
}
```

### Human Breast Cell Atlas
```bibtex
@article{kumar2024spatially,
  title={A spatially resolved single cell genomic atlas of the adult human breast},
  author={Kumar, Tapsi and others},
  journal={Nature},
  year={2024},
  note={In press}
}
```

### HEST-1k
```bibtex
@article{jaume2024hest,
  title={HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis},
  author={Jaume, Guillaume and others},
  journal={arXiv preprint arXiv:2406.16192},
  year={2024}
}
```

---

## Contributing

We welcome contributions to expand this resource! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

To suggest a new dataset:
1. Open an issue with the dataset details
2. Include: publication link, data access URL, sample count, platform used
3. Provide a brief description of what makes it valuable for TNBC research

---

## Related Resources

### Analysis Tools
- **Seurat:** R toolkit for single-cell and spatial analysis
- **Scanpy/Squidpy:** Python tools for spatial transcriptomics
- **Stereoscope:** Probabilistic deconvolution
- **Beyondcell:** Drug response prediction
- **SCEVAN:** Clonal composition analysis
- **Loupe Browser:** 10x Genomics visualization tool

### Other Breast Cancer Spatial Datasets
- **HER2+ datasets:** Zenodo 3957257 (Andersson et al.)
- **ER+ samples:** Included in Wu et al. and CNIO studies
- **Breast cancer progression:** Various GEO datasets

### Complementary Data Types
- **scRNA-seq:** SCP1039, Human Cell Atlas
- **Mass Cytometry:** Mendeley gb83sywsjc (Wagner et al.)
- **Multiplex Imaging:** CODEX, MERSCOPE data from HBCA

---

## Acknowledgments

This repository aggregates information about publicly available datasets. We thank all the researchers, institutions, and patients who contributed to making this data available for the scientific community.

**Key Contributing Institutions:**
- University of Southern California (USC)
- Université Libre de Bruxelles (ULB)
- Spanish National Cancer Research Centre (CNIO)
- Garvan Institute of Medical Research
- MD Anderson Cancer Center
- University of California, Irvine
- Harvard Medical School / Mahmood Lab

---

## License

This repository documentation is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Individual datasets have their own licenses:
- Most spatial transcriptomics data: CC BY 4.0
- HEST-1k: CC BY-NC-SA 4.0 (Non-Commercial)
- Check individual dataset licenses before use

---

## Contact

For questions or suggestions about this repository:
- Open an issue on GitHub
- Contact: [Your contact information]

**Last Updated:** January 2026

---

**Keywords:** Triple-Negative Breast Cancer, TNBC, Spatial Transcriptomics, Visium, Tumor Microenvironment, Cancer Genomics, Single-Cell, Breast Cancer Atlas, H&E Imaging, Spatial Omics
