# Dataset Comparison Guide

Detailed comparison of TNBC spatial transcriptomics datasets to help researchers choose the right data for their studies.

## Overview Table

| Feature | GSE210616 | Zenodo 14204217 | Zenodo 14247036 | Zenodo 4739739 | HEST-1k |
|---------|-----------|-----------------|-----------------|----------------|---------|
| **TNBC Samples** | 22 patients | Multiple | 4 samples | 4 samples | Subset available |
| **Total Samples** | 43 sections | Multiple | 9 | 6 | 367 cancer samples |
| **Platform** | Visium | Custom ST | Visium | Visium | Multi-platform |
| **Year** | 2022 | 2024 | 2024 | 2021 | 2024 |
| **Institution** | USC | ULB Brussels | CNIO Spain | Garvan AU | Harvard |
| **Size** | 35 GB | 58 GB | 7 GB | 920 MB | >100B |
| **Raw Data** | ❌ | ✅ | ✅ | ✅ | ✅ |
| **H&E Images** | ✅ | ✅ | ✅ | ✅ | ✅ (WSI) |
| **Annotations** | Basic | Extensive (18 types) | Moderate | Pathologist | Variable |
| **Clinical Data** | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Code Available** | ❌ | ✅ | ✅ | ❌ | ✅ |
| **scRNA-seq Match** | ❌ | ❌ | ✅ | ✅ (SCP1039) | Variable |

## Detailed Comparisons

### 1. Sample Characteristics

#### GSE210616 (USC TNBC Cohort)
- **Unique Feature:** Focus on racial disparities (African-American vs Caucasian)
- **Clinical Annotations:** Limited
- **Best For:** Studying racial differences in TNBC
- **Limitations:** No raw sequencing data, basic annotations

#### Zenodo 14204217 (Belgian TNBC Atlas)
- **Unique Feature:** Most comprehensive annotations (18 tissue categories)
- **Clinical Annotations:** Extensive survival data, treatment information
- **Best For:** Detailed spatial analysis, TLS studies, heterogeneity
- **Includes:**
  - Clusterings at multiple resolutions
  - Megaclusters across patients
  - Deconvolution results
  - External validation datasets (SCAN-B, METABRIC, I-SPY2)
  - IHC images (CD3/CD20)
- **Limitations:** Custom platform (not standard Visium)

#### Zenodo 14247036 (CNIO Drug Response)
- **Unique Feature:** Drug sensitivity predictions (>1,200 drugs)
- **Clinical Annotations:** Moderate
- **Best For:** Drug response studies, TME analysis, therapeutic heterogeneity
- **Includes:**
  - Beyondcell drug sensitivity scores
  - Functional pathway enrichment
  - Matched scRNA-seq data
  - Clonal composition (SCEVAN)
- **Limitations:** Only 4 TNBC samples (mixed with ER+ and HER2+)

#### Zenodo 4739739 (Wu et al. Nature Genetics)
- **Unique Feature:** Integration with comprehensive scRNA-seq atlas
- **Clinical Annotations:** Good
- **Best For:** Single-cell spatial integration, ecotype analysis
- **Includes:**
  - Pathologist annotations per spot
  - 9 ecotypes identified
  - Extensive scRNA-seq reference (SCP1039)
  - Immune landscape characterization
- **Limitations:** Smaller dataset (6 total samples)

#### HEST-1k (Multi-organ Atlas)
- **Unique Feature:** Largest collection, WSI quality, multi-cancer
- **Clinical Annotations:** Variable across cohorts
- **Best For:** Foundation model training, cross-cancer comparisons
- **Includes:**
  - 25 cancer types
  - High-resolution WSI (<1.15 µm/px)
  - Nuclei segmentation
  - Standardized format
- **Limitations:** Gated access, computational requirements, subset is TNBC

---

### 2. Spatial Resolution Comparison

| Dataset | Technology | Spot Size | Resolution | Genes Detected |
|---------|-----------|-----------|------------|----------------|
| GSE210616 | Visium | 55 µm | ~1-10 cells/spot | ~18K |
| Zenodo 14204217 | Custom ST | ~100 µm | ~10-30 cells/spot | Variable |
| Zenodo 14247036 | Visium | 55 µm | ~1-10 cells/spot | ~18K |
| Zenodo 4739739 | Visium | 55 µm | ~1-10 cells/spot | ~18K |
| HEST-1k | Multiple | Variable | Variable | Variable |

---

### 3. Available Data Types

#### Raw Data
| Dataset | Raw Counts | Filtered Counts | Normalized | Deconvolved |
|---------|-----------|----------------|------------|-------------|
| GSE210616 | ❌ | ✅ | ❌ | ❌ |
| Zenodo 14204217 | ✅ | ✅ | ✅ (batch-corrected) | ✅ |
| Zenodo 14247036 | ✅ | ✅ | ✅ (SCTransform) | ✅ |
| Zenodo 4739739 | ✅ | ✅ | ❌ | ❌ |
| HEST-1k | ✅ | ✅ | Variable | Variable |

#### Image Data
| Dataset | H&E | Annotations | IHC | WSI | Image Format |
|---------|-----|-------------|-----|-----|--------------|
| GSE210616 | ✅ | Basic | ❌ | ❌ | JPG, PNG |
| Zenodo 14204217 | ✅ | Detailed (18 types) | ✅ (CD3/CD20) | ✅ | NDPI, JPG, PNG |
| Zenodo 14247036 | ✅ | Moderate | ❌ | ❌ | Standard |
| Zenodo 4739739 | ✅ | Pathologist | ❌ | ❌ | PDF |
| HEST-1k | ✅ | Variable | ❌ | ✅ | High-res |

#### Metadata
| Dataset | Clinical | Survival | Treatment | Demographics | Batch Info |
|---------|----------|----------|-----------|--------------|-----------|
| GSE210616 | ✅ | ❌ | ❌ | ✅ (Race) | ✅ |
| Zenodo 14204217 | ✅ | ✅ | ✅ | ✅ | ✅ |
| Zenodo 14247036 | ✅ | ❌ | ❌ | ✅ | ✅ |
| Zenodo 4739739 | ✅ | ❌ | ❌ | ✅ | ✅ |
| HEST-1k | ✅ | Variable | Variable | Variable | ✅ |

---

### 4. Analysis Resources

#### Available Code
| Dataset | Preprocessing | Analysis | Visualization | Language |
|---------|---------------|----------|---------------|----------|
| GSE210616 | ❌ | ❌ | ❌ | - |
| Zenodo 14204217 | ✅ | ✅ | ✅ | R |
| Zenodo 14247036 | ✅ | ✅ | ✅ | R/Python |
| Zenodo 4739739 | ❌ | ❌ | ❌ | - |
| HEST-1k | ✅ | ✅ | ✅ | Python |

**Code Repositories:**
- **Zenodo 14204217:** [BCTL-Bordet/ST](https://github.com/BCTL-Bordet/ST)
- **Zenodo 14247036:**
  - [cnio-bu/SSc-breast](https://github.com/cnio-bu/SSc-breast)
  - [cnio-bu/ST-preprocess](https://github.com/cnio-bu/ST-preprocess)
  - [cnio-bu/breast-bcspatial](https://github.com/cnio-bu/breast-bcspatial)
- **HEST-1k:** [mahmoodlab/HEST](https://github.com/mahmoodlab/HEST)

---

### 5. Scientific Focus

#### Research Applications Matrix

| Research Question | Best Dataset(s) | Why? |
|-------------------|-----------------|------|
| **Racial disparities in TNBC** | GSE210616 | Only dataset with AA/Caucasian comparison |
| **Tumor microenvironment heterogeneity** | Zenodo 14204217 | Extensive annotations, clustering |
| **Drug response prediction** | Zenodo 14247036 | Beyondcell analysis, >1,200 drugs |
| **Immune landscape** | Wu et al. (4739739) | Detailed immune phenotyping, ecotypes |
| **TLS identification** | Zenodo 14204217 | Dedicated TLS analysis pipeline |
| **Foundation model training** | HEST-1k | Large-scale, diverse, standardized |
| **Clonal heterogeneity** | Zenodo 14247036 | SCEVAN clonal composition |
| **Spatial niches** | Wu et al. (4739739) | Ecotype analysis, niche characterization |
| **Cross-cancer comparison** | HEST-1k | Multi-cancer, standardized |
| **Method development** | Zenodo 14204217 | Multiple analysis pipelines, benchmarking |

---

### 6. Integration Potential

#### Datasets with Matched omics

**Zenodo 14247036:**
- ✅ Matched scRNA-seq reference
- ✅ Same samples profiled
- Best for: Cell type deconvolution validation

**Wu et al. (Zenodo 4739739):**
- ✅ Extensive scRNA-seq atlas (SCP1039)
- ✅ CITE-seq protein data
- ✅ Same patient samples
- Best for: Multi-modal integration

**Zenodo 14204217:**
- ✅ Bulk RNA-seq from same samples
- ✅ External validation cohorts (SCAN-B, METABRIC)
- ✅ IHC validation (CD3/CD20)
- Best for: Cross-platform validation

---

### 7. Quality Metrics

| Dataset | Documentation | Reproducibility | Community Use | Updates |
|---------|--------------|----------------|---------------|---------|
| GSE210616 | Fair | Limited | Low | Static |
| Zenodo 14204217 | Excellent | High | Growing | Active (v3) |
| Zenodo 14247036 | Good | High | New | Recent (v2) |
| Zenodo 4739739 | Good | Moderate | High | Static |
| HEST-1k | Excellent | High | Growing | Active |

**Documentation Quality:**
- **Excellent:** Comprehensive README, detailed metadata, code
- **Good:** README with essential info, some code
- **Fair:** Basic information, minimal documentation

---

## Use Case Recommendations

### For Beginners
**Recommended:** Wu et al. (Zenodo 4739739)
- Smaller size (920 MB)
- Standard Visium format
- Well-documented publication
- Clear pathologist annotations

### For Comprehensive Analysis
**Recommended:** Zenodo 14204217
- Most complete dataset
- Multiple analysis types included
- Extensive clinical data
- Survival outcomes available

### For Drug Discovery
**Recommended:** Zenodo 14247036
- Drug sensitivity predictions
- Functional pathway scores
- TME characterization
- Therapeutic heterogeneity focus

### For ML/AI Development
**Recommended:** HEST-1k
- Large scale
- Standardized format
- High-quality WSI
- Multi-cancer diversity

### For Population Studies
**Recommended:** GSE210616
- Racial diversity
- Uniform cohort (all TNBC)
- Multiple patients

---

## Download Size & Time Estimates

| Dataset | Compressed | Uncompressed | Download (100 Mbps) | Storage Required |
|---------|-----------|--------------|---------------------|------------------|
| GSE210616 | 35 GB | ~50 GB | ~47 min | 60 GB |
| Zenodo 14204217 | 58 GB | ~100 GB | ~78 min | 120 GB |
| Zenodo 14247036 | 7 GB | ~15 GB | ~9 min | 25 GB |
| Zenodo 4739739 | 920 MB | ~2 GB | ~1 min | 5 GB |
| HEST-1k | Variable | Variable | Hours | TBs |

**Recommendations:**
- Start with smaller datasets (Wu et al.) to test workflows
- Use university/institution connections for faster downloads
- Consider downloading to HPC/cluster storage directly
- Plan for 2-3x storage space for compressed + uncompressed + analysis

---

## Citation Impact

| Dataset | Publication | Citations (approx.) | Journal | Impact Factor |
|---------|-------------|---------------------|---------|---------------|
| GSE210616 | PMID: 36283023 | Low | N/A | N/A |
| Zenodo 14204217 | Nature Comm 2024 | New | Nature Commun | 16.6 |
| Zenodo 14247036 | Preprint 2024 | New | N/A | N/A |
| Zenodo 4739739 | Nat Genet 2021 | >500 | Nature Genetics | 30.8 |
| HEST-1k | arXiv 2024 | Growing | N/A | N/A |

*Citation counts as of January 2026

---

## Quick Decision Tree

```
Do you need TNBC-specific data?
├─ Yes
│  ├─ Need >20 samples? → GSE210616 or Zenodo 14204217
│  ├─ Need drug predictions? → Zenodo 14247036
│  ├─ Need extensive annotations? → Zenodo 14204217
│  └─ Want matched scRNA-seq? → Zenodo 4739739 or 14247036
│
├─ No (reference/normal tissue)
│  └─ Human Breast Cell Atlas
│
└─ Multi-cancer/large-scale?
   └─ HEST-1k

Need to start small/test pipeline?
└─ Start with Wu et al. (920 MB)

Building ML models?
└─ HEST-1k

Studying tumor heterogeneity?
└─ Zenodo 14204217

Studying immune landscape?
└─ Wu et al. (with scRNA-seq integration)
```

---

## Technical Requirements

### Computational Resources (Minimum)

| Analysis Type | CPU | RAM | Storage | Software |
|---------------|-----|-----|---------|----------|
| Basic exploration | 4 cores | 16 GB | 100 GB | R/Python |
| Deconvolution | 8 cores | 32 GB | 200 GB | R + packages |
| Full analysis | 16 cores | 64 GB | 500 GB | HPC access |
| ML training (HEST) | 32+ cores | 128+ GB | 2+ TB | GPU required |

### Software Dependencies

**R Packages:**
- Seurat (≥4.0)
- Giotto, SpatialExperiment
- ggplot2, patchwork
- For Belgian dataset: custom STstuff package

**Python Packages:**
- Scanpy, Squidpy
- pandas, numpy
- matplotlib, seaborn
- For HEST: torch, timm

---

This comparison should help you choose the most appropriate dataset for your research needs!
