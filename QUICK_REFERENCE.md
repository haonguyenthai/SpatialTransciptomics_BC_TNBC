# Quick Reference Table

Fast lookup table for TNBC spatial transcriptomics datasets.

## Dataset Summary

| Dataset | Link | Samples | Platform | Size | License |
|---------|------|---------|----------|------|---------|
| **USC TNBC** | [GSE210616](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616) | 22 TNBC patients (43 sections) | Visium | 35 GB | Public |
| **Belgian TNBC** | [Zenodo 14204217](https://zenodo.org/records/14204217) | Multiple TNBC patients | Custom ST | 58 GB | CC BY 4.0 |
| **CNIO Drug** | [Zenodo 14247036](https://zenodo.org/records/14247036) | 4 TNBC (9 total) | Visium | 7 GB | CC BY 4.0 |
| **Wu et al.** | [Zenodo 4739739](https://zenodo.org/records/4739739) | 4 TNBC (6 total) | Visium | 920 MB | CC BY 4.0 |
| **HBCA** | [CellxGene](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637) | 0 (normal tissue) | Multi-modal | Variable | CC BY 4.0 |
| **HEST-1k** | [Hugging Face](https://huggingface.co/datasets/MahmoodLab/hest) | Subset of 1,255 | Multi-platform | >100B | CC BY-NC-SA 4.0 |

## Download Commands

### Quick Downloads (Copy & Paste)

```bash
# Wu et al. (recommended for beginners - smallest dataset)
wget https://zenodo.org/records/4739739/files/filtered_count_matrices.tar.gz
wget https://zenodo.org/records/4739739/files/spatial.tar.gz
wget https://zenodo.org/records/4739739/files/metadata.tar.gz

# CNIO Drug Response Study
wget https://zenodo.org/records/14247036/files/Breast-bcSpatial.zip

# Belgian TNBC Atlas (large - 58 GB)
wget https://zenodo.org/records/14204217/files-archive
```

### Using Helper Script

```bash
# Download our helper script
wget https://raw.githubusercontent.com/[your-repo]/scripts/download_datasets.py

# Install dependencies
pip install requests

# Download Wu et al. dataset
python download_datasets.py --dataset wu --output ./data

# List all available datasets
python download_datasets.py --list
```

## Key Publications

| Dataset | Publication | Journal | Year | PMID/DOI |
|---------|-------------|---------|------|----------|
| USC TNBC | Bassiouni et al. | N/A | 2022 | [36283023](https://pubmed.ncbi.nlm.nih.gov/36283023/) |
| Belgian TNBC | Venet et al. | Nature Commun | 2024 | [10.5281/zenodo.14204217](https://doi.org/10.5281/zenodo.14204217) |
| CNIO Drug | Jiménez-Santos et al. | Preprint | 2024 | [10.5281/zenodo.14247036](https://doi.org/10.5281/zenodo.14247036) |
| Wu et al. | Wu et al. | Nature Genetics | 2021 | [s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1) |
| HBCA | Kumar et al. | Nature | 2024 | In press |
| HEST-1k | Jaume et al. | arXiv | 2024 | [2406.16192](https://arxiv.org/abs/2406.16192) |

## Data Types Available

| Dataset | Raw Counts | H&E Images | Annotations | Clinical | scRNA-seq | Code |
|---------|-----------|-----------|-------------|----------|-----------|------|
| USC TNBC | ❌ | ✅ | Basic | ✅ | ❌ | ❌ |
| Belgian TNBC | ✅ | ✅ | Extensive (18) | ✅ | ❌ | ✅ |
| CNIO Drug | ✅ | ✅ | Moderate | ✅ | ✅ | ✅ |
| Wu et al. | ✅ | ✅ | Pathologist | ✅ | ✅ (SCP1039) | ❌ |
| HBCA | ✅ | ✅ | Variable | ✅ | ✅ | ✅ |
| HEST-1k | ✅ | ✅ (WSI) | Variable | ✅ | Variable | ✅ |

## Key Features Comparison

### Best For...

| Use Case | Recommended Dataset | Reason |
|----------|---------------------|--------|
| Beginners | Wu et al. | Small (920 MB), standard format |
| Racial disparities | USC TNBC | AA vs Caucasian comparison |
| Tumor heterogeneity | Belgian TNBC | 18 annotation categories, extensive clustering |
| Drug response | CNIO Drug | >1,200 drug predictions |
| Immune landscape | Wu et al. | Immune ecotypes, CITE-seq integration |
| TLS studies | Belgian TNBC | Dedicated TLS identification |
| ML/AI development | HEST-1k | Large-scale, standardized |
| Normal reference | HBCA | 714K cells, 58 cell states |

### Sample Sizes

| Dataset | TNBC Patients | Total Sections | Sections per Patient |
|---------|---------------|----------------|---------------------|
| USC TNBC | 22 | 43 | 2 (mostly) |
| Belgian TNBC | Multiple | Multiple | Variable |
| CNIO Drug | 4 | 4 | 1 |
| Wu et al. | 4 | 4 | 1 |

### Platform Details

| Dataset | Technology | Spot Size | Approx Cells/Spot | Genes |
|---------|-----------|-----------|-------------------|-------|
| USC TNBC | Visium | 55 µm | 1-10 | ~18K |
| Belgian TNBC | Custom ST | ~100 µm | 10-30 | Variable |
| CNIO Drug | Visium | 55 µm | 1-10 | ~18K |
| Wu et al. | Visium | 55 µm | 1-10 | ~18K |
| HEST-1k | Multiple | Variable | Variable | Variable |

## Code Repositories

| Dataset | GitHub Repository | Language |
|---------|------------------|----------|
| Belgian TNBC | [BCTL-Bordet/ST](https://github.com/BCTL-Bordet/ST) | R |
| CNIO Drug | [cnio-bu/breast-bcspatial](https://github.com/cnio-bu/breast-bcspatial) | R/Python |
| HBCA | [navinlabcode/HumanBreastCellAtlas](https://github.com/navinlabcode/HumanBreastCellAtlas) | R/Python |
| HEST-1k | [mahmoodlab/HEST](https://github.com/mahmoodlab/HEST) | Python |

## Processing Tools

### Required Software

| Analysis | R Packages | Python Packages |
|----------|-----------|----------------|
| Basic QC | Seurat | scanpy, squidpy |
| Visualization | ggplot2, patchwork | matplotlib, seaborn |
| Deconvolution | Stereoscope, SPOTlight | cell2location, RCTD |
| Spatial stats | spatstat | squidpy |
| Image processing | EBImage | opencv, scikit-image |

### Recommended Workflows

**For Visium Data:**
```r
# R
library(Seurat)
obj <- Load10X_Spatial("path/to/data")
obj <- SCTransform(obj)
obj <- RunPCA(obj)
SpatialFeaturePlot(obj, features = "gene_name")
```

```python
# Python
import scanpy as sc
import squidpy as sq
adata = sc.read_visium("path/to/data")
sc.pp.normalize_total(adata)
sq.pl.spatial_scatter(adata, color="gene_name")
```

## File Formats

### Common File Types

| Extension | Description | Load With |
|-----------|-------------|-----------|
| .h5 | HDF5 matrix | Seurat, scanpy |
| .mtx | Matrix Market | Seurat, scanpy |
| .tsv/.csv | Text matrix | read.table, pandas |
| .rds | R object | readRDS |
| .tar.gz | Compressed archive | tar -xzf |
| .json | Spatial metadata | jsonlite, json |
| .png/.jpg | Images | EBImage, PIL |

## Disk Space Requirements

| Dataset | Download | Extracted | Working | Total Recommended |
|---------|----------|-----------|---------|-------------------|
| Wu et al. | 920 MB | 2 GB | 3 GB | 6 GB |
| CNIO Drug | 7 GB | 15 GB | 10 GB | 35 GB |
| USC TNBC | 35 GB | 50 GB | 30 GB | 120 GB |
| Belgian TNBC | 58 GB | 100 GB | 50 GB | 210 GB |

*Working space includes temporary files and analysis outputs

## Download Time Estimates

| Dataset | 10 Mbps | 100 Mbps | 1 Gbps |
|---------|---------|----------|--------|
| Wu et al. | 12 min | 1 min | <1 min |
| CNIO Drug | 94 min | 9 min | <1 min |
| USC TNBC | 467 min | 47 min | 5 min |
| Belgian TNBC | 773 min | 78 min | 8 min |

## Citation Templates

### BibTeX (Quick Copy)

**Wu et al.:**
```bibtex
@article{wu2021single,
  title={A single-cell and spatially resolved atlas of human breast cancers},
  author={Wu, Sunny Z and others},
  journal={Nature Genetics},
  volume={53},
  pages={1334--1347},
  year={2021}
}
```

**Belgian TNBC:**
```bibtex
@article{venet2024spatial,
  title={Spatial transcriptomics reveals substantial heterogeneity in TNBC},
  author={Venet, David and others},
  journal={Nature Communications},
  year={2024}
}
```

## Contact & Support

| Dataset | Support Contact | GitHub Issues |
|---------|----------------|---------------|
| USC TNBC | [GEO Support](mailto:geo@ncbi.nlm.nih.gov) | N/A |
| Belgian TNBC | [GitHub Issues](https://github.com/BCTL-Bordet/ST/issues) | ✅ |
| CNIO Drug | [GitHub Issues](https://github.com/cnio-bu/breast-bcspatial/issues) | ✅ |
| Wu et al. | [SCP Support](https://singlecell.broadinstitute.org/single_cell) | N/A |
| HEST-1k | [HF Discussions](https://huggingface.co/datasets/MahmoodLab/hest/discussions) | ✅ |

## Access Requirements

| Dataset | Registration | Approval | Access Time |
|---------|--------------|----------|-------------|
| USC TNBC | ❌ | ❌ | Immediate |
| Belgian TNBC | ❌ | ❌ | Immediate |
| CNIO Drug | ❌ | ❌ | Immediate |
| Wu et al. | ❌ | ❌ | Immediate |
| HBCA | ❌ | ❌ | Immediate |
| HEST-1k | ✅ | ✅ | 1-2 days |

## Data Update Status

| Dataset | Version | Last Updated | Active Development |
|---------|---------|--------------|-------------------|
| USC TNBC | v1 | Nov 2022 | ❌ |
| Belgian TNBC | v3 | Sep 2025 | ✅ |
| CNIO Drug | v2 | Nov 2024 | ✅ |
| Wu et al. | v1 | Jun 2021 | ❌ |
| HBCA | Latest | Ongoing | ✅ |
| HEST-1k | Latest | Jan 2026 | ✅ |

---

## Additional Notes

### Data Privacy
- All datasets contain de-identified patient data
- Follow institutional IRB requirements
- Respect data use agreements
- Do not attempt to re-identify patients

### Computational Needs
- Minimum: 16 GB RAM, 4 cores
- Recommended: 32+ GB RAM, 8+ cores
- For HEST-1k: GPU with 16+ GB VRAM

### Community Resources
- **Seurat Forum:** [satijalab.org/seurat](https://satijalab.org/seurat)
- **Scanpy Discourse:** [scanpy.discourse.group](https://scanpy.discourse.group/)
- **10x Support:** [support.10xgenomics.com](https://support.10xgenomics.com/)

---

Last updated: January 2026
