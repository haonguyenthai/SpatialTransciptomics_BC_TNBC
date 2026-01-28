# Quickstart Guide

Practical guide to get started with TNBC spatial transcriptomics datasets, including download instructions, loading data, and basic analysis examples.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Downloading Data](#downloading-data)
3. [Loading Data](#loading-data)
4. [Basic Analysis Examples](#basic-analysis-examples)
5. [Visualization](#visualization)
6. [Common Issues](#common-issues)

---

## Prerequisites

### Software Installation

#### R (Recommended for Visium data)

```r
# Install Seurat and dependencies
install.packages("Seurat")
install.packages("hdf5r")  # For reading h5 files
install.packages("ggplot2")
install.packages("patchwork")
install.packages("dplyr")

# Optional: Spatial-specific packages
install.packages("BiocManager")
BiocManager::install("SpatialExperiment")
BiocManager::install("ggspavis")

# For deconvolution
devtools::install_github("almaan/stereoscope")
```

#### Python (Alternative)

```bash
pip install scanpy squidpy anndata
pip install matplotlib seaborn
pip install pandas numpy scipy
pip install h5py

# For spatial analysis
pip install spatial-omics
```

---

## Downloading Data

### Option 1: Direct Download (Zenodo datasets)

```bash
# Create data directory
mkdir -p ~/tnbc_spatial_data
cd ~/tnbc_spatial_data

# Wu et al. dataset (smallest, good for testing)
wget https://zenodo.org/records/4739739/files/filtered_count_matrices.tar.gz
wget https://zenodo.org/records/4739739/files/spatial.tar.gz
wget https://zenodo.org/records/4739739/files/metadata.tar.gz

tar -xzf filtered_count_matrices.tar.gz
tar -xzf spatial.tar.gz
tar -xzf metadata.tar.gz

# Or use zenodo_get for easier batch download
pip install zenodo_get
zenodo_get 10.5281/zenodo.4739739
```

### Option 2: GEO Download (GSE210616)

```bash
# Using wget
wget -r -np -nd -A "*.gz,*.csv,*.h5" \
  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210616/
```

```r
# Or using GEOquery in R
library(GEOquery)
gse <- getGEO("GSE210616")
```

### Option 3: Programmatic Download

```python
# Python script for automated download
import requests
import os

def download_file(url, output_path):
    """Download file with progress tracking"""
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    
    with open(output_path, 'wb') as f:
        downloaded = 0
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                print(f"\rDownloaded {downloaded/1e6:.1f}/{total_size/1e6:.1f} MB", end='')
    print("\nDownload complete!")

# Example: Download Wu et al. data
base_url = "https://zenodo.org/records/4739739/files/"
files = [
    "filtered_count_matrices.tar.gz",
    "spatial.tar.gz",
    "metadata.tar.gz"
]

for file in files:
    download_file(base_url + file, file)
```

---

## Loading Data

### Loading Visium Data (10x format)

#### R (Seurat)

```r
library(Seurat)
library(ggplot2)

# Load a single sample
sample_path <- "path/to/sample_directory"
spatial_obj <- Load10X_Spatial(
    data.dir = sample_path,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1"
)

# Check data
print(spatial_obj)
dim(spatial_obj)  # Genes x Spots

# View spatial plot
SpatialPlot(spatial_obj, features = "nCount_Spatial")
```

#### Python (Scanpy/Squidpy)

```python
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

# Load Visium data
adata = sc.read_visium("path/to/sample_directory")

# Check data
print(adata)
print(adata.shape)  # Spots x Genes

# Basic quality control
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Spatial plot
sq.pl.spatial_scatter(adata, color="total_counts")
plt.show()
```

### Loading Multiple Samples

```r
# R: Load and merge multiple samples
library(Seurat)

sample_dirs <- c(
    "CID44971_section1",
    "CID44971_section2",
    "CID4465_section1"
)

# Load all samples
spatial_list <- lapply(sample_dirs, function(dir) {
    Load10X_Spatial(
        data.dir = dir,
        filename = "filtered_feature_bc_matrix.h5",
        slice = basename(dir)
    )
})

# Merge samples
spatial_merged <- merge(
    spatial_list[[1]], 
    spatial_list[2:length(spatial_list)],
    add.cell.ids = sample_dirs
)
```

### Loading Custom Formats (Belgian Dataset)

```r
# Load Belgian TNBC dataset
library(data.table)

# Load count matrix
counts <- fread("rawCountsMatrices/sample_001.tsv.gz")
counts_matrix <- as.matrix(counts[, -1])
rownames(counts_matrix) <- counts$V1

# Load spot coordinates
spots <- readRDS("byArray/sample_001/allSpots.RDS")

# Load annotations
annotations <- readRDS("byArray/sample_001/annotBySpot.RDS")

# Create Seurat object
spatial_obj <- CreateSeuratObject(
    counts = counts_matrix,
    assay = "Spatial"
)

# Add spatial coordinates
spatial_coords <- data.frame(
    x = spots$pixel_x,
    y = spots$pixel_y
)
rownames(spatial_coords) <- colnames(counts_matrix)

# Add coordinates to object
spatial_obj@images$slice1 <- new(
    Class = "VisiumV1",
    coordinates = spatial_coords
)
```

---

## Basic Analysis Examples

### Example 1: Quality Control

```r
library(Seurat)

# Calculate QC metrics
spatial_obj[["percent.mt"]] <- PercentageFeatureSet(
    spatial_obj, 
    pattern = "^MT-"
)

# Visualize QC metrics
VlnPlot(
    spatial_obj, 
    features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
    ncol = 3, 
    pt.size = 0
)

# Spatial QC plots
SpatialFeaturePlot(
    spatial_obj,
    features = c("nCount_Spatial", "nFeature_Spatial")
)

# Filter spots
spatial_obj <- subset(
    spatial_obj,
    subset = nFeature_Spatial > 200 & 
             nFeature_Spatial < 7000 & 
             percent.mt < 20
)
```

### Example 2: Normalization and Clustering

```r
# SCTransform normalization (recommended for spatial data)
spatial_obj <- SCTransform(spatial_obj, assay = "Spatial")

# PCA and clustering
spatial_obj <- RunPCA(spatial_obj, assay = "SCT")
spatial_obj <- FindNeighbors(spatial_obj, reduction = "pca", dims = 1:30)
spatial_obj <- FindClusters(spatial_obj, resolution = 0.5)
spatial_obj <- RunUMAP(spatial_obj, reduction = "pca", dims = 1:30)

# Visualize clusters
DimPlot(spatial_obj, reduction = "umap", label = TRUE)
SpatialDimPlot(spatial_obj, label = TRUE, label.size = 3)
```

### Example 3: Marker Gene Visualization

```r
# Define TNBC-relevant markers
markers <- list(
    Tumor = c("KRT8", "KRT18", "EPCAM"),
    Immune = c("PTPRC", "CD3E", "CD68"),
    Stroma = c("COL1A1", "DCN", "ACTA2"),
    Proliferation = c("MKI67", "TOP2A")
)

# Spatial visualization
SpatialFeaturePlot(
    spatial_obj,
    features = c("KRT8", "PTPRC", "COL1A1", "MKI67"),
    alpha = c(0.1, 1)
)

# Violin plots by cluster
VlnPlot(
    spatial_obj,
    features = c("KRT8", "CD3E", "COL1A1"),
    pt.size = 0
)
```

### Example 4: Spot Deconvolution

```r
library(SPOTlight)  # Or stereoscope

# Assuming you have scRNA-seq reference
# Load scRNA reference
sc_ref <- readRDS("path/to/scrnaseq_reference.rds")

# Run deconvolution
deconv_results <- spotlight_deconvolution(
    se_sc = sc_ref,
    se_st = spatial_obj,
    groups = "celltype",
    mgs = marker_genes,
    hvg = 3000
)

# Add cell type proportions to spatial object
cell_types <- deconv_results$mat

# Visualize cell type proportions
for (ct in colnames(cell_types)) {
    spatial_obj[[ct]] <- cell_types[, ct]
}

SpatialFeaturePlot(
    spatial_obj,
    features = colnames(cell_types)[1:4]
)
```

### Example 5: Spatial Statistics

```python
import scanpy as sc
import squidpy as sq

# Load data
adata = sc.read_visium("path/to/sample")

# Preprocess
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

# Compute spatial neighbors
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

# Calculate spatial autocorrelation (Moran's I)
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=["KRT8", "PTPRC", "COL1A1"]
)

# Neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key="clusters")
sq.pl.nhood_enrichment(adata, cluster_key="clusters")

# Co-occurrence score
sq.gr.co_occurrence(adata, cluster_key="clusters")
sq.pl.co_occurrence(
    adata,
    cluster_key="clusters",
    clusters=["Tumor", "Immune", "Stroma"]
)
```

---

## Visualization

### Basic Spatial Plots

```r
# Gene expression overlay
SpatialFeaturePlot(
    spatial_obj,
    features = "KRT8",
    pt.size.factor = 1.6,
    alpha = c(0.1, 1)
) + 
    scale_fill_viridis_c() +
    theme(legend.position = "right")

# Multiple genes side-by-side
SpatialFeaturePlot(
    spatial_obj,
    features = c("KRT8", "PTPRC", "COL1A1", "MKI67"),
    ncol = 2,
    alpha = c(0.1, 1)
)

# Cluster overlay
SpatialDimPlot(
    spatial_obj,
    cells.highlight = CellsByIdentities(
        spatial_obj, 
        idents = c("0", "1", "2")
    ),
    facet.highlight = TRUE
)
```

### Advanced Visualizations

```python
import squidpy as sq
import matplotlib.pyplot as plt

# Create multi-panel spatial plot
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

# Gene expression
sq.pl.spatial_scatter(
    adata,
    color="KRT8",
    ax=axes[0, 0],
    title="Tumor (KRT8)"
)

sq.pl.spatial_scatter(
    adata,
    color="PTPRC",
    ax=axes[0, 1],
    title="Immune (PTPRC)"
)

# Clusters
sq.pl.spatial_scatter(
    adata,
    color="clusters",
    ax=axes[1, 0],
    title="Spatial Clusters"
)

# Neighborhood graph
sq.pl.spatial_scatter(
    adata,
    connectivity_key="spatial_connectivities",
    ax=axes[1, 1],
    title="Spatial Network"
)

plt.tight_layout()
plt.savefig("spatial_overview.pdf", dpi=300)
```

### Comparing Multiple Samples

```r
# Plot same gene across multiple samples
SpatialFeaturePlot(
    spatial_merged,
    features = "KRT8",
    images = c("slice1", "slice2", "slice3"),
    ncol = 3
)

# Compare cluster distributions
table(spatial_merged$orig.ident, spatial_merged$seurat_clusters)

# Visualization
library(ggplot2)
cluster_freq <- as.data.frame(
    table(spatial_merged$orig.ident, spatial_merged$seurat_clusters)
)
colnames(cluster_freq) <- c("Sample", "Cluster", "Freq")

ggplot(cluster_freq, aes(x = Sample, y = Freq, fill = Cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_minimal() +
    labs(y = "Proportion", title = "Cluster Distribution by Sample")
```

---

## Common Issues

### Issue 1: Memory Problems

```r
# Problem: Large spatial objects crash R
# Solution: Process samples individually or use disk-backed storage

# Save intermediate results
saveRDS(spatial_obj, "spatial_processed.rds")

# Or use on-disk storage
options(Seurat.object.assay.version = "v5")  # Seurat v5 on-disk
```

### Issue 2: Image Alignment

```r
# Problem: H&E image not aligned with spots
# Solution: Adjust scale factors

# Check current scale factors
spatial_obj@images$slice1@scale.factors

# Manually adjust if needed
spatial_obj@images$slice1@scale.factors$lowres <- 0.05
```

### Issue 3: Missing Metadata

```r
# Problem: No clinical metadata in downloaded files
# Solution: Load from supplementary files

# Load clinical data
clinical <- read.csv("clinical_data.csv")
rownames(clinical) <- clinical$sample_id

# Add to Seurat object
spatial_obj <- AddMetaData(
    spatial_obj,
    metadata = clinical[spatial_obj$orig.ident, ]
)
```

### Issue 4: Different Coordinate Systems

```python
# Problem: Coordinates don't match between data sources
# Solution: Transform coordinates

import numpy as np

# Flip y-axis
adata.obsm["spatial"][:, 1] = -adata.obsm["spatial"][:, 1]

# Scale coordinates
adata.obsm["spatial"] = adata.obsm["spatial"] * scale_factor

# Rotate (if needed)
angle = np.pi / 4  # 45 degrees
rotation_matrix = np.array([
    [np.cos(angle), -np.sin(angle)],
    [np.sin(angle), np.cos(angle)]
])
adata.obsm["spatial"] = adata.obsm["spatial"] @ rotation_matrix.T
```

---

## Example Workflow: Complete Analysis

```r
library(Seurat)
library(ggplot2)
library(dplyr)

# 1. Load data
spatial_obj <- Load10X_Spatial("path/to/sample")

# 2. QC
spatial_obj[["percent.mt"]] <- PercentageFeatureSet(spatial_obj, pattern = "^MT-")
spatial_obj <- subset(spatial_obj, 
                      subset = nFeature_Spatial > 200 & 
                              nFeature_Spatial < 7000 & 
                              percent.mt < 20)

# 3. Normalization
spatial_obj <- SCTransform(spatial_obj)

# 4. Dimensionality reduction
spatial_obj <- RunPCA(spatial_obj)
spatial_obj <- RunUMAP(spatial_obj, dims = 1:30)

# 5. Clustering
spatial_obj <- FindNeighbors(spatial_obj, dims = 1:30)
spatial_obj <- FindClusters(spatial_obj, resolution = 0.5)

# 6. Find markers
markers <- FindAllMarkers(spatial_obj, only.pos = TRUE, min.pct = 0.25)
top_markers <- markers %>% 
    group_by(cluster) %>% 
    top_n(n = 10, wt = avg_log2FC)

# 7. Visualize
SpatialDimPlot(spatial_obj, label = TRUE)
SpatialFeaturePlot(spatial_obj, 
                  features = c("KRT8", "PTPRC", "COL1A1", "MKI67"))

# 8. Save results
saveRDS(spatial_obj, "spatial_analyzed.rds")
write.csv(markers, "cluster_markers.csv")

# 9. Export for further analysis
counts <- GetAssayData(spatial_obj, slot = "counts")
metadata <- spatial_obj@meta.data
coords <- spatial_obj@images$slice1@coordinates

writeMM(counts, "counts_matrix.mtx")
write.csv(metadata, "metadata.csv")
write.csv(coords, "spatial_coords.csv")
```

---

## Next Steps

After completing this quickstart:

1. **Explore Advanced Analysis**
   - Cell-cell interaction analysis
   - Spatially variable genes
   - Trajectory inference
   - Integration with other modalities

2. **Read Publications**
   - Study the methods used in each dataset
   - Understand biological interpretations
   - Learn about TNBC heterogeneity

3. **Join Community**
   - Seurat forums
   - Scanpy discourse
   - Twitter/X: #spatialomics #spatialgenomics

4. **Contribute**
   - Share your analyses
   - Report issues
   - Suggest improvements

---

## Resources

### Documentation
- [Seurat Spatial Vignettes](https://satijalab.org/seurat/articles/spatial_vignette.html)
- [Squidpy Tutorials](https://squidpy.readthedocs.io/)
- [10x Genomics Support](https://support.10xgenomics.com/spatial-gene-expression)

### Papers
- Wu et al., Nature Genetics 2021
- Venet et al., Nature Communications 2024
- Various methods papers for tools used

### Code Repositories
- See DATASET_COMPARISON.md for specific repos

Good luck with your TNBC spatial transcriptomics research!
