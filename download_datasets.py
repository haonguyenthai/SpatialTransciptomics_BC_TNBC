#!/usr/bin/env python3
"""
TNBC Spatial Transcriptomics Dataset Downloader

This script helps download TNBC spatial transcriptomics datasets
from various sources (Zenodo, GEO, etc.)

Usage:
    python download_datasets.py --dataset wu --output ./data
    python download_datasets.py --dataset all --output ./data
"""

import argparse
import os
import sys
import requests
from pathlib import Path
import subprocess
from typing import List, Dict

# Dataset configurations
DATASETS = {
    "wu": {
        "name": "Wu et al. (Zenodo 4739739)",
        "size_gb": 0.9,
        "tnbc_samples": "4/6 samples",
        "files": {
            "filtered_matrices": "https://zenodo.org/records/4739739/files/filtered_count_matrices.tar.gz",
            "spatial": "https://zenodo.org/records/4739739/files/spatial.tar.gz",
            "metadata": "https://zenodo.org/records/4739739/files/metadata.tar.gz",
            "images": "https://zenodo.org/records/4739739/files/images.pdf"
        }
    },
    "cnio": {
        "name": "CNIO Drug Response (Zenodo 14247036)",
        "size_gb": 7.0,
        "tnbc_samples": "4/9 samples",
        "files": {
            "all_data": "https://zenodo.org/records/14247036/files/Breast-bcSpatial.zip"
        }
    },
    "belgian": {
        "name": "Belgian TNBC Atlas (Zenodo 14204217)",
        "size_gb": 58.0,
        "tnbc_samples": "Multiple TNBC",
        "note": "Large dataset - ensure sufficient storage",
        "files": {
            "clinical": "https://zenodo.org/records/14204217/files/Clinical.tar",
            "clustering": "https://zenodo.org/records/14204217/files/clustering.tar",
            "images": "https://zenodo.org/records/14204217/files/Images.tar",
            "robjects": "https://zenodo.org/records/14204217/files/Robjects.tar",
            "raw_counts": "https://zenodo.org/records/14204217/files/rawCountsMatrices.tar",
            "classification": "https://zenodo.org/records/14204217/files/classification.tar",
            "deconvolution": "https://zenodo.org/records/14204217/files/deconvolution.tar"
        }
    },
    "her2": {
        "name": "Andersson et al. HER2+ (Zenodo 3957257)",
        "size_gb": 0.63,
        "tnbc_samples": "0 (HER2+ reference)",
        "note": "HER2+ dataset - useful for comparative analysis",
        "files": {
            "count_matrices": "https://zenodo.org/records/3957257/files/count-matrices.zip",
            "images": "https://zenodo.org/records/3957257/files/images.zip",
            "metadata": "https://zenodo.org/records/3957257/files/meta.zip",
            "spot_selections": "https://zenodo.org/records/3957257/files/spot-selections.zip"
        }
    },
    "gse210616": {
        "name": "USC TNBC Cohort (GEO GSE210616)",
        "size_gb": 35.0,
        "tnbc_samples": "22 patients",
        "note": "Requires GEO FTP access",
        "method": "geo"
    }
}


def download_file(url: str, output_path: Path, chunk_size: int = 8192) -> bool:
    """
    Download a file with progress tracking.
    
    Args:
        url: URL to download from
        output_path: Path to save file
        chunk_size: Size of chunks to download at a time
        
    Returns:
        True if successful, False otherwise
    """
    try:
        print(f"\nDownloading: {output_path.name}")
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\rProgress: {percent:.1f}% ({downloaded/1e6:.1f}/{total_size/1e6:.1f} MB)", end='')
        
        print("\n✓ Download complete!")
        return True
        
    except requests.exceptions.RequestException as e:
        print(f"\n✗ Download failed: {e}")
        return False
    except KeyboardInterrupt:
        print("\n✗ Download interrupted by user")
        if output_path.exists():
            output_path.unlink()
        return False


def download_geo_dataset(geo_id: str, output_dir: Path) -> bool:
    """
    Download dataset from GEO using wget.
    
    Args:
        geo_id: GEO accession ID (e.g., GSE210616)
        output_dir: Directory to save files
        
    Returns:
        True if successful, False otherwise
    """
    series_prefix = geo_id[:7]  # e.g., GSE2106 from GSE210616
    ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{series_prefix}nnn/{geo_id}/"
    
    print(f"\nDownloading GEO dataset {geo_id}...")
    print(f"This may take a while...")
    
    cmd = [
        "wget",
        "-r",           # recursive
        "-np",          # no parent
        "-nd",          # no directories
        "-P", str(output_dir),
        ftp_url
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print("✓ GEO download complete!")
        return True
    except subprocess.CalledProcessError:
        print("✗ GEO download failed. Make sure wget is installed.")
        return False
    except FileNotFoundError:
        print("✗ wget not found. Please install wget or download manually from GEO.")
        return False


def extract_archive(archive_path: Path) -> bool:
    """
    Extract tar.gz or zip archive.
    
    Args:
        archive_path: Path to archive file
        
    Returns:
        True if successful, False otherwise
    """
    import tarfile
    import zipfile
    
    print(f"\nExtracting: {archive_path.name}")
    
    try:
        if archive_path.suffix == '.zip':
            with zipfile.ZipFile(archive_path, 'r') as zip_ref:
                zip_ref.extractall(archive_path.parent)
        elif archive_path.suffix == '.gz' or archive_path.suffix == '.tar':
            with tarfile.open(archive_path, 'r:*') as tar_ref:
                tar_ref.extractall(archive_path.parent)
        
        print("✓ Extraction complete!")
        return True
        
    except Exception as e:
        print(f"✗ Extraction failed: {e}")
        return False


def download_dataset(dataset_id: str, output_dir: Path, extract: bool = True) -> bool:
    """
    Download a specific dataset.
    
    Args:
        dataset_id: Dataset identifier (wu, cnio, belgian, gse210616)
        output_dir: Directory to save files
        extract: Whether to extract archives after download
        
    Returns:
        True if successful, False otherwise
    """
    if dataset_id not in DATASETS:
        print(f"✗ Unknown dataset: {dataset_id}")
        print(f"Available datasets: {', '.join(DATASETS.keys())}")
        return False
    
    config = DATASETS[dataset_id]
    dataset_dir = output_dir / dataset_id
    dataset_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Dataset: {config['name']}")
    print(f"Size: {config['size_gb']} GB")
    print(f"TNBC samples: {config['tnbc_samples']}")
    if 'note' in config:
        print(f"Note: {config['note']}")
    print(f"{'='*60}")
    
    # Check if GEO dataset
    if config.get('method') == 'geo':
        return download_geo_dataset(dataset_id.upper(), dataset_dir)
    
    # Download files
    success = True
    downloaded_files = []
    
    for file_key, url in config.get('files', {}).items():
        filename = url.split('/')[-1]
        output_path = dataset_dir / filename
        
        if output_path.exists():
            print(f"\n✓ File already exists: {filename}")
            downloaded_files.append(output_path)
            continue
        
        if download_file(url, output_path):
            downloaded_files.append(output_path)
        else:
            success = False
    
    # Extract archives
    if extract and success:
        for file_path in downloaded_files:
            if file_path.suffix in ['.gz', '.tar', '.zip']:
                extract_archive(file_path)
    
    return success


def list_datasets():
    """Print information about available datasets."""
    print("\nAvailable TNBC Spatial Transcriptomics Datasets:\n")
    print(f"{'ID':<12} {'Name':<40} {'Size (GB)':<10} {'TNBC Samples':<15}")
    print("=" * 85)
    
    for dataset_id, config in DATASETS.items():
        print(f"{dataset_id:<12} {config['name']:<40} {config['size_gb']:<10} {config['tnbc_samples']:<15}")
    
    print("\nUsage examples:")
    print("  python download_datasets.py --dataset wu --output ./data")
    print("  python download_datasets.py --dataset belgian --output ./data --no-extract")
    print("  python download_datasets.py --list")


def check_disk_space(required_gb: float) -> bool:
    """
    Check if sufficient disk space is available.
    
    Args:
        required_gb: Required space in GB
        
    Returns:
        True if sufficient space, False otherwise
    """
    import shutil
    
    total, used, free = shutil.disk_usage("/")
    free_gb = free / (2**30)
    
    if free_gb < required_gb * 1.5:  # 1.5x for extraction
        print(f"⚠ Warning: Only {free_gb:.1f} GB free. Need ~{required_gb*1.5:.1f} GB")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Download TNBC spatial transcriptomics datasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List available datasets
  python download_datasets.py --list
  
  # Download Wu et al. dataset (smallest)
  python download_datasets.py --dataset wu --output ./data
  
  # Download Belgian dataset without extracting
  python download_datasets.py --dataset belgian --output ./data --no-extract
  
  # Download multiple datasets
  python download_datasets.py --dataset wu cnio --output ./data
        """
    )
    
    parser.add_argument(
        '--dataset',
        nargs='+',
        choices=list(DATASETS.keys()) + ['all'],
        help='Dataset(s) to download'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('./tnbc_data'),
        help='Output directory (default: ./tnbc_data)'
    )
    parser.add_argument(
        '--no-extract',
        action='store_true',
        help='Do not extract archives after download'
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help='List available datasets and exit'
    )
    
    args = parser.parse_args()
    
    if args.list:
        list_datasets()
        return 0
    
    if not args.dataset:
        print("Error: Please specify --dataset or use --list")
        parser.print_help()
        return 1
    
    # Expand 'all' to all datasets
    datasets_to_download = args.dataset
    if 'all' in datasets_to_download:
        datasets_to_download = list(DATASETS.keys())
    
    # Calculate total size
    total_size = sum(DATASETS[d]['size_gb'] for d in datasets_to_download)
    
    print(f"\nTotal download size: ~{total_size:.1f} GB")
    
    # Check disk space
    if not check_disk_space(total_size):
        response = input("Continue anyway? (y/N): ")
        if response.lower() != 'y':
            return 1
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Download datasets
    success_count = 0
    for dataset in datasets_to_download:
        if download_dataset(dataset, args.output, extract=not args.no_extract):
            success_count += 1
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Download Summary:")
    print(f"  Successful: {success_count}/{len(datasets_to_download)}")
    print(f"  Output directory: {args.output.absolute()}")
    print(f"{'='*60}\n")
    
    return 0 if success_count == len(datasets_to_download) else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\nDownload interrupted by user.")
        sys.exit(1)
