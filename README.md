# Axonal Cluster Workflow  

**Automated ImageJ and Python/Jupyter Notebook pipeline for segmentation and spatial analysis of BASP1 and PLPPR3 protein clusters in hippocampal neuron axons.**  

This workflow, developed for spatial analysis of **PLPPR3**, **BASP1**, **Synaptophysin-1**, and **Drebrin** in hippocampal neurons, was published in *bioRxiv* by Kroon *et al.* ([DOI:10.1101/2024.03.11.584206](https://doi.org/10.1101/2024.03.11.584206)).  
For detailed parameter settings and processing rationale, see the Methods section of Kroon et al. .

---

## Overview  
This repository provides a reproducible workflow for analyzing axonal protein clusters from confocal microscopy images.  
It combines:  

- **ImageJ macro** for automated segmentation of axons and protein clusters  
- **Python/Jupyter Notebook scripts** for spatial measurements and nearest-neighbor analysis  

---

## ImageJ Macro  

The macro processes **raw 2D fluorescence images** or **maximum-intensity projections** and outputs segmented masks and cluster measurements for axons and clusters.  
The macro’s GUI allows you to:  

- Select input images  
- Specify channels to process  
- Set parameters for background correction, smoothing, thresholding, and morphological filtering  

![Macro GUI](https://github.com/ngimber/Axonal_Cluster_Workflow/blob/main/ImageJ_Macro/GUI.png)  

### Processing Steps  

1. **Axon Segmentation**  
   - Normalize local contrast  
   - Gaussian blur smoothing  
   - Binarization (Triangle method)  
   - Morphological closing (1.41 µm)  
   - Ridge detection  
   - Skeletonization  

2. **Cluster Segmentation**  
   - Rolling ball background subtraction   
   - Gaussian blur (140 nm)  
   - Noise removal (Opening)  
   - Thresholding (Intermodes)  
   - Watershed splitting of connected clusters  

3. **Filtering**  
   - Retain only clusters overlapping with segmented axons  

---

## Python/Jupyter Notebook Analysis  

The provided Python scripts and Jupyter Notebooks perform:  

- Nearest neighbor distance calculations  
- Cluster localization analysis (e.g., PLPPR3 within presynaptic sites)  
- Data visualization and summary statistics  

---

## Requirements  

- [ImageJ/Fiji](https://imagej.net/software/fiji/)  
- Python 3.x  
- Jupyter Notebook  
- Python packages:  
  - `pandas`  
  - `numpy`  
  - `scipy`  
  - `matplotlib`  
  - `seaborn`  

---

## Usage  

### 1. Run the ImageJ Macro  
- Open **Fiji/ImageJ**  
- Load the macro from `ImageJ_Macro/`  
- Provide **raw 2D images** or **maximum-intensity projections** as input  (e.g. .tif)
- Adjust parameters in the GUI  
- Run segmentation  

### 2. Run the Python/Jupyter Analysis  
- Install required packages:  
  ```bash
  pip install pandas numpy scipy matplotlib seaborn
  ```  
- Open the Jupyter Notebook in `Python_Scripts/`  
- Load the macro output masks and run analysis cells  

---

## Citation  
If you use this workflow, please cite:  
> Kroon, T. et al. (2024). *Phosphorylation of PLPPR3 membrane proteins as signaling integrator at neuronal synapses*. *bioRxiv*. [https://doi.org/10.1101/2024.03.11.584206](https://doi.org/10.1101/2024.03.11.584206)  
