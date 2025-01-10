# Axon_Skel_Analyzer

* **Developed by:** Thomas & Héloïse
* **Developed for:** Sandra
* **Team:** Fuchs
* **Date:** January 2025
* **Software:** Fiji

### Images description

3D images taken with 63x objective on a confocal microscope.

3 channels:
  1. *405:* DAPI nuclei
  2. *561:* Tug filaments
  3. *642:* ORF1p cell bodies

### Macro description

* Segment nuclei using max intensity projection + background subtraction + median filtering + Huang thresholding + holes filling + watershed splitting + size and circularity sorting
* Segment filaments using max intensity projection + background subtraction + median filtering + Huang thresholding + median filtering closing + holes filling + size sorting
* Segment cell bodies using sum slices projection + median filtering + Huang thresholding + opening + dilation + holes filling + size sorting
* Clear cell bodies in filaments mask
* Skeletonize filaments, filter out small branches, and analyze skeleton
* Analyze filaments local thickness

### Output

**1 ZIP file:**
* *..._nuclei.zip*: nuclei ROIs

**3 TIF images:**
* *..._filaments.tif*: filaments binary mask 
* *..._filaments_skel.tif*: filaments skeleton
* *..._filaments_locThk.tif*: filaments local thickness

**3 CSV files:**
* *results_global.csv*: one row per image: nuclei nb, filaments total area, filaments branches nb, filaments branches total length, filaments branches mean diam, filaments junctions nb
* *results_branches_length.csv*: one row per branch: branch length, branch starting point position, branch end point position
* *results_branches_diam.csv*: one row per branch: branch diam

### Dependencies

None

### Version history

Version 1 released on January 10, 2025.
