# Axon_Skel_Analyzer

* **Developed by:** Thomas & Héloïse
* **Developed for:** Sandra
* **Team:** Fuchs
* **Date:** January 2025
* **Software:** Fiji

### Images description

3D images taken with x40 objective on a confocal microscope.

3 channels:
  1. *405:* DAPI nuclei
  2. *561:* TUJ neurites
  3. *642:* ORF1p cell bodies

### Macro description

* Segment nuclei using max intensity projection + background subtraction + median filtering + Huang thresholding + holes filling + watershed splitting + size and circularity sorting
* Segment neurites using max intensity projection + background subtraction + median filtering + Huang thresholding + median filtering closing + holes filling + size sorting
* Segment cell bodies using sum slices projection + median filtering + Huang thresholding + opening + dilation + holes filling + size sorting
* Clear cell bodies in neurites binary mask
* Skeletonize neurites, filter out small branches, and analyze skeleton
* Analyze neurites local thickness

### Output

**1 ZIP file:**
* *..._nuclei.zip*: nuclei ROIs

**3 TIF images:**
* *..._filaments.tif*: neurites binary mask 
* *..._filaments_skel.tif*: neurites skeleton
* *..._filaments_locThk.tif*: neurites local thickness

**3 CSV files:**
* *results_global.csv*: nuclei nb, filaments total area, filaments branches nb, filaments branches total length, filaments branches mean diam, filaments junctions nb (one row per image)
* *results_branches_length.csv*: branch length, branch starting point position, branch end point position (one row per filaments branch)
* *results_branches_diam.csv*: branch diam (one row per filaments branch)

### Dependencies

GDSC Fiji plugin (installable via Update Sites)

### Version history

*DAPI_TUJ_ORF1p.ijm*: version 1 released on January 13, 2025.

*DAPI_TAU_TUJ.ijm*: version 2 released on May 22, 2025.

Adapted to another set of 3D images taken with x40 objective on a confocal microscope.

3 channels:
  1. DAPI nuclei
  2. TAU axons
  3. TUJ neurites
     
Key differences from version 1:
 - Cell bodies and axons segmented in TAU channel
 - Local thickness analysis of axons performed on TUJ channel, which labels the full thickness of neuronal processes
