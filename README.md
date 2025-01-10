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

**channel 1:**
* Max intensity projection and give the number of nucleus with simple analyse workflow (background substraction, median, size and circularity filtering)
**channel 3:**
* Sum slices and detect the body cell and remove it (background substraction, median, post processing(Open and Dilate), size filtering)
**channel 2:**
* Max intensity projection, preprocessing, 2 Tubeness Filters (based on DoG) one small to keep small axon and one bigger to keep bigger axon. Concatenate both images, max intensity projection, Huang threshold.
* Run Local thickness and normalize results by dividing pixel value by thickness.
* Run skeletonize and keep only filaments bigger than 1000 pixels

### Output

**1 *zip* file:**
* ROIs corresponding to nuclei

**3 *tif* images:**
* mask_labeled-skeleton: mask of all filaments bigger than 1000 pixels
* mask_LocThk: mask of the Local Thickness in all the detected axons
* mask_nucleus: mask and label of segmented nucleus
* mask_tagged_skeleton: skeleton over 1000 pixels with junctions and end-points
* skeleton_locThk: skeleton with the local thickness values

**3 *csv* files:**
* results_Branches: all informations corresponding o the branches
* results_Thickness: all thickness value normalized and nucleus number but also axons area

### Dependencies

None

### Version history

Version 1 released on January 10, 2025.
