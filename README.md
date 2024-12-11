# Axon_Skel_Analyzer

* **Developed by:** Thomas & Héloïse
* **Developed for:** Sandra
* **Team:** Prochiantz
* **Date:** December 2024
* **Software:** Fiji

### Images description

**channel 1:** Nucleus (DAPI)
**channel 2:** axons 
**channel 3:** cell body 

### Macros description

**Axon_Skel_Analyzer.ijm:**

**channel 1:**
* Max intensity projection and give the number of nucleus with simple analyse workflow (background substraction, median, size and circularity filtering)
**channel 3:**
* Sum slices and detect the body cell and remove it (background substraction, median, post processing(Open and Dilate), size filtering)
**channel 2:**
* Max intensity projection, preprocessing, 2 Tubeness Filters (based on DoG) one small to keep small axon and one bigger to keep bigger axon. Concatenate both images, max intensity projection, Huang threshold.
* Run Local thickness and normalize results by dividing pixel value by thickness.
* Run skeletonize and keep only filaments bigger than 1000 pixels

### Output 

**5 images:**
* mask_labeled-skeleton: mask of all filaments bigger than 1000 pixels
* mask_LocThk: mask of the Local Thickness in all the detected axons
* mask_nucleus: mask and label of segmented nucleus
* mask_tagged_skeleton: skeleton over 1000 pixels with junctions and end-points
* skeleton_locThk: skeleton with the local thickness values

**2 csv files:**
* results_Branches: all informations corresponding o the branches
* results_Thickness: all thickness value normalized and nucleus number but also axons area

**1 zip file:**
* Roi corresponding to nucleus

### Dependencies

**Axon_Skel_Analyzer.ijm:** Fiji: up to date

### Version history

Version 1.0 released on December 11, 2024.
