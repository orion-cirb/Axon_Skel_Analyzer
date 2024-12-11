/////////////////////////////////////////////////////////////////
//      Authors Thomas Caille & Héloïse Monnet @ ORION-CIRB    //
//       https://github.com/orion-cirb/MorphOocyte_Nuclei      //
/////////////////////////////////////////////////////////////////


// Ask for the images directory
inputDir = getDirectory("Please select a directory containing images to analyze");

// Create results directory
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
//resultDir = inputDir + "Results_" + year + "-" + (month+1) + "-" + dayOfMonth + "_" + hour + "-" + minute + "-" + second + File.separator();
resultDir = inputDir + "Results" + File.separator();
if (!File.isDirectory(resultDir)) {
	File.makeDirectory(resultDir);
}


// Get all files in the input directory
inputFiles = getFileList(inputDir);




// Create a file named "results.csv" and write headers in it


// Loop through all files with .TIF extension
for (i = 5; i < inputFiles.length; i++) {
    if (endsWith(inputFiles[i], ".nd")) {
    	// get the file name
    	nameNoExt = substring(inputFiles[i],0,lengthOf(inputFiles[i])-3);
    	// create the final result folder for each image
		nameResultDir = inputDir + "Results" + File.separator()+ nameNoExt +File.separator ;
		if (!File.isDirectory(nameResultDir)) {
			File.makeDirectory(nameResultDir);
		}
    	fileResults = File.open(nameResultDir+ "results_Thickness.csv");
		print(fileResults, "Nucleus_Nb,Thickness,Count,Area_axon\n");
		// open the first channel and count the number of nucleus
		run("Bio-Formats Importer", "open=["+inputDir + inputFiles[i]+"] autoscale color_mode=Default specify_range split_channels view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
    	rename("nucleus_number");
    	slices = nSlices;
    	run("Z Project...", "projection=[Max Intensity]");
    	getStatistics(area,mean,min,max,std,histogram);
    	// substract the background
		run("Subtract...", "value="+mean+(std/2));
		run("Median...", "radius=10");
		setAutoThreshold("Huang dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Fill Holes");
		run("Watershed");
		// count only particles larger than 2000 pixels and 0.70 circularity
    	run("Analyze Particles...", "size=2000-Infinity circularity=0.70-1.00 pixel clear add");
    	// save the image and the associated ROIs
    	saveAs("Tiff", nameResultDir + "mask_nucleus");
    	roiManager("Save", nameResultDir + "Roi_nucleus.zip");
		// print into the "results_thickness" file the number of nucleus found (1st line)
    	print(fileResults  , nResults+",,");
		
    	// open the channel 2 and do a max intensity projection
    	run("Bio-Formats Importer", "open=["+inputDir + inputFiles[i]+"] autoscale color_mode=Default specify_range split_channels view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");
    	run("Z Project...", "projection=[Max Intensity]");
    	rename("Zsum_ves");
    	// open the channel 3 and sum the slices, the channel 3 will be use to delete the nucleus from the image 
    	run("Bio-Formats Importer", "open=["+inputDir + inputFiles[i]+"] autoscale color_mode=Default specify_range split_channels view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");
    	rename("rawImage");

    	run("Z Project...", "projection=[Sum Slices]");
    	run("Set Measurements...", "mean standard redirect=None decimal=0");
		List.setMeasurements;
		
		// get the 3rd quartile bgnoise of the 3 channel
		mean = List.getValue("Mean");
		stdDev = List.getValue("StdDev");
		bgNoise = (mean + stdDev) / nSlices;
		List.clear();
		// Remove background noise
		run("Subtract...", "value=" + bgNoise +" stack");
		// apply a big median filter to smooth the image without smoothing edges
		run("Median...", "radius=15");
		// otsu threshold to the channel 3
		setAutoThreshold("Otsu dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		// some post precessing opérations allowing a better segmentation of the nucleus
		run("Options...", "iterations=15 count=1 black do=Open");
		run("Options...", "iterations=5 count=1 black do=Dilate");
		run("Fill Holes");
		// get all the cell body, combine them and add them to selection
		run("Analyze Particles...", "size=18000-Infinity pixel clear add");
		roiManager("Combine");

		run("Add Selection...");		
		
		// go back to the channel 2 and delete the previou cell body selection 
		selectImage("Zsum_ves");
		run("Restore Selection");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
		run("Select None");
		// do some preprocessing step to have a better segmentation 
		selectImage("Zsum_ves");
		run("Subtract Background...", "rolling=50");
		run("Median...", "radius=2");
		// use a small tubeness filter to keep only small axon
		run("Tubeness", "sigma=2");
		rename("tub_small");
		selectImage("Zsum_ves");
		// use a bigger one to keep big axon
		run("Tubeness", "sigma=10");
		rename("tub_big");
		// concatenate both image 
		run("Concatenate...", "keep image1=[tub_small] image2=[tub_big] image3=[-- None --]");
		run("Z Project...", "projection=[Max Intensity]");
		run("Set Measurements...", "area mean standard limit redirect=None decimal=0");
		
  		//threshold the new image 
		setAutoThreshold("Huang dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		// some post-precessing allowing a better segmentation of the axon
		run("Options...", "iterations=4 count=5 black do=Open");
		run("Options...", "iterations=5 count=2 black do=Close");
		rename("bin_mask");
		// measure and register the axon area in the "results_Thickness" file (2nd line)
		List.setMeasurements;
		area_bin = List.getValue("Area");
		print(fileResults  , ", , ,"+area_bin +"\n");
		
		close("\\Others");
		
		// run the local thickness analyse 
		run("Local Thickness (complete process)", "threshold=128");
   		// get some stats 
		getStatistics(area, mean, min, max, std, histogram);
		nBins=max;
		newMax=max-1;
		getHistogram(values, counts, nBins);
		vabs=newArray(newMax);
		norm=newArray(newMax);
		// loops threw the histogram and print results in a .csv file
		for (z=1 ; z<(newMax) ; z++){
			vabs[z] = Math.floor(values[z]);
		// normalize the value by the thickness
			norm[z]= counts[z] / vabs[z];
			
			print(fileResults  ," ,"+ vabs[z]+","+norm[z]+"\n");
		}
		File.close(fileResults);
    
		// make a skeleton of the binary vessel image
		selectImage("bin_mask");
		run("Skeletonize");
		run("Analyze Skeleton (2D/3D)", "prune=none show display");
		
		getStatistics(area, mean, min, max, std, histogram);

		nBins=max;
		getHistogram(values, counts, nBins);
		
		//only take filaments larger than 1000 to eliminate errors
		for (j=1 ; j<nBins ; j++){
			if (counts[j] > 1000) {
				selectImage("bin_mask-labeled-skeletons");
				run("Duplicate...", " ");
				run("Manual Threshold...", "min="+values[j]+" max="+values[j]);
				setOption("BlackBackground", true);
				run("Convert to Mask");
			}
		}
		// check the number of filaments larger than 1000 pixels and concatenate them
		title=Image.title;
		titleNumber=substring(title, 27);
		titleNumber=parseFloat(titleNumber);
		concatCommand = "" ;
		defaultConcatCommand = "keep image1=bin_mask-labeled-skeletons-1";
		if (titleNumber > 1) {
			for (t = 2; t <= titleNumber; t++) {	
				concatCommand = concatCommand + " image" + t + "=bin_mask-labeled-skeletons-"+t;	
			}
			endConcatCommand =defaultConcatCommand +concatCommand + " image" +(titleNumber + 1)+ "=[-- None --]";
			run("Concatenate...", endConcatCommand);
			run("Z Project...", "projection=[Max Intensity]");
		} 
		
		// select filaments and create  a selection on the local thickness image
		run("Create Selection");
		selectImage("bin_mask_LocThk");
		run("Duplicate...", " ");
		//save the local thickness image
		saveAs("Tiff", nameResultDir + "mask_locThk");
		run("Restore Selection");
		setBackgroundColor(0, 0, 0);
		run("Clear Outside");
		
		setMinAndMax(0, 255);
		run("8-bit");
		run("Select None");
		// save the filament selection in the local thickness image
		saveAs("Tiff", nameResultDir + "skeleton_locThk");
		close("\\Others");
		run("Analyze Skeleton (2D/3D)", "prune=none show display");
		
    	// open the new result file containing the related branches results
		fileResultsBranches = File.open(nameResultDir+ "results_Branches.csv");
		print(fileResultsBranches, "Branches_Nb,Count,Branches_length,ID,endPoint_x,endPoint_y,junctionPoint_x,junctionPoint_y,Intensity,Junctions_Nb\n");
		// save the mask-labeled-skeleton image
		saveAs("Tiff", nameResultDir + "mask_labeled-skeleton");
		selectImage("Tagged skeleton");
		// save the mask tagged skeleton 
		saveAs("Tiff", nameResultDir + "mask_tagged-skeleton");
		
		prev_Branches = 0;
		selectWindow("Branch information");
		// do a loop going throught and saving all the results corresponding to the banche informations
		for (n = 0; n < nResults; n++) {
			branches = getResult("# Branches", n);
			branchingPoints = getResult("# Junctions", n);
			// write informations about branches and branchingPoints for each skeleton found by the program 
			print(fileResultsBranches, branches+",,,,,,,,,"+branchingPoints+"\n");
			selectWindow("Branch information");
			// loop througth the result table and get some information about end point position and corresponding branching point position
			for (m = prev_Branches; m < (branches+prev_Branches); m++) {
				Branch_Length = Table.get("Branch length", m);
				skeleton_ID = Table.get("Skeleton ID", m);
				x1 = Table.get("V1 x", m);
				y1 = Table.get("V1 y", m);
				x2 = Table.get("V2 x", m);
				y2 = Table.get("V2 y", m);
				intensity = Table.get("average intensity (inner 3rd)", m);
				if (Branch_Length > 0) {
					// prevention loop, the number can be modifie to only register the information about branches > the length you enter
					print(fileResultsBranches," ,"+(m-prev_Branches)+","+Branch_Length+","+skeleton_ID+","+x1+","+y1+","+x2+","+y2+","+intensity+"\n");
				}
				selectWindow("Branch information");
			}	 
			prev_Branches = (prev_Branches+branches);
		}
		File.close(fileResultsBranches);
		close("*");
    }
}

