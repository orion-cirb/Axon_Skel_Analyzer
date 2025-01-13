/*
 * Description: Segment nuclei and filaments, perform skeleton and local thickness analysis of filaments
 * Developed for: Sandra, Fuchs' team
 * Author: Thomas Caille & Héloïse Monnet @ ORION-CIRB 
 * Date: January 2025
 * Repository: https://github.com/orion-cirb/Axon_Skel_Analyzer
 * Dependencies: None
*/



// Hide images during macro execution
setBatchMode(true);

// Ask for the images directory
inputDir = getDirectory("Please select a directory containing images to analyze");
print("Analysis started");

// Create results directory
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
resultDir = inputDir + "Results_" + year + "-" + (month+1) + "-" + dayOfMonth + "_" + hour + "-" + minute + "-" + second + File.separator();
if (!File.isDirectory(resultDir)) {
	File.makeDirectory(resultDir);
}

// Get all files in the input directory
inputFiles = getFileList(inputDir);

// Create global/branches_length/branches_diam results files and write headers in them
fileResultsGlobal = File.open(resultDir + "results_global.csv");
print(fileResultsGlobal, "Image name,Nuclei nb,Filaments total area (µm2),Filaments branches nb,Filaments branches total length (µm),Filaments branches mean diam (µm),Filaments junctions nb\n");
File.close(fileResultsGlobal);
fileResultsBranchesLength = File.open(resultDir + "results_branches_length.csv");
print(fileResultsBranchesLength, "Image name,Branch length (µm),V1 x,V1 y,V2 x,V2 y\n");
File.close(fileResultsBranchesLength);
fileResultsBranchesDiam = File.open(resultDir + "results_branches_diam.csv");
print(fileResultsBranchesDiam, "Image name,Branch mean diam (µm)\n");
File.close(fileResultsBranchesDiam);

// Loop through all files with .nd extension
for (i = 0; i < inputFiles.length; i++) {
    if (endsWith(inputFiles[i], ".nd")) {
    	print("Analyzing image " + inputFiles[i] + "...");
    	imgName = replace(inputFiles[i],".nd","");
		
		// Open DAPI nuclei channel (channel 1)
		run("Bio-Formats Importer", "open=["+inputDir + inputFiles[i]+"] autoscale color_mode=Default specify_range split_channels view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
    	// Perform max intensity z-projection
    	run("Z Project...", "projection=[Max Intensity]");
    	
    	// Detect and count nuclei
    	// Preprocessing
    	getStatistics(area,mean,min,max,std,histogram);
		run("Subtract...", "value=" + (mean + std/2));
		run("Median...", "radius=10");
		// Huang thresholding
		setAutoThreshold("Huang dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		// Postprocessing
		run("Fill Holes");
		run("Watershed");
		// Filter out nuclei with area < 40 µm2 and circularity < 0.7
		run("Set Measurements...", "redirect=None decimal=2");
    	run("Analyze Particles...", "size=40-Infinity circularity=0.70-1.00 clear add");
    	nbNuclei = nResults;
    	// Save ROIs
    	roiManager("Save", resultDir + imgName + "_nuclei.zip");
    	roiManager("reset");
    	close("*");
		
    	// Open Tug filaments channel (channel 2)
    	run("Bio-Formats Importer", "open=["+inputDir + inputFiles[i]+"] autoscale color_mode=Default specify_range split_channels view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");
    	// Perform max intensity z-projection
    	run("Z Project...", "projection=[Max Intensity]");
		
		// Segment filaments
		// Preprocessing
		run("Subtract Background...", "rolling=50 sliding");
		run("Median...", "radius=2");
  		// Huang thresholding
		setAutoThreshold("Huang dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		// Postprocessing
		run("Options...", "iterations=4 count=1 black pad do=Close");
		run("Median...", "radius=2");
		// Fill holes with area < 10 µm2
		fillHoles(0, 10);
		rename("filamentsMask");
		close("\\Others");
		
		// Open ORF1p cell bodies channel (channel 3)
    	run("Bio-Formats Importer", "open=["+inputDir + inputFiles[i]+"] autoscale color_mode=Default specify_range split_channels view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");
    	nbSlices = nSlices;
		// Perform sum slices z-projection
    	run("Z Project...", "projection=[Sum Slices]");
    	
    	// Segment cell bodies
		// Preprocessing
		run("Median...", "radius=50");
		// Huang thresholding
		setAutoThreshold("Huang dark");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		// Postprocessing
		run("Options...", "iterations=20 count=1 black do=Open");
		run("Options...", "iterations=5 count=1 black do=Dilate");
		run("Fill Holes");
		// Filter out cell bodies with area < 200 µm2
		run("Analyze Particles...", "size=200-Infinity show=Masks");
		run("Invert LUT");
	
		// Clear cell bodies in filaments mask
		run("Create Selection");
		selectImage("filamentsMask");
		run("Restore Selection");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
		run("Select None");
		close("\\Others");
		
		// Filter out filaments with area < 20 µm2
		run("Analyze Particles...", "size=20-Infinity show=Masks");
		run("Invert LUT");
		// Save filaments mask
		saveAs("Tiff", resultDir + imgName + "_filaments");
		rename("filamentsMask");
		close("\\Others");
		// Measure filaments total area
		run("Create Selection");
		run("Set Measurements...", "area redirect=None decimal=2");
		List.setMeasurements;
		areaFilaments = List.getValue("Area");
		run("Select None");
    
		// Skeletonize filaments
		selectImage("filamentsMask");
		run("Duplicate...", "title=filamentsSkel");
		run("Skeletonize");
		// Filter out branches with area < 0.4 µm2 (length < 4 µm approximately)
		filterBranches(0.4);
		saveAs("Tiff", resultDir + imgName + "_filaments_skel");
		rename("filamentsSkel");
		
		// Analyze skeleton
		run("Analyze Skeleton (2D/3D)", "prune=none show");
		nbBranches = 0;
		nbJunctions = 0;
		for (n = 0; n < nResults; n++) {
    		nbBranches += getResult("# Branches", n);
    		nbJunctions += getResult("# Junctions", n);
    	}
		lengthBranches = 0;
		for (n = 0; n < Table.size; n++) {
			selectWindow("Branch information");
    		lengthBranches += Table.get("Branch length", n);
    	}
    	
		// Save parameters in branches_length results file
		selectWindow("Branch information");
		for (n = 0; n < Table.size; n++) {
			File.append(imgName+","+Table.get("Branch length", n)+","+Table.get("V1 x", n)+","+Table.get("V1 y", n)+","+Table.get("V2 x", n)+","+Table.get("V2 y", n), resultDir+"results_branches_length.csv");
		}
		close("Results");
		close("Branch information");
		close("Tagged skeleton");
		
		// Compute filaments local thickness
		selectImage("filamentsMask");
		run("Local Thickness (masked, calibrated, silent)");
		setMinAndMax(0, 20);
		saveAs("Tiff", resultDir + imgName + "_filaments_locThk");
		rename("filamentsLocThk");
		
		// Get skeleton mean intensity on the local thickness image
		selectImage("filamentsSkel");
		run("Create Selection");
		selectImage("filamentsLocThk");
		run("Restore Selection");
		run("Set Measurements...", "mean redirect=None decimal=2");
		List.setMeasurements;
		diamBranches = List.getValue("Mean");
		run("Select None");
		
		// Save parameters in global results file
    	File.append(imgName+","+nbNuclei+","+areaFilaments+","+nbBranches+","+lengthBranches+","+diamBranches+","+nbJunctions, resultDir+"results_global.csv");
		
		// Get skeleton without junctions and endpoints
		selectImage("filamentsSkel");
		run("Select None");
		getSkeletonWithoutJunctions(false);
		
		// Save parameters in branches_diam results file
		run("Set Measurements...", "mean redirect=[filamentsLocThk] decimal=2");
		run("Analyze Particles...", "size=0-Infinity show=Nothing display clear add");
		for (n = 0; n < nResults; n++) {
			File.append(imgName+","+getResult("Mean", n), resultDir+"results_branches_diam.csv");
		}	
		
		close("*");
		close("Results");
		roiManager("reset");
    }
}

setBatchMode(false);

print("Analysis done!");



/******************** UTILS ********************/

function fillHoles(minHoleArea, maxHoleArea) {
	run("Invert");
	run("Analyze Particles...", "size=" + minHoleArea + "-" + maxHoleArea + " add");
	run("Invert");
	roiManager("deselect");
	roiManager("combine");
	setForegroundColor(255, 255, 255);
	run("Fill", "slice");
	roiManager("reset");
	run("Select None");
}


function filterBranches(minBranchArea) {
	skelImgName = getTitle();
	
	getSkeletonWithoutJunctions(true);
	run("Set Measurements...", "area min redirect=[Tagged skeleton] decimal=2");
	run("Analyze Particles...", "size=0-Infinity show=Nothing display clear add");
	
	selectWindow(skelImgName);
	setBackgroundColor(0, 0, 0);
	for (n = 0; n < nResults; n++) {		
		if(getResult("Area", n) < minBranchArea && getResult("Min", n) == 30) {
			roiManager("select", n);
			run("Clear", "slice");
		}
	}
	run("Skeletonize");
	
	close(skelImgName+"NoJct");
	close("Tagged skeleton");
	close("Results");
	roiManager("reset");
}


function getSkeletonWithoutJunctions(keepTaggedSkel) {
	outputImgName = getTitle()+"NoJct";
	run("Duplicate...", "title="+outputImgName);
	
	run("Analyze Skeleton (2D/3D)", "prune=none");
	run("Duplicate...", "title=jct");
	setThreshold(60, 80, "raw");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	imageCalculator("Subtract", outputImgName, "jct");
	
	selectWindow(outputImgName);
	close("jct");
	close("Results");
	if(!keepTaggedSkel) {
		close("Tagged skeleton");
	}
}

/***********************************************/
