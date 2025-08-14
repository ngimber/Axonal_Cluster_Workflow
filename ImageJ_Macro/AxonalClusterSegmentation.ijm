//******************************
//***** Enter parameters here **
//******************************


// File type extension for input images
extention=".tif"

// Channel index used to derive the axonal skeleton
skeletonChannel=1;
// Enable to pause at thresholding for manual review
debugging=false;

//**************** prepare ImageJ for analysis ****************
dir=getDirectory("Choose Directory");
files=getFileList(dir);
subfolder=dir+"\\analysis_IJ\\";
File.makeDirectory(subfolder);
File.makeDirectory(subfolder+"ROIs\\");
File.makeDirectory(subfolder+"Overlays\\");
File.makeDirectory(subfolder+"Skeletons\\");
File.makeDirectory(subfolder+"Binaries\\");
File.makeDirectory(subfolder+"Tables_bgsubtracted\\");
File.makeDirectory(subfolder+"Tables\\");
File.makeDirectory(subfolder+"AxonLength\\");


// ***segmentation parameters***

// Size filters for Analyze Particles (in pixels unless calibrated)
minParticleSize=0;
maxParticleSize=999999;
crop=false;// true restricts segmentation to the skeleton area of channel 1
run("Options...", "iterations=1 count=1 black do=Nothing");

// ***Preparation for the parameter dialog***
stop=false;
i=0;
close("*");
while (stop==false)
	{
	i++;
	if ((indexOf(files[i], extention)) >= 0) 
		{	
		run("Bio-Formats", "  open="+dir+files[i]+" color_mode=Default view=Hyperstack stack_order=XYCZT");
		run("Grays");
		getDimensions(width, height,channels , slices, frames);
		print(channels);
		stop=true;
		}
	}
	close("*");
Threshs=newArray(channels-1);
// Available automatic thresholding methods (dark foreground variants)
methods=newArray("Default Dark", "Huang Dark", "Intermodes Dark", "IsoData Dark", "IJ_IsoData Dark", "Li Dark", "MaxEntropy Dark", "Mean Dark", "MinError Dark", "Minimum Dark", "Moments Dark", "Otsu Dark", "Percentile Dark", "RenyiEntropy Dark", "Shanbhag Dark", "Triangle Dark", "Yen Dark");
Thresh_default= "Moments Dark";
Balls=newArray(channels-1);
Balls_default=3;
Sigmas=newArray(channels-1);
Sigmas_default=2;
Opens=newArray(channels-1);
Opens_default=0;
chosenFactors=newArray(channels);
factor_default=2;
// Initialize per-channel defaults
for (i=0;i<channels;i++)
		{
		Threshs[i]=Thresh_default;
		chosenFactors[i]=factor_default;
		Balls[i]=Balls_default;
		Sigmas[i]=Sigmas_default;
		Opens[i]=Opens_default;
		}
// Parameter dialog
Dialog.create("Axonal Cluster Segmentation");
  	
  	Dialog.addMessage("File Type");
  	Dialog.addString("extention", extention);
  	Dialog.addMessage("\n");
  	Dialog.addMessage("\n");	
  	Dialog.addNumber("Channel for Axon Segmentation (1-"+channels+")",skeletonChannel);
  	Dialog.addMessage("\n");
  	Dialog.addMessage("\n");	  	
  	Dialog.addCheckbox("Debug Mode (stops during thresholding)",debugging);
  	Dialog.addMessage("\n");
  	Dialog.addMessage("\n");	
  
	print(channels);
  	for (i=1;i<=channels;i++)
		{
		Dialog.addMessage("\n");
		Dialog.addMessage("\n");
		Dialog.addNumber("Ballsize channel "+i+1+":  ", Balls[i-1]);
		Dialog.addNumber("sigma channel "+i+1+": ", Sigmas[i-1]);
		Dialog.addNumber("opening channel "+i+1+": ", Opens[i-1]);
		// Choose auto-threshold method per channel; optional skip
		//Dialog.addString("Threshold channel "+i+1+": Lower Threshold (0 for no binarization):  ", Threshs[i-1]);
		Dialog.addChoice("Threshold channel "+i+1+": Lower Threshold", methods, Thresh_default);
		Dialog.addCheckbox("Skip Binarization", false);
		Dialog.addNumber("Decrease Auto-Threshold by Factor ", chosenFactors[i-1]);
		}
 


  	Dialog.show();


	
	extention = Dialog.getString(); 
	skeletonChannel=Dialog.getNumber();
	debugging=Dialog.getCheckbox();
	
	
	
	
  	for (i=1;i<=channels;i++)
		{
		Balls[i-1]=Dialog.getNumber();
		Sigmas[i-1]=Dialog.getNumber();
		Opens[i-1]=Dialog.getNumber();
		//Threshs[i-1]=Dialog.getString();
		Threshs[i-1]=Dialog.getChoice();
		skip=Dialog.getCheckbox();
		if(skip==true){Threshs[i-1]=0;}
		chosenFactors[i-1]=Dialog.getNumber();
		}




// Generate a log file with all chosen parameters
print("\\Clear");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print("\n");
print(year+"-"+month+"-"+dayOfWeek+": "+hour+"-"+minute+"-"+second);
print("extention "+extention);  
print("axon channel "+skeletonChannel);


 	for (i=1;i<=channels;i++)
		{
		print("\n");
		print("channel", i);
		
		print("ballsize",Balls[i-1]);
		print("sigma",Sigmas[i-1]);
		print("opening",Opens[i-1]);
		print("threshold method", Threshs[i-1]);
		print("factor",chosenFactors[i-1]);
		}
selectWindow("Log");
saveAs("Text", ""+subfolder+"log.txt");

print("\\Clear");
roiManager("reset");

// Set measurement outputs and CSV export behavior
run("Set Measurements...", "area mean standard modal min centroid center perimeter fit shape feret's integrated median display redirect=None decimal=9");
run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file copy_row save_column save_row");

//**************** Main Loop
for (i=0;i<files.length;i++)
	{
	if ((indexOf(files[i], extention)) >= 0) 

		{
		if (File.exists(subfolder+"Binaries\\Dapi_"+files[i]+".tif")==true)
	{print(files[i]+" already exists");}
		else 

			{
			print(files[i]+" does not exist");
		
			while (nImages>0) 
			{ 
				  selectImage(nImages); 
				  close(); 
			  	} 

				if(isOpen("Results")==true)
				{
			  	close("Results"); 
				}
				
			run("Collect Garbage");
			call("java.lang.System.gc");
	
			run("Bio-Formats", "  open="+dir+files[i]+" color_mode=Default view=Hyperstack stack_order=XYCZT");
			title=getTitle();
	
	
			
			//**********************************************
			//********** Segmentations and Analysis ********
			//**********************************************


			//********** Get Skeleton **********
				roiManager("reset");
				run("Select All");	
				getPixelSize(unit, pixelWidth, pixelHeight);
				getDimensions(width, height, channels, slices, frames);
				setSlice(skeletonChannel);

				run("Duplicate...", "title=tmp channels=1");
				run("Normalize Local Contrast", "block_radius_x=20 block_radius_y=20 standard_deviations=10 center");
				run("Gaussian Blur...", "sigma=2");
				
				setThreshold(10, 9999999, "raw");
				run("Create Selection");
				
				setAutoThreshold("Yen");
				run("Convert to Mask");
				run("Invert");
			
				run("Morphological Filters", "operation=Opening element=Disk radius=2");
				run("Morphological Filters", "operation=Closing element=Disk radius=6");
				run("Divide...", "value=255");
				rename("Skeleton");		
				run("Clear Results");
				roiManager("reset");
				run("Select All");
				run("Measure");
				print(width);

				axonLength=getResult("Mean", 0)*width*height*pixelWidth;// this is correct: skeleton is one pixel wide
				print("\\Clear");
				print("axon Length in "+unit+": "+axonLength);
				selectWindow("Log");
				saveAs("Text", ""+subfolder+"AxonLength\\"+title+"_axonLength.txt");
				print("\\Clear");
				selectWindow("Skeleton");
				run("Multiply...", "value=255");
				saveAs("Tiff", subfolder+"Skeletons\\Skeleton_"+title+".tif");
				rename("Skeleton");

			//********** Segmentation **********
				selectWindow(title);
				getDimensions(width, height, channels, slices, frames);
				for (c=0;c<channels;c++)
				{	
				selectWindow(title);
				run("Select All");
				setSlice(c+1);
				run("Duplicate...", "  channels="+c+1+"");
				wait(100);
				run("Subtract Background...", "rolling="+Balls[c]+" stack");
				run("Gaussian Blur...", "sigma="+Sigmas[c]+"");
				run("Morphological Filters", "operation=Opening element=Square radius="+Opens[c]+"");

				run("Multiply...", "value=10");
				resetThreshold;
				setAutoThreshold(Threshs[c]);
				getThreshold(lower, upper);
				setThreshold(round(upper*chosenFactors[c]), 99999999999999999);

		

	
				if(debugging==true){waitForUser ("adjust Threshold");}
				

				run("Convert to Mask");
				//run("Invert");
				run("Watershed");
				rename("binary");
				run("Select All");
				roiManager("reset");
				selectWindow("binary");
				run("Analyze Particles...", "display clear add");
				run("Analyze Particles...", "size="+minParticleSize+"-"+maxParticleSize+" pixel display clear add");
				selectWindow("binary");
				saveAs("Tiff", subfolder+"Binaries\\"+c+"_"+title+".tif");
				close(""+c+"_"+title+".tif");
				roiManager("Save", subfolder+"ROIs\\"+c+"_"+title+".zip");
				selectImage(title);

				
				
				
				
				
			//********** Optional: crop to axonal area only ********** >>> not tested
				if(crop==true)
				
				{
				selectImage("Skeleton");
				run("Gaussian Blur...", "sigma=2 scaled");
				run("Clear Results");
				for (ri = roiManager("count")-1; ri < 0; ri--) 
					{				
					roiManager("Select", ri);
					run("Clear Results");
					run("Measure");
					if (getResult("IntDen", 0)<=0);
							{
							roiManager("Select", ri);
							roiManager("delete");
							print("deleted"+ri);
							}
					
					}
				}

				
			//********** Measurements **********				
				
				
				for (ch=1;ch<=channels-1;ch++)
				{
				selectImage(title);	
				setSlice(ch);
				run("Clear Results");
				for (ri = 0; ri < roiManager("count")-1; ri++) 
					{				
					roiManager("Select", ri);
					roiManager("Select", ri);
					run("Measure");
					}
					saveAs("Results", ""+subfolder+"Tables\\"+title+"segmentCh"+c+1+"measureCh"+ch+".csv");
					close(title+"segmentCh"+c+"measureCh"+ch+".csv");					
					
					
				run("Clear Results");
				
				selectImage(title);	
				setSlice(ch);
				run("Select All");
				run("Measure");
				mode=getResult("Mode", 0);
				run("Clear Results");		
				run("Subtract...", "value="+mode+" slice");
				for (ri = 0; ri < roiManager("count")-1; ri++) 
					{				
					roiManager("Select", ri);
					roiManager("Select", ri);
					run("Measure");
					}
					saveAs("Results", ""+subfolder+"Tables_bgSubtracted\\noBG"+title+"segmentCh"+c+1+"+measureCh"+ch+".csv");
					close("noBG"+title+"segmentCh"+c+1+"+measureCh"+ch+".csv");
					run("Clear Results");
					
				}

				
			//********** Create Overlay **********		
			selectImage(title);	
			run("Select All");
			setSlice(c+1);			
			run("Duplicate...", "  channels="+c+1+"");	


			roiManager("Deselect");
			roiManager("Combine");
			setForegroundColor(255, 255, 255);
			setLineWidth(1);
			run("Draw", "slice");	
			saveAs("Tiff", subfolder+"Overlays\\"+c+"_"+title+".tif");				


				}
			}
		}}
		
		close("*");

				
					
	


// Utility: read the current Results table column into an array
function readArray(columnname)
	{
	storageArray=newArray(nResults);
	for(row=0;row<nResults;row++)
		{
		storageArray[row]=getResult(columnname, row);
		}
		return storageArray;
	}
