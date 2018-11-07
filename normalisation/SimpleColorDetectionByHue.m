% Demo macro to very, very simple color detection in 
% HSV (Hue, Saturation, Value) color space.
% Requires the Image Processing Toolbox.  Developed under MATLAB R2014a.
%
function SimpleColorDetectionByHue()
clc;	% Clear command window.
clear;	% Delete all variables.
close all;	% Close all figure windows except those created by imtool.
imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.

% Change the current folder to the folder of this m-file.
% (The "cd" line of code below is from Brett Shoelson of The Mathworks.)
if(~isdeployed)
	cd(fileparts(which(mfilename))); % From Brett
end

ver % Display user's toolboxes in their command window.

% Introduce the demo, and ask user if they want to continue or exit.
message = sprintf('This demo will illustrate very simple color detection\nin HSV color space.\nIt requires the Image Processing Toolbox.\nDo you wish to continue?');
reply = questdlg(message, 'Run Demo?', 'OK','Cancel', 'OK');
if strcmpi(reply, 'Cancel')
	% User canceled so exit.
	return;
end

try
	% Check that user has the Image Processing Toolbox installed.
	hasIPT = license('test', 'image_toolbox');
	if ~hasIPT
		% User does not have the toolbox installed.
		message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
		reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
		if strcmpi(reply, 'No')
			% User said No, so exit.
			return;
		end
	end

	% Continue with the demo.  Do some initialization stuff.
	close all;
	fontSize = 16;
	figure;
	% Maximize the figure. 
	set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 

	% Change the current folder to the folder of this m-file.
	% (The line of code below is from Brett Shoelson of The Mathworks.)
	if(~isdeployed)
		cd(fileparts(which(mfilename)));
	end

	% Ask user if they want to use a demo image or their own image.
	message = sprintf('Do you want use a standard demo image,\nOr pick one of your own?');
	reply2 = questdlg(message, 'Which Image?', 'Demo','My Own', 'Demo');
	% Open an image.
	if strcmpi(reply2, 'Demo')
		% Read standard MATLAB demo image.
	% 	fullImageFileName = 'peppers.png';
		message = sprintf('Which demo image do you want to use?');
		selectedImage = questdlg(message, 'Which Demo Image?', 'Onions', 'Peppers', 'Kids', 'Onions');
		if strcmp(selectedImage, 'Onions')
			fullImageFileName = 'onion.png';
		elseif strcmp(selectedImage, 'Peppers')
			fullImageFileName = 'peppers.png';
		else
			fullImageFileName = 'kids.tif';
		end
	else
		% They want to pick their own.
		% Change default directory to the one containing the standard demo images for the MATLAB Image Processing Toolbox. 
		originalFolder = pwd; 
		folder = 'C:\Program Files\MATLAB\R2010a\toolbox\images\imdemos'; 
		if ~exist(folder, 'dir') 
			folder = pwd; 
		end 
		cd(folder); 
		% Browse for the image file. 
		[baseFileName, folder] = uigetfile('*.*', 'Specify an image file'); 
		fullImageFileName = fullfile(folder, baseFileName); 
		% Set current folder back to the original one. 
		cd(originalFolder);
		selectedImage = 'My own image'; % Need for the if threshold selection statement later.

	end

	% Check to see that the image exists.  (Mainly to check on the demo images.)
	if ~exist(fullImageFileName, 'file')
		message = sprintf('This file does not exist:\n%s', fullImageFileName);
		uiwait(msgbox(message));
		return;
	end

	% Read in image into an array.
	[rgbImage, storedColorMap] = imread(fullImageFileName); 
	[rows, columns, numberOfColorBands] = size(rgbImage); 
	% If it's monochrome (indexed), convert it to color. 
	% Check to see if it's an 8-bit image needed later for scaling).
	if strcmpi(class(rgbImage), 'uint8')
		% Flag for 256 gray levels.
		eightBit = true;
	else
		eightBit = false;
	end
	if numberOfColorBands == 1
		if isempty(storedColorMap)
			% Just a simple gray level image, not indexed with a stored color map.
			% Create a 3D true color image where we copy the monochrome image into all 3 (R, G, & B) color planes.
			rgbImage = cat(3, rgbImage, rgbImage, rgbImage);
		else
			% It's an indexed image.
			rgbImage = ind2rgb(rgbImage, storedColorMap);
			% ind2rgb() will convert it to double and normalize it to the range 0-1.
			% Convert back to uint8 in the range 0-255, if needed.
			if eightBit
				rgbImage = uint8(255 * rgbImage);
			end
		end
	end 
	
	% Display the original image.
	subplot(3, 4, 1);
	hRGB = imshow(rgbImage);
	% Set up an infor panel so you can mouse around and inspect the value values.
	hrgbPI = impixelinfo(hRGB);
	set(hrgbPI, 'Units', 'Normalized', 'Position',[.15 .69 .15 .02]);
	drawnow; % Make it display immediately. 
	if numberOfColorBands > 1 
		title('Original Color Image', 'FontSize', fontSize); 
	else 
		caption = sprintf('Original Indexed Image\n(converted to true color with its stored colormap)');
		title(caption, 'FontSize', fontSize);
	end

	% Convert RGB image to HSV
	hsvImage = rgb2hsv(rgbImage);
	% Extract out the H, S, and V images individually
	hImage = hsvImage(:,:,1);
	sImage = hsvImage(:,:,2);
	vImage = hsvImage(:,:,3);
	
	% Display the hue image.
	subplot(3, 4, 2);
	h1 = imshow(hImage);
	title('Hue Image', 'FontSize', fontSize);
	% Set up an infor panel so you can mouse around and inspect the hue values.
	hHuePI = impixelinfo(h1);
	set(hHuePI, 'Units', 'Normalized', 'Position',[.34 .69 .15 .02]);
	
	% Display the saturation image.
	h2 = subplot(3, 4, 3);
	imshow(sImage);
	title('Saturation Image', 'FontSize', fontSize);
	% Set up an infor panel so you can mouse around and inspect the saturation values.
	hSatPI = impixelinfo(h2);
	set(hSatPI, 'Units', 'Normalized', 'Position',[.54 .69 .15 .02]);
	
	% Display the value image.
	h3 = subplot(3, 4, 4);
	imshow(vImage);
	title('Value Image', 'FontSize', fontSize);
	% Set up an infor panel so you can mouse around and inspect the value values.
	hValuePI = impixelinfo(h3);
	set(hValuePI, 'Units', 'Normalized', 'Position',[.75 .69 .15 .02]);

	message = sprintf('These are the individual HSV color bands.\nNow we will compute the image histograms.');
	reply = questdlg(message, 'Continue with Demo?', 'OK','Cancel', 'OK');
	if strcmpi(reply, 'Cancel')
		% User canceled so exit.
		return;
	end

	% Compute and plot the histogram of the "hue" band.
	hHuePlot = subplot(3, 4, 6); 
	[hueCounts, hueBinValues] = imhist(hImage); 
	maxHueBinValue = find(hueCounts > 0, 1, 'last'); 
	maxCountHue = max(hueCounts); 
	% Display with area() rather than bar, due to bug in bar().  The bug, and workaround of using area(), are discussed in
	% http://www.mathworks.com/matlabcentral/answers/103538-why-does-a-bar-subplot-change-when-i-create-another-bar-subplot-on-the-same-figure-in-matlab-8-0-r
	% Supposedly it's been fixed in R2014b.
% 	bar(hueBinValues, hueCounts, 'r'); 
	area(hueBinValues, hueCounts, 'FaceColor', 'r'); 
	grid on; 
	xlabel('Hue Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Hue Image', 'FontSize', fontSize);

	% Compute and plot the histogram of the "saturation" band.
	hSaturationPlot = subplot(3, 4, 7); 
	[saturationCounts, saturationBinValues] = imhist(sImage); 
	maxSaturationBinValue = find(saturationCounts > 0, 1, 'last'); 
	maxCountSaturation = max(saturationCounts); 
% 	bar(saturationBinValues, saturationCounts, 'g', 'BarWidth', 0.95); 
	area(saturationBinValues, saturationCounts, 'FaceColor', 'g'); 
	grid on; 
	xlabel('Saturation Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Saturation Image', 'FontSize', fontSize);

	% Compute and plot the histogram of the "value" band.
	hValuePlot = subplot(3, 4, 8); 
	[valueCounts, valueBinValues] = imhist(vImage); 
	maxValueBinValue = find(valueCounts > 0, 1, 'last'); 
	maxCountValue = max(valueCounts); 
% 	bar(valueBinValues, valueCounts, 'b'); 
	area(valueBinValues, valueCounts, 'FaceColor', 'b'); 
	grid on; 
	xlabel('Value Value'); 
	ylabel('Pixel Count'); 
	title('Histogram of Value Image', 'FontSize', fontSize);

	% Set all axes to be the same width and height.
	% This makes it easier to compare them.
	maxCount = max([maxCountHue,  maxCountSaturation, maxCountValue]); 
	axis([hHuePlot hSaturationPlot hValuePlot], [0 1 0 maxCount]); 

	% Plot all 3 histograms in one plot.
	subplot(3, 4, 5); 
	plot(hueBinValues, hueCounts, 'r', 'LineWidth', 2); 
	grid on; 
	xlabel('Values'); 
	ylabel('Pixel Count'); 
	hold on; 
	plot(saturationBinValues, saturationCounts, 'g', 'LineWidth', 2); 
	plot(valueBinValues, valueCounts, 'b', 'LineWidth', 2); 
	title('Histogram of All Bands', 'FontSize', fontSize); 
	maxGrayLevel = max([maxHueBinValue, maxSaturationBinValue, maxValueBinValue]); % Just for our information....
	% Make x-axis to just the max gray level on the bright end. 
	xlim([0 1]); 

	% Now select thresholds for the 3 color bands.
	message = sprintf('Now we will select some color threshold ranges\nand display them over the histograms.');
	reply = questdlg(message, 'Continue with Demo?', 'OK','Cancel', 'OK');
	if strcmpi(reply, 'Cancel')
		% User canceled so exit.
		return;
	end

	% Assign the low and high thresholds for each color band.
	if strcmpi(reply2, 'My Own') || strcmpi(selectedImage, 'Kids') > 0
		% Take a guess at the values that might work for the user's image.
		hueThresholdLow = 0;
		hueThresholdHigh = graythresh(hImage);
		saturationThresholdLow = graythresh(sImage);
		saturationThresholdHigh = 1.0;
		valueThresholdLow = graythresh(vImage);
		valueThresholdHigh = 1.0;
	else
		% Use values that I know work for the onions and peppers demo images.
		[hueThresholdLow, hueThresholdHigh, saturationThresholdLow, saturationThresholdHigh, valueThresholdLow, valueThresholdHigh] = SetThresholds()
	end
	
	% Interactively and visually set/adjust thresholds using custom thresholding application.
	% Available on the File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/29372-thresholding-an-image
% 	[hueThresholdLow, hueThresholdHigh] = threshold(hueThresholdLow, hueThresholdHigh, hImage);
% 	[saturationThresholdLow, saturationThresholdHigh] = threshold(saturationThresholdLow, saturationThresholdHigh, sImage);
% 	[valueThresholdLow, valueThresholdHigh] = threshold(valueThresholdLow, valueThresholdHigh, vImage);

	% Show the thresholds as vertical magenta bars on the histograms.
	PlaceThresholdBars(6, hueThresholdLow, hueThresholdHigh);
	PlaceThresholdBars(7, saturationThresholdLow, saturationThresholdHigh);
	PlaceThresholdBars(8, valueThresholdLow, valueThresholdHigh);

	message = sprintf('Next we will apply each color band threshold range to its respective color band.');
	reply = questdlg(message, 'Continue with Demo?', 'OK','Cancel', 'OK');
	if strcmpi(reply, 'Cancel')
		% User canceled so exit.
		return;
	end

	% Now apply each color band's particular thresholds to the color band
	hueMask = (hImage >= hueThresholdLow) & (hImage <= hueThresholdHigh);
	saturationMask = (sImage >= saturationThresholdLow) & (sImage <= saturationThresholdHigh);
	valueMask = (vImage >= valueThresholdLow) & (vImage <= valueThresholdHigh);

	% Display the thresholded binary images.
	fontSize = 16;
	subplot(3, 4, 10);
	imshow(hueMask, []);
	title('=   Hue Mask', 'FontSize', fontSize);
	subplot(3, 4, 11);
	imshow(saturationMask, []);
	title('&   Saturation Mask', 'FontSize', fontSize);
	subplot(3, 4, 12);
	imshow(valueMask, []);
	title('&   Value Mask', 'FontSize', fontSize);
	% Combine the masks to find where all 3 are "true."
	% Then we will have the mask of only the red parts of the image.
	coloredObjectsMask = uint8(hueMask & saturationMask & valueMask);
	subplot(3, 4, 9);
	imshow(coloredObjectsMask, []);
	caption = sprintf('Mask of Only Regions\nof The Specified Color');
	title(caption, 'FontSize', fontSize);

	% Tell user that we're going to filter out small objects.
	smallestAcceptableArea = 100; % Keep areas only if they're bigger than this.
	message = sprintf('Note the small regions in the image in the lower left.\nNext we will eliminate regions smaller than %d pixels.', smallestAcceptableArea);
	reply = questdlg(message, 'Continue with Demo?', 'OK','Cancel', 'OK');
	if strcmpi(reply, 'Cancel')
		% User canceled so exit.
		return;
	end

	% Open up a new figure, since the existing one is full.
	figure;  
	% Maximize the figure. 
	set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

	% Get rid of small objects.  Note: bwareaopen returns a logical.
	coloredObjectsMask = uint8(bwareaopen(coloredObjectsMask, smallestAcceptableArea));
	subplot(3, 3, 1);
	imshow(coloredObjectsMask, []);
	fontSize = 13;
	caption = sprintf('bwareaopen() removed objects\nsmaller than %d pixels', smallestAcceptableArea);
	title(caption, 'FontSize', fontSize);

	% Smooth the border using a morphological closing operation, imclose().
	structuringElement = strel('disk', 4);
	coloredObjectsMask = imclose(coloredObjectsMask, structuringElement);
	subplot(3, 3, 2);
	imshow(coloredObjectsMask, []);
	fontSize = 16;
	title('Border smoothed', 'FontSize', fontSize);

	% Fill in any holes in the regions, since they are most likely red also.
	coloredObjectsMask = imfill(logical(coloredObjectsMask), 'holes');
	subplot(3, 3, 3);
	imshow(coloredObjectsMask, []);
	title('Regions Filled', 'FontSize', fontSize);

	message = sprintf('This is the filled, size-filtered mask.\nNext we will apply this mask to the original RGB image.');
	reply = questdlg(message, 'Continue with Demo?', 'OK','Cancel', 'OK');
	if strcmpi(reply, 'Cancel')
		% User canceled so exit.
		return;
	end

	% You can only multiply integers if they are of the same type.
	% (coloredObjectsMask is a logical array.)
	% We need to convert the type of coloredObjectsMask to the same data type as hImage.
	coloredObjectsMask = cast(coloredObjectsMask, 'like', rgbImage); 
% 	coloredObjectsMask = cast(coloredObjectsMask, class(rgbImage));

	% Use the colored object mask to mask out the colored-only portions of the rgb image.
	maskedImageR = coloredObjectsMask .* rgbImage(:,:,1);
	maskedImageG = coloredObjectsMask .* rgbImage(:,:,2);
	maskedImageB = coloredObjectsMask .* rgbImage(:,:,3);
	% Show the masked off red image.
	subplot(3, 3, 4);
	imshow(maskedImageR);
	title('Masked Red Image', 'FontSize', fontSize);
	% Show the masked off saturation image.
	subplot(3, 3, 5);
	imshow(maskedImageG);
	title('Masked Green Image', 'FontSize', fontSize);
	% Show the masked off value image.
	subplot(3, 3, 6);
	imshow(maskedImageB);
	title('Masked Blue Image', 'FontSize', fontSize);
	% Concatenate the masked color bands to form the rgb image.
	maskedRGBImage = cat(3, maskedImageR, maskedImageG, maskedImageB);
	% Show the masked off, original image.
	subplot(3, 3, 8);
	imshow(maskedRGBImage);
	fontSize = 13;
	caption = sprintf('Masked Original Image\nShowing Regions of Only the Specified Color');
	title(caption, 'FontSize', fontSize);
	% Show the original image next to it.
	subplot(3, 3, 7);
	imshow(rgbImage);
	title('The Original Image (Again)', 'FontSize', fontSize);

	% Measure the mean HSV and area of all the detected blobs.
	[meanHSV, areas, numberOfBlobs] = MeasureBlobs(coloredObjectsMask, hImage, sImage, vImage);
	if numberOfBlobs > 0
		fprintf(1, '\n----------------------------------------------\n');
		fprintf(1, 'Blob #, Area in Pixels, Mean H, Mean S, Mean V\n');
		fprintf(1, '----------------------------------------------\n');
		for blobNumber = 1 : numberOfBlobs
			fprintf(1, '#%5d, %14d, %6.2f, %6.2f, %6.2f\n', blobNumber, areas(blobNumber), ...
				meanHSV(blobNumber, 1), meanHSV(blobNumber, 2), meanHSV(blobNumber, 3));
		end
	else
		% Alert user that no colored blobs were found.
		message = sprintf('No blobs of the specified color were found in the image:\n%s', fullImageFileName);
		fprintf(1, '\n%s\n', message);
		uiwait(msgbox(message));
	end

	subplot(3, 3, 9);
	ShowCredits();
	message = sprintf('Done!\n\nThe demo has finished.\n\nLook the MATLAB command window for\nthe area and color measurements of the %d regions.', numberOfBlobs);
	uiwait(msgbox(message));
	
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	fprintf(1, '%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
return; % from SimpleColorDetection()
% ---------- End of main function ---------------------------------


%----------------------------------------------------------------------------
function [meanHSV, areas, numberOfBlobs] = MeasureBlobs(maskImage, hImage, sImage, vImage)
try
	[labeledImage, numberOfBlobs] = bwlabel(maskImage, 8);     % Label each blob so we can make measurements of it
	if numberOfBlobs == 0
		% Didn't detect any blobs of the specified color in this image.
		meanHSV = [0 0 0];
		areas = 0;
		return;
	end
	% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
	blobMeasurementsHue = regionprops(labeledImage, hImage, 'area', 'MeanIntensity');   
	blobMeasurementsSat = regionprops(labeledImage, sImage, 'area', 'MeanIntensity');   
	blobMeasurementsValue = regionprops(labeledImage, vImage, 'area', 'MeanIntensity');   
	
	meanHSV = zeros(numberOfBlobs, 3);  % One row for each blob.  One column for each color.
	meanHSV(:,1) = [blobMeasurementsHue.MeanIntensity]';
	meanHSV(:,2) = [blobMeasurementsSat.MeanIntensity]';
	meanHSV(:,3) = [blobMeasurementsValue.MeanIntensity]';
	
	% Now assign the areas.
	areas = zeros(numberOfBlobs, 3);  % One row for each blob.  One column for each color.
	areas(:,1) = [blobMeasurementsHue.Area]';
	areas(:,2) = [blobMeasurementsSat.Area]';
	areas(:,3) = [blobMeasurementsValue.Area]';
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	fprintf(1, '%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
return; % from MeasureBlobs()
	
	
%----------------------------------------------------------------------------
% Function to show the low and high threshold bars on the histogram plots.
function PlaceThresholdBars(plotNumber, lowThresh, highThresh)
try
	% Show the thresholds as vertical red bars on the histograms.
	subplot(3, 4, plotNumber); 
	hold on;
	yLimits = ylim;
	line([lowThresh, lowThresh], yLimits, 'Color', 'r', 'LineWidth', 3);
	line([highThresh, highThresh], yLimits, 'Color', 'r', 'LineWidth', 3);
	% Place a text label on the bar chart showing the threshold.
	fontSizeThresh = 14;
	annotationTextL = sprintf('%d', lowThresh);
	annotationTextH = sprintf('%d', highThresh);
	% For text(), the x and y need to be of the data class "double" so let's cast both to double.
	text(double(lowThresh + 5), double(0.85 * yLimits(2)), annotationTextL, 'FontSize', fontSizeThresh, 'Color', [0 .5 0], 'FontWeight', 'Bold');
	text(double(highThresh + 5), double(0.85 * yLimits(2)), annotationTextH, 'FontSize', fontSizeThresh, 'Color', [0 .5 0], 'FontWeight', 'Bold');
	
	% Show the range as arrows.
	% Can't get it to work, with either gca or gcf.
% 	annotation(gca, 'arrow', [lowThresh/maxXValue(2) highThresh/maxXValue(2)],[0.7 0.7]);

catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	fprintf(1, '%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
return; % from PlaceThresholdBars()


%---------------------------------------------------------------------------------------------------------------------------------
% Ask user what color they want for the onions and peppers images and set up pre-defined threshold values.
function [hueThresholdLow, hueThresholdHigh, saturationThresholdLow, saturationThresholdHigh, valueThresholdLow, valueThresholdHigh] = SetThresholds()
try
% 	button = menu('What color do you want to find?', 'yellow', 'green', 'red', 'white');
	% Menu with purple commented out because it's all around and the regionfill just ends up selecting the whole image.
	button = menu('What color do you want to find?', 'yellow', 'green', 'red', 'white', 'purple');
	% Use values that I know work for the onions and peppers demo images.
	switch button
		case 1
			% Yellow
			hueThresholdLow = 0.10;
			hueThresholdHigh = 0.14;
			saturationThresholdLow = 0.4;
			saturationThresholdHigh = 1;
			valueThresholdLow = 0.8;
			valueThresholdHigh = 1.0;
		case 2
			% Green
			hueThresholdLow = 0.15;
			hueThresholdHigh = 0.60;
			saturationThresholdLow = 0.36;
			saturationThresholdHigh = 1;
			valueThresholdLow = 0;
			valueThresholdHigh = 0.8;
		case 3
			% Red.
			% IMPORTANT NOTE FOR RED.  Red spans hues both less than 0.1 and more than 0.8.
			% We're only getting one range here so we will miss some of the red pixels - those with hue less than around 0.1.
			% To properly get all reds, you'd have to get a hue mask that is the result of TWO threshold operations.
			hueThresholdLow = 0.80;
			hueThresholdHigh = 1;
			saturationThresholdLow = 0.58;
			saturationThresholdHigh = 1;
			valueThresholdLow = 0.55;
			valueThresholdHigh = 1.0;
		case 4
			% White
			hueThresholdLow = 0.0;
			hueThresholdHigh = 1;
			saturationThresholdLow = 0;
			saturationThresholdHigh = 0.36;
			valueThresholdLow = 0.7;
			valueThresholdHigh = 1.0;
		otherwise
			% Purple
			hueThresholdLow = 0.76;
			hueThresholdHigh = 0.94;
			saturationThresholdLow = 0.33;
			saturationThresholdHigh = 0.67;
			valueThresholdLow = 0.1;
			valueThresholdHigh = 0.7;
	end
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	fprintf(1, '%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
return; % From SetThresholds()

%---------------------------------------------------------------------------------------------------------------------------------
% Display the MATLAB logo.
function ShowCredits()
try
% 	xpklein;
% 	surf(peaks(30));
	logoFig = subplot(3,3,9);
	caption = sprintf('A MATLAB Demo');
	text(0.5,1.15, caption, 'Color','r', 'FontSize', 18, 'FontWeight','b', 'HorizontalAlignment', 'Center') ;
	positionOfLowerRightPlot = get(logoFig, 'position');
	L = 40*membrane(1,25);
	logoax = axes('CameraPosition', [-193.4013 -265.1546  220.4819],...
		'CameraTarget',[26 26 10], ...
		'CameraUpVector',[0 0 1], ...
		'CameraViewAngle',9.5, ...
		'DataAspectRatio', [1 1 .9],...
		'Position', positionOfLowerRightPlot, ...
		'Visible','off', ...
		'XLim',[1 51], ...
		'YLim',[1 51], ...
		'ZLim',[-13 40], ...
		'parent',gcf);
	s = surface(L, ...
		'EdgeColor','none', ...
		'FaceColor',[0.9 0.2 0.2], ...
		'FaceLighting','phong', ...
		'AmbientStrength',0.3, ...
		'DiffuseStrength',0.6, ... 
		'Clipping','off',...
		'BackFaceLighting','lit', ...
		'SpecularStrength',1.0, ...
		'SpecularColorReflectance',1, ...
		'SpecularExponent',7, ...
		'Tag','TheMathWorksLogo', ...
		'parent',logoax);
	l1 = light('Position',[40 100 20], ...
		'Style','local', ...
		'Color',[0 0.8 0.8], ...
		'parent',logoax);
	l2 = light('Position',[.5 -1 .4], ...
		'Color',[0.8 0.8 0], ...
		'parent',logoax);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	fprintf(1, '%s\n', errorMessage);
	uiwait(warndlg(errorMessage));
end
return; % from ShowCredits()