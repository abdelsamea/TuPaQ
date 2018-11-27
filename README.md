# TuPaQ
Tumour Parcellation and Quantification (TuPaQ) software for rapid and automated segmentation of tumour epithelium

what is TuPaQ?

Tumour Parcellation and Quantification (TuPaQ) is a standalone Windows-compatible application implemented in Matlab for tumour epithelium identification and quantification in Colorectal Cancer (CRC) tissue images scanned by Light microscopy (e.g., dark objects of interests on a bright background).

Software Installation

You need to simply install the APP ‘MATLAB App installer’ (e.g., .mlappinstall) and it will automatically appear APPs tab in ‘MY APPS’.

If you would like to further develop TuPaQ using Matlab, please run TIDemo.m as the main file of TuPaQ but in this case please make sure that you have MATLAB 2016a or higher installed and also make sure that the following toolboxes are installed as well: •	Mapping Toolbox •	Computer Vision system Toolbox •	Statistics and Machine learning Toolbox •	Fuzzy logic Toolbox • Image processing Toolbox

Software Components

The TuPaQ contains the following components:

Epithelium segmentation
Stroma segmentation
Tumour epithelium identification
Actual tumour size extraction
Nuclei cluster division in tumour epithelium
Nuclei cluster division in Stroma
Usage instructions and software settings:

Phase 1: tumour segmentation

1- Upload your input RGB (MxNx3) image, the input image should be in one of the following format: “jpg, tiff, tif, png, bmp, or gif”.

2- The recommended value for downscaling parameter is 0.1 for an image/section of 40X level of magnification for accurate Nuclei detection.

3-	The recommended “number of classes” parameter is 3 (e.g., epithelium, stroma, and background). However, it can be set to 2 in some cases for accurate epithelium extraction and integration of stroma and background in one class.

4-	“Smoothness parameter” used to control the smoothness of the contour of the segmentation model.

5-	Set the number of pixels for 1mm for actual tumour size calculation to produce the actual tumour epithelium size in mm.

Phase 2: Tumour Quantification (Nuceli counts)

1-	Tick to automatically calculate the minimum depth limit (it is recommended for noisy-free and high quality images) or set the minimum depth limit manually (0.01 is recommend but it can be increased to reduce the number of cell cluster division into individual cells).

2-	Set minimum size of the nucleus in pixels.

3-	Press to tumour quantification, for example to calculate the ratio between the number of tumour Nuclei and total number of Nuclei in input image.

Stain Normalization phase:

Press “Reference Image” to select a new reference image for stain normalization.

Training phase:

1-	Tick to start a new training session based on your own datasets including some tumour epithelim and stroma cases.

2-	Set 2D map width of Self Organizing Neural network Classifier (3 is recommended for robustness behaviour)

3-	Set 2D map height of Self Organizing Neural network Classifier.

4-	Set the learning rate of the Classifier (0.1 is recommended to reduce the possibility of trapping to local minima).

5-	Press “SOMnormal” to select one folder that contains your normal cases. Please make sure to store your training images original image followed by its associate binary mask (that indicates the normal regions), consecutively. The same should be done with “SOMtumour” button for a separate folder contains your tumour training samples.

Disclaimer

If you find TuPaQ useful and use it in your publications, please cite the following articles below.

M. M. Abdelsamea, A. Pitiot, R. Barbora, J. Besusparis, A. Laurinavičius, M. Ilyas, A cascade-learning approach for automated segmentation of tumor epithelium in colorectal cancer, Expert Systems with Applications, Volume 118, 2019, Pages 539-552.

M. M. Abdelsamea, R. Barbora, J. Besusparis, S. Cham, A. Pitiot, A. Laurinavičius, M. Ilyas, Tumour Parcellation and Quantification (TuPaQ): a tool for refining biomarker analysis through rapid and automated segmentation of tumour epithelium, submitted to Clinical Cancer Research, 2018.

Contact

Please do not hesitate to send your comment/feedback to Mohammed Abdelsamea (m.abdelsamea@aun.edu.eg) if you experienced any problems.
