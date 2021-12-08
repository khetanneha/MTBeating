/* Neha Khetan, June - Sep 2021
This folder contains the scripts used for the analysis of the experimental tracks obtained for the Single filament beating/movement project
*/


%% ============================================================================


||       ------------------------------------------                   READ ME      --------------------------------------------------------------------------------------------------     ||

To run the script you need to provide the following:

1. filename
2. path
3. Specify  "selectType" -  to choose for the source and method 
					Source and Method by which filament contour coordinates were obtained
							1. Neuron J by Shivani and 
							2. Custom code by Dhruv Khatri 
					In the script choose the following option "selectType" accordingly
					% Choose 0 = fro data from SY: excel - sheet (.xlsx) format from SY
					%        1 = for data from DK: Filament tracking - .mat format

					selectType = 0;
4. If using File "run_ExpAnalysis.m" please add the folder "export_fig-master.tar.xz" to your path
    else,
    please use the file "run_ExpAnalysis_V2.m"
    

5. If interpolation of the contour points is required, please assign "1" to the following
   						% To Interpolate points: for improving data quality ?
						%   Methods: Linear, polynomial, spline
									InterPolatePolyfit  = 0  ; % for polynomial based
									InterPolateSpline   = 0  ;
									InterPolateLinear   = 0 ;


%% ============================================================================




FILES that deals with these various options and methods.

1. coordinates_V2.mat is the data shared by Dhruv


2. I have uploaded 3 pdf files in the sub-folder 'example' to illustrate the #3 and #5
           
           a. Comparison between data from RIGOR, ANTIBODY and BIOTIN -steptavidin. With and without polynomial interpolation.
               Filename: ComparingExptDataTypes_MethodAnalysis.pdf : 
               
           b. Comparison between data from Shivani yadav ( Neuron J) and Dhruv Khatri (Custom script) for the Biotin-streptavidin data.
               Filename: ComparingExptAnalysis_DK_SYdata.pdf
               
           c. Comparison between different method for Rigor data from ( Neuron J) - raw data and interpolation methods of linear, spline and polynomial
               Filename: ComparingSimExpt_v4B.pdf

%% ============================================================================
