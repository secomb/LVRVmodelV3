# LVRVmodelV3
**Fast computational model for analysis of interventricular interactions**  

For details, see: Analysis of interventricular interactions using a fast computational model, by Vivian Nguyen, Timothy W. Secomb and Michael J. Moulton (in review)  
This version of the program includes parameter estimation, based on target values of output variables.  
Multiple sets of input files, labeled 000, 001, 002 etc. can be processed. In these files, lines starting // are ignored.  
A number of switches are specified in parameters000.txt. If switches are off, a single forward simulation is performed. The switches are as follows:  
Optimize = on: Run parameter optimization using the Levenberg-Marquardt algorithm  
Skip = on: Skip this file, useful when processing multiple files  
odeSystem = LVRVmodel2024 (or LVmodel2024 for comparison with axisymmetric LV model)  
septalDisp = on (or off for comparison with model with no septal displacement)  
RVhypertrophy = on (this is needed to give a RV free wall with realistic thickness)  
LVAD = off (not fully implemented)  
pericardium = off (option to include effects of constrictive pericardium)  
 
Initial parameter values are specified in parameters000.txt. If parameter "param" is to be optimized, then the line  
ipar_param = n  
should be included in the file, where the nth multiplicative factor xx[n] is optimized. This allows for more than one input parameter to be varied by a single factor. The factors are initially set to 1 and their upper and lower bounds are specified in the last section of the file.  
Target values of output properties are specified in targets000.txt. Target "tar" is included in the objective function if the line  
wt_tar = w  
is included in the file, where w is a relative weight to be applied to target variable tar, for example 1 or 10.  
Several output files are generated, as follows.    
graphs000.ps: A set of graphs of relevant quantities in PostScript  
outputFile_All000.txt: All variables as functions of time, can be copied into a spreadsheet  
residuals000.txt: Fitted and target values of output variables, showing deviations  
chisqVals000.txt: Variation of objective function during optimization  

We have tested this package using Microsoft Visual C++ 2022 under Windows 11. For error reporting and suggestions please contact Dr. Timothy W. Secomb, secomb@arizona.edu. We welcome your comments and suggestions.  

This program is freely available for non-commercial use, provided appropriate acknowledgement is given. Commercial users please contact us before using this program. No assurance is given that it is free of errors and any use is at the user’s risk.  
