# Bioengineering Master's Thesis Project
## Crowdsourcing Seizure Detection Algorithms Using Kaggle and ieeg.org
### University of Pennsylvania


## Upload annotations to ieeg.org
1. Export annotations from Persyst to csv format.
2. Store the csv files in Upload Annotations Pipeline/csv files.
3. Connect to ieeg.org using the toolbox and open the dataset you wish to annotate.
4. Run the upload annots.m function in MATLAB command window.

**Notes**
- csv files must be in correct format for proper annotation on ieeg.org. Check example csv file.

## Preprocessing
  The datasets must be clipped into one-second segments for classification by the seizure detection algorithms. Be prepared for the generation of thousands of .mat files. They will amount to several gigabytes in size.

**Steps**

1. Run the clipIEEGsegs.m file in matlab. Be sure to have your ieeg.org pw file in the current working directory. Alter the ieeg.org access command so that it reflects *your* username and pw file.
2. Store the EEG .mat files in a folder labeled 'Patient_#' where the # indicates the patient ID (1-8).

## Seizure Detection
**Run algorithms on clips**

  Each Kaggle algorithm is accompanied by a README file that explains where the clipped files need to be placed. Some of the code detailing clip location (directory name) may need to be changed. In addition, make certain you have all the necessary dependencies installed. The standard detector must be run first, however, as it also generates the key in the process of loading in the files.  

## Thesis Postprocessing
  
  The Thesis\_Processing.m file is very specific and the later graphs will not be accurate on new datasets. However, the ROC curves and AUC values should be accurate. Be sure that all submissions are in the current working directory and run Thesis\_Processing.m.
