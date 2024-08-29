# 4.-Extracting-Features-from-Radiance-Gray-Temperatures-and-Emissivity-Time-Series

## Radiance_Temp_Analyzer ##
The script processes and analyzes radiance and temperature data from text files with the following steps:

### 1. File Handling: ###

    Reads radiance and temperature data from text files.

### 2. Data Processing: ###

    - Fitting and Derivatives:

      * Fits the radiance data using PchipInterpolator.
      
      * Computes derivatives and cumulative integrals of the fitted data.
      
    - Peak and Saddle Point Analysis:

      * Identifies significant peaks and saddle points in the radiance data.
      
      * Estimates temperatures at these points and computes related metrics (e.g., peak values, saddle points, temperature errors).
      
    - Area Calculations:

      * Computes the area under the radiance curve in different sections (shock rise, decay, and growth).
      
      * Computes normalized areas based on cumulative integrals.
      
### 3. Data Storage: ###

    Saves the processed results into an Excel file with separate sheets for each data file.
    
### 4.Plotting (commented out): ###

    Functions for plotting the data and results are defined but not executed in the current script.
    
The script automates the analysis of multiple datasets, calculates various metrics, and organizes the results for further examination.

## Radiance_Analyzer ##
Processes only the radiance data 

## Emissivity_Temp_Analyzer ##
This Python script processes and analyzes radiance, temperature, and emissivity data from multiple text files, performing several tasks and generating output for further examination. Here's a summary of its functionality:

### 1. Data Import: ###
    
    - User inputs the folder path containing the text files.
    
    - Radiance, temperature, and emissivity data are read from text files.

### 2. Data Processing: ###

    - Each dataset is processed to clean and prepare it for analysis.
    
    - Functions are used to fit, differentiate, and integrate data.
    
### 3. Interactive Analysis: ###

    - User reviews and decides whether to process or skip each dataset.
    
    - The script prompts for time intervals and identifies extrema and peak values.
    
### 4.Visualization: ###

    Plots are generated to help visualize the data and its characteristics.
    
### 5.Results Export: ###

    Processed data and results are saved into Excel files for further use.
    
This script is designed for detailed analysis of radiance, temperature, and emissivity data, providing visualization and statistical insights into combustion or deflagration experiments.
