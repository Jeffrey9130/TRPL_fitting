# Introduction
This program is used to extract material parameters from the TRPL-type
data. Some features are summarized below:
- It can analyze the surface recombination velocity, bulk lifetime and diffusion coefficient;
- It supports 1D, 2D and 3D diffusion model;
- It supports single curve fitting and global fitting among several curves;
- All fitting results will be saved automatically in CSV files

# Dependence
You need to install matlab toolbox `v2struct`.

# How to use
## 1. Raw data preparation
### 1.1 TRPL data
The TRPL raw data should be saved as matlab .mat file at *rawdata/SAMPLE_NAME.mat* following the structure below.


`SAMPLE_NAME.mat`  
+-`lambda_alpha.mat`  
+-`CURVE_1.mat`  
+-`CURVE_2.mat (if you want to run a global fitting)`  
+- ...

The lambda_alpha.mat saves the excitation wavelength (nm), absorption coefficient at the excitation wavelength (cm^-1), emission wavelength(nm), absorption coefficient at the emission wavelength (cm^-1) as the structure below:

`lambda_alpha.mat`

|         |column1                     |column2                            | ...|
| --------|----------------------------|-----------------------------------|----|
|row1     |excitation_lambda_of_CURVE1 | excitation_lambda_of_CURVE1       | ...|
|row2     |excitation_alpha_of_CURVE1  | excitation_alpha_of_CURVE1        | ...|
|row3     |emission_lambda_of_CURVE1   | emission_lambda_of_CURVE1         | ...|
|row2     |emission_alpha_of_CURVE1    | emission_alpha_of_CURVE1          | ...|

In each of the CURVE.mat file, experimental data are saved as below

`CURVE.mat`

|column1  |  column2 |
|---------| ---------|
|time_step1| data_1 |
|time_step2| data_2 |
| ...| ...|

## 2. Main program
### 2.1 data_processing(cutoff_time, plot_tag, sampling_ratio)
- `cutoff_time`: The end time for the fitting curves, use the same time unit of that in the raw data;
- `plot_tage`: Set 1 to plot the rawdata;
- `sampling_ratio`: Uniformaly spaced sample the raw data by a ratio to reduce the data amount and save fitting time. Set 1 to use all the data, set 0.5 to drop half of the data (fetch one data at every other point)
### 2.2 Control parameters
    See the comments for more details
### 2.3 Initialze fitting parameters
All parameters below will be used in the fitting, the frist num_global paramters will be used as the global fitting paramters.
- `tau`: bulk liftime
- `D`: diffusion coeffient
- `S1`: surface recombination velocity of the excitated surface
- `S2`: surface recombination velocity of the backside surface
- `Amp`: amplitude of the curve
- `offset`: y-offset of the curve
- `small_alpha`: smooth parameter of the initialized carrier distribution
                     See pdex1ic.m for more details.
## 3. Output
Fitting results are save in the *output/* folder.
- *p.mat*                     original fitting results
- *fitting_result.csv*        fitting results
- *curve_exp.csv*             experimental data (time v.s. pl)
- *curve_sim.csv*             fitted data       (time v.s. pl)
