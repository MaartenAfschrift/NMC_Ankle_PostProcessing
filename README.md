# NMC_Ankle_PostProcessing



This repository contains all the matlab scripts to analyze the data collected in the study (add ref to publication). I ran the code on Matlab2021b.

The code contains two main scripts

- *CreateDataMatrix.m:* downloads the data (.csv files) from dropbox and creates some large matrices that store all the outputs (muscle activity, COM movement, ...). These matrices are saved .mat files in a folder ./Data/ResultsFiles
- PlotDataAndRunStatistics: Loads all the outcomes, runs a statistical test and plots the presented in the paper.

You will notice that the results generated with this script are slightly different compared to the values presented in the paper. These small differences are caused by the interpolation of some signals (EMG, Forceplate, Exoskeleton sensor data) to reduce the size of the csv files. You can find more details in the readme file in the datafolder (Note that you can download the data by running the matlab function GetDataAfschrift2022() or using this URL https://www.dropbox.com/s/pwe0hkrwexhpqej/Data_NMC_Balance_Ankle.zip?dl=1).



