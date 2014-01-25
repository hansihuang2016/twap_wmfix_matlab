clear all
close all

load reg1.mat reg2.mat 
load prices_all_ts.mat prices_equities_adjclose_ts.mat
load predictors_ts.mat

%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';

%% Specify start/end years and currency
    
currency = 'EURUSD' %ALLCAPS please

startyear = 2008
endyear = 2013

yearvector = startyear:endyear;

for j = 1:size(yearvector,2)
    year = num2str(yearvector(j));

    %Get Downloaded Histdata.com FX Data
    prices_multiyear = ...
        csvread(strcat(workingdir,...
        'DAT_ASCII_',currency,'_M1_',year,'_cleaned.csv'));
    
    %Create a time-series object using said data
    prices_multiyear_ts = fints(datenum(prices_multiyear(:,1:6)),...
        prices_multiyear(:,7:10),fx_names);
    
    