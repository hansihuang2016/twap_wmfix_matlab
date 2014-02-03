clear all
close all

load reg1.mat
load reg1_d1.mat

%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';

%Get Downloaded Histdata.com FX Data
prices_all = ...
    csvread(strcat(workingdir,'DAT_ASCII_EURUSD_M1_201401_cleaned.csv'));

fx_names = {'Open'; 'High'; 'Low'; 'Close'};

%Create a time-series object using said data
prices_all_ts = fints(datenum(prices_all(:,1:6)),prices_all(:,7:10),fx_names);

%Extract start and end dates for equities extraction
datebounds = ftsbound(prices_all_ts,2);

%% MOVING ON TO THE EQUITIES STUFF NOW
%% Download equities data to get trading dates for equities. 
%   We will use these to get the dates for our FX data

tickers_americas = {'^MERV', '^MXX', '^BVSP', '^GSPTSE', '^GSPC', 'DOW',...
    '^RUT'};
names_americas = {'BuenosAires', 'Mexico', 'SaoPaulo', 'Toronto',...
    'SandP500', 'DowJones', 'Russell2000'};

tickers_europe = {'^ATX', '^BFX', '^FCHI', '^GDAXI', '^SSMI', '^FTSE'};
names_europe = {'Vienna', 'Brussels', 'Paris', 'Frankfurt', 'Swiss',...
    'London'};

tickers_asiapac = {'^AORD', '^HSI', '^BSESN', '^JKSE', '^KLSE',...
    '^N225', '^NZ50', '^STI', '^KS11', '^TWII'};
names_asiapac = {'Australia', 'HangSeng', 'Bombay', 'Jakarta', 'Malaysia',...
    'Nikkei225', 'NewZealand', 'Singapore', 'Korea', 'Taiwan'};

tickers_equities = [tickers_americas tickers_europe tickers_asiapac];
names_equities = [names_americas names_europe names_asiapac];

%connecting to yahoo! once to get basic info on one of the tickers
%This will populate the first column in the time series and help with
%creating a structure for the other equities time series
connect = yahoo
    equities_tmp = fetch(connect, tickers_equities(1),'Adj Close',...
            datebounds(1,:),datebounds(2,:),'d'); 

    prices_equities_adjclose_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
        names_equities(1));
        
    equities_tmp = fetch(connect, tickers_equities(1),'Open',...
            datebounds(1,:),datebounds(2,:),'d');
    
    prices_equities_open_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
        names_equities(1));
close(connect)

%connecting to yahoo! to get equities info
connect = yahoo
    %creating a loop (*sigh*) -- we have to do this because matlab does not
    %allow us to download more than one ticker at a time from yahoo...:P
    for ii = 1:size(tickers_equities,2)
        %creating a temporary vector for the yahoo download for ADJ CLOSE
        equities_tmp = fetch(connect, tickers_equities(ii),'Adj Close',...
            datebounds(1,:),datebounds(2,:),'d');        
        %converting to FTS object
        equities_tmp_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
            names_equities(ii));
        prices_equities_adjclose_ts = merge(prices_equities_adjclose_ts, equities_tmp_ts,...
            'DateSetMethod', 'intersection');
        
        %creating a temporary vector for the yahoo download for OPEN
        equities_tmp = fetch(connect, tickers_equities(ii),'Open',...
            datebounds(1,:),datebounds(2,:),'d');        
        %converting to FTS object
        equities_tmp_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
            names_equities(ii));
        prices_equities_open_ts = merge(prices_equities_adjclose_ts, equities_tmp_ts,...
            'DateSetMethod', 'intersection');                
    end
close(connect)

prices_equities_adjclose = fts2mat(prices_equities_adjclose_ts);

clear ii

%% Specify start- and end-hour here

%Specify start and end hours in 24-hour format 
%(e.g., 5 for 5am, 17 for 5pm)
starthour = 10
startminute = 00

endhour = 11
endminute = 30

starthr = starthour+(startminute/60);
endhr = endhour+(endminute/60);

%Extract start and end dates for equities extraction
datebounds = ftsbound(prices_all_ts,2);


%% Start data extraction for FX

%Extracting the end-hour and start-hour price for each day in the year

%Doing it for the end-hour first
%If the end-hour is 11am, this is our proxy for the 11am WM Fix
dv = cellstr(datestr(endhr/24:1/60/24:endhr/24));
prices_endhour_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
    datebounds(2,:), [], 1, 'd',dv);

%Doing it for the start-hour
dv = cellstr(datestr(starthr/24:1/60/24:starthr/24));
prices_starthour_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
    datebounds(2,:), [], 1, 'd',dv);

%Extract start and end dates
datebounds = ftsbound(prices_all_ts,2);

% prices_endhour_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
%     datebounds(2,:), [], 1, 'd',{'11:00'});

%% Data extraction for Equities

%get the series of the equities datenums
datenums_equities_prices = getfield(prices_equities_adjclose_ts,'dates');
datevec_equities_prices = datevec(datenums_equities_prices);

%%Now extracting from the end-hour prices only the days for which the datenums 
%%series is the smallest
prices_endhour_ts = prices_endhour_ts(datestr(datenums_equities_prices));

%Adding a timestamp to the equities time series because Matlab does not
%behave well without that timestamp...

%convert to non-FTS object
prices_equities_adjclose = fts2mat(prices_equities_adjclose_ts);
prices_equities_open = fts2mat(prices_equities_open_ts);

%create the timestamp
time = [endhour endminute 00];

equitiestime_prices = repmat(time, size(prices_equities_adjclose,1), 1);

%add the timestamp
prices_equities_adjclose = [datevec_equities_prices(:,1:3) equitiestime_prices ... 
    prices_equities_adjclose];
prices_equities_open = [datevec_equities_prices(:,1:3) equitiestime_prices ... 
    prices_equities_open];

%reconvert to time series object
prices_equities_adjclose_ts = fints(datenum(prices_equities_adjclose(:,1:6)),...
    prices_equities_adjclose(:,7:size(prices_equities_adjclose,2)), names_equities);
prices_equities_open_ts = fints(datenum(prices_equities_open(:,1:6)),...
    prices_equities_open(:,7:size(prices_equities_open,2)), names_equities);

%clearing vars for mem reasons
clear datevec_equities
   
%% Creating a separate time series object with the variables (including 
% interaction terms) that have been used in the stepwise regression

% [Note that I know I will probably find a one-step Matlab function that
% does this as soon as I have successfully created this time series
% object...]

%creating a times series of *just* the predictor variables that have come 
%out of the stepwise regression (NOT including interaction terms)
predictors_ts = extfield(prices_equities_adjclose_ts, ...
    reg1.PredictorNames);

%creating a separate variable for the NAMES of all the predictor variables
%including the interaction terms
coeff_names = reg1.CoefficientNames;

%Now trying to add the interaction terms to the original predictors time
%series by using the fact that Matlab spits out interaction terms using a
%colon (':')...this is terribly inefficient but it gets the job done
for jj = 1:size(coeff_names,2)
    %get an index value for the position of the colon (':') if it exists.
    %if it does not exist then k will be null
    k = strfind(coeff_names{jj},':');
    %check that k is not null (this means we have an interaction term)
    if isequal(k,'') == 0
        %create a separate variable since coeff_names is a cell array
        strtmp = coeff_names{jj};
        %replace the colon with an underscore ('_') because time series
        %object names cannot contain any other symbol
        strtmp(k) = '_';
        %get the first interaction term
        firstterm = strtmp(1:k-1);
        %get the second interaction term
        secondterm = strtmp(k+1:size(strtmp,2));
        %extract *first* term from the main time series object containing 
        %ALL lagged equities data
        getfirst = extfield(prices_equities_adjclose_ts, ...
            firstterm);
        %extract *second* term from the main time series object containing 
        %ALL lagged equities data
        getsecond = extfield(prices_equities_adjclose_ts, ...
            secondterm);
        %create the interaction variable
        interac = getfirst.*getsecond;
        %by default it gets named from the first variable so need to change
        %the name...
        interac = chfield(interac, strtmp(1:k-1), strtmp);
        %merging it with the time series of predictor variables
        predictors_ts = merge(predictors_ts, interac);
    end
    clear idx strtmp firstterm secondterm getfirst getsecond interac
end
clear idx strtmp firstterm secondterm getfirst getsecond interac        

%converting it to non-ts obj and adding a series of '1's to 
%to account for the intercept
predictors = fts2mat(predictors_ts);
predictors = [ones(size(predictors,1),1) predictors];

%% Making the predictions here

%Extract coefficient values from regression
coeffs = reg1.Coefficients.Estimate;

%Get predictions
predicTIONS = coeffs'*predictors';

%Get actuals
prices_endhour = fts2mat(prices_endhour_ts);

%Plot predictions vs. actuals
scatter(prices_endhour(:,4),predicTIONS)


    