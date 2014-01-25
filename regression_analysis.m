
clear all
close all

%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';

load prices_all_ts.mat
% load prices_equities_open_ts.mat
load prices_equities_adjclose_ts.mat

year = datevec(prices_all_ts.dates(1));

year(1)

%% Do you want to take logs of equity prices?

takelogs = 0 %1 = yes; 0 = no

%% Do you want to run a lagged regression?

lagdecision = 1 %1 = yes; 0 = no

%% Do you want to calculate the TWAP between the start- and the end-hour?

twapdecision = 0 %1 = yes; 0 = no

%% Specify start- and end-hour here

%Specify start and end hours in 24-hour format 
%(e.g., 5 for 5am, 17 for 5pm)
starthour = 10
startminute = 00

endhour = 12
endminute = 00

starthour = starthour+(startminute/60);
endhour = endhour+(endminute/60);

%Extract start and end dates for equities extraction
datebounds = ftsbound(prices_all_ts,2);


%% Start data extraction for FX

%Extracting the end-hour and start-hour price for each day in the year

%Doing it for the end-hour first
%If the end-hour is 11am, this is our proxy for the 11am WM Fix
dv = cellstr(datestr(endhour/24:1/60/24:endhour/24));
prices_endhour_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
    datebounds(2,:), [], 1, 'd',dv);

%Doing it for the start-hour
dv = cellstr(datestr(starthour/24:1/60/24:starthour/24));
prices_starthour_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
    datebounds(2,:), [], 1, 'd',dv);

%Extract start and end dates
datebounds = ftsbound(prices_all_ts,2);


% prices_endhour_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
%     datebounds(2,:), [], 1, 'd',{'11:00'});

%% Data extraction for Equities

returns_equities_adjclose_ts = tick2ret(prices_equities_adjclose_ts);
returns_equities_open_ts = tick2ret(prices_equities_open_ts);

%get the equity index names
names_equities = fieldnames(prices_equities_adjclose_ts,1);

%get the series of the equities datenums
datenums_equities_rtns = getfield(returns_equities_adjclose_ts,'dates');
datevec_equities_rtns = datevec(datenums_equities_rtns);

datenums_equities_prices = getfield(prices_equities_adjclose_ts,'dates');
datevec_equities_prices = datevec(datenums_equities_prices);

%%Now extracting from the end-hour prices only the days for which the datenums 
%%series is the smallest
prices_endhour_ts = prices_endhour_ts(datestr(datenums_equities_rtns));

%Adding a timestamp to the equities time series because Matlab does not
%behave well without that timestamp...

%convert to non-FTS object
returns_equities_adjclose = fts2mat(returns_equities_adjclose_ts);
returns_equities_open = fts2mat(returns_equities_open_ts);

prices_equities_adjclose = fts2mat(prices_equities_adjclose_ts);
prices_equities_open = fts2mat(prices_equities_open_ts);

%%%%%%%%%%%%%%%%%%
%%Take logs of prices to make them comparable to FX data...
%%%%%%%%%%%%%%%%%%
if takelogs == 1
    prices_equities_adjclose = log(prices_equities_adjclose);
    prices_equities_open = log(prices_equities_open);
end

%create the timestamp
if endminute == 0
    time = [(endhour-1) 59 00];
else
    time = [(endhour-1) (endminute-1) 00];
end
equitiestime_rtns = repmat(time, size(returns_equities_adjclose,1), 1);
equitiestime_prices = repmat(time, size(prices_equities_adjclose,1), 1);

%add the timestamp
returns_equities_adjclose = [datevec_equities_rtns(:,1:3) equitiestime_rtns ...
    returns_equities_adjclose];
returns_equities_open = [datevec_equities_rtns(:,1:3) equitiestime_rtns ...
    returns_equities_open];

prices_equities_adjclose = [datevec_equities_prices(:,1:3) equitiestime_prices ... 
    prices_equities_adjclose];
prices_equities_open = [datevec_equities_prices(:,1:3) equitiestime_prices ... 
    prices_equities_open];

%reconvert to time series object
prices_equities_adjclose_ts = fints(datenum(prices_equities_adjclose(:,1:6)),...
    prices_equities_adjclose(:,7:size(prices_equities_adjclose,2)), names_equities);
prices_equities_open_ts = fints(datenum(prices_equities_open(:,1:6)),...
    prices_equities_open(:,7:size(prices_equities_open,2)), names_equities);

returns_equities_adjclose_ts = fints(datenum(returns_equities_adjclose(:,1:6)),...
    returns_equities_adjclose(:,7:size(returns_equities_adjclose,2)), names_equities);
returns_equities_open_ts = fints(datenum(returns_equities_open(:,1:6)),...
    returns_equities_open(:,7:size(returns_equities_open,2)), names_equities);

%clearing vars for mem reasons
clear datevec_equities prices_equities_adjclose prices_equities_close...
    prices_equities_open prices_equities_diff returns_equities_adjclose...
    returns_equities_close returns_equities_open

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATING THE TWAP HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if twapdecision == 1
    
    %Create vector of times and extract all prices from
    %starthour to starthour + 59mins daily
    dv = cellstr(datestr(starthour/24:1/60/24:endhour/24-1/60/24));
    prices_hourly_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
        datebounds(2,:), [], 1, 'd',dv);

    %% Calculating the TWAP between start-hour and end-hour

    %%FIRST: Convert to non-FTS object to be able to use the 
    %accumarray() function
    prices_hourly = fts2mat(prices_hourly_ts,1);
    prices_endhour = fts2mat(prices_endhour_ts,1);

    %%SECOND: Get a matrix of datevectors in the form [y, m, d, h, m, s] to use
    %%to use the unique() function and get a unique date for each time period.
    %%In other words, there will be a unique value for each set of values
    %%from 10:00 to 10:59am, 9:00 to 9:59am, and so on...
    dates_hourly = datevec(prices_hourly(:,1));

    %Now get the unique dates using the unique() function so that we can 
    %eventually get the TWAP for each of those dates
    [uniquedates_hourly, ~, subs_hourly]    = ...
        unique(dates_hourly(:,1:3),'rows');

    %%THIRD: Calculate the TWAP for each date using the accumarray() function
    %Doing it by column (open, high, low close) 
    prices_TWAP_open = [uniquedates_hourly accumarray(subs_hourly, ...
        prices_hourly(:,2), [], @mean)];
    prices_TWAP_high = [uniquedates_hourly accumarray(subs_hourly, ...
        prices_hourly(:,3), [], @mean)];
    prices_TWAP_low = [uniquedates_hourly accumarray(subs_hourly, ...
        prices_hourly(:,4), [], @mean)];
    prices_TWAP_close = [uniquedates_hourly accumarray(subs_hourly, ....
        prices_hourly(:,5), [], @mean)];

    %%FOURTH: We have to add a time column as we have an average per date
    %We are doing this because once we convert it back to a FTS object we will
    %need an identifier for each time
    TWAPtime = repmat(time, size(prices_TWAP_open,1), 1);
    prices_TWAP = [prices_TWAP_open(:,1:3) TWAPtime prices_TWAP_open(:,4)...
        prices_TWAP_high(:,4) prices_TWAP_low(:,4) prices_TWAP_close(:,4)];

    %Adding date+time columns for the end-hour price matrix for completeness
    prices_endhour = [datevec(prices_endhour(:,1)) prices_endhour(:,2:5)];

    %clearing up variables for mem reasons
    clear prices_all time prices_TWAP_open prices_TWAP_high...
        prices_TWAP_low prices_TWAP_close subs_hourly dv dates_hourly TWAPtime

    %%FIFTH: Merge all data series based on dates where necessary

    %Create a vector of datenums for the start- and end-hour prices
    datenums_endhour = getfield(prices_endhour_ts,'dates');
    datenums_starthour = getfield(prices_starthour_ts,'dates');
    
    %Create a vector of datenums for the TWAP prices
    datenums_TWAP = datenum(prices_TWAP(:,1:3));
    
    %Merge all the datenums
    datenums_fx = intersect(datenums_endhour, datenums_starthour, ...
        datenums_TWAP, 'rows');    

    %Converting the hourly TWAPs to FTS object
    prices_TWAP_ts = fints(datenum(prices_TWAP(:,1:6)), ...
        prices_TWAP(:,7:10),fx_names);
    
else
    %% Create a vector of datenums for the end-hour prices
    
    %This is all we need to do if we are not calculating the TWAP
    datenums_endhour = getfield(prices_endhour_ts,'dates');
    datenums_starthour = getfield(prices_starthour_ts,'dates');
    datenums_fx = intersect(datenums_endhour, datenums_starthour, 'rows');
end

%% Merge all data series based on dates where necessary

%Keep the intersection of all the fx and the equities dates
common_datenums_all = intersect(datenums_fx, datenums_equities_prices,...
    'rows');

%Now that we have the common dates, let's make everything match up to those
%dates
prices_endhour_ts = prices_endhour_ts(datestr(common_datenums_all));
prices_starthour_ts = prices_starthour_ts(datestr(common_datenums_all));

returns_equities_adjclose_ts = returns_equities_adjclose_ts(datestr(common_datenums_all));
returns_equities_open_ts = returns_equities_open_ts(datestr(common_datenums_all));

prices_equities_adjclose_ts = prices_equities_adjclose_ts(datestr(common_datenums_all));
prices_equities_open_ts = prices_equities_open_ts(datestr(common_datenums_all));

%clearing up variables for mem reasons
clear missingdates datenums_hourly datenums_endhour rowidx_nan...
    check_missing nan_entry nan_entry_ts common_datenums datenums_starthour...
    datenums_equities_prices datenums_equities_rtns

%% Converting FX FTS objects back to non-FTS objects to use in regression

prices_endhour = fts2mat(prices_endhour_ts);
prices_starthour = fts2mat(prices_starthour_ts);


%% Doing the lagging and lining up of variables

%Lagging equities adj close by 1 period
prices_equities_adjclose_ts_lagged1 = ...
    lagts(prices_equities_adjclose_ts,1,NaN);
prices_equities_adjclose_lagged1 = ...
    fts2mat(prices_equities_adjclose_ts_lagged1);

%lining up the one-period lagged equities and removing the last period 
%of end-hour prices to line up previous day equity prices with current 
%day's FX values

%removing the first 'NaN' value in the lagged equities data
prices_equities_adjclose_lagged1 = ...
    prices_equities_adjclose_lagged1(2:size(prices_equities_adjclose_lagged1,1),:);

%removing the last day of the year for the FX values
prices_lagged_endhour = prices_endhour(1:size(prices_endhour,1)-1,:);
prices_lagged_starthour = prices_starthour(1:size(prices_starthour,1)-1,:);

prices_equities_adjclose = fts2mat(prices_equities_adjclose_ts);
prices_equities_open = fts2mat(prices_equities_open_ts);

%% Stepwise regression on Close of End-Hour FX and Adj Close of equities

if lagdecision == 1
    %running the regression on lagged equity/unlagged FX
    reg1 = stepwiselm(prices_equities_adjclose_lagged1, ...
        prices_lagged_endhour(:,4), ...
        'VarNames',[names_equities' 'EndHourPrices'])
%     reg4 = stepwiselm(prices_equities_adjclose_lagged1, ...
%         prices_lagged_starthour(:,4), ...
%         'VarNames',[names_equities 'StartHourPrices'])

else
    %running regression on unlagged equity/unlagged FX
    reg1 = stepwiselm(prices_equities_adjclose, prices_endhour(:,4))
end

%% Creating a separate time series object with the variables (including 
% interaction terms) that have been used in the stepwise regression

% [Note that I know I will probably find a one-step Matlab function that
% does this as soon as I have successfully created this time series
% object...]

%creating a times series of *just* the predictor variables that have come 
%out of the stepwise regression (NOT including interaction terms)
predictors_ts = extfield(prices_equities_adjclose_ts_lagged1, ...
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
        getfirst = extfield(prices_equities_adjclose_ts_lagged1, ...
            firstterm);
        %extract *second* term from the main time series object containing 
        %ALL lagged equities data
        getsecond = extfield(prices_equities_adjclose_ts_lagged1, ...
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

%converting it to non-ts obj to use in regression
predictors = fts2mat(predictors_ts);
%getting all the names
names_predictors = fieldnames(predictors_ts,1);

%removing the first 'NaN' value in the predictors data to line it up with
%the FX data
predictors = predictors(2:size(predictors,1),:);

%% Linear multivariate regression on Close of 11am FX and Adj Close 
% of equities

%this is a check to see that our linear multivariate regression matches up
%with the stepwise regression
reg2 = fitlm(predictors, prices_lagged_endhour(:,4), ...
    'VarNames',[names_predictors' 'EndHourPrices'])

save predictors_ts.mat predictors_ts
save reg1.mat reg1
save reg2.mat reg2

