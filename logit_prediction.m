clear all
close all

load B.mat 
load stats.mat

%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';

%Get Downloaded Histdata.com FX Data
prices_all = ...
    csvread(strcat(workingdir,'DAT_ASCII_EURUSD_M1_2014_cleaned.csv'));

fx_names = {'Open'; 'High'; 'Low'; 'Close'};

%Create a time-series object using said data
prices_all_ts = fints(datenum(prices_all(:,1:6)),prices_all(:,7:10),fx_names);

%Extract start and end dates for equities extraction
datebounds = ftsbound(prices_all_ts,2);

lagperiod = 1 %set number of days you want to lag equities here

%% MOVING ON TO THE EQUITIES STUFF NOW
%% Download equities data to get trading dates for equities. 
%   We will use these to get the dates for our FX data

tickers_americas = {'^MERV', '^MXX', '^BVSP', '^GSPTSE', '^GSPC',...
    '^RUT'};
names_americas = {'BuenosAires', 'Mexico', 'SaoPaulo', 'Toronto',...
    'SandP500', 'Russell2000'};

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

    prices_equities_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
        names_equities(1));
        
    equities_tmp = fetch(connect, tickers_equities(1),'Open',...
            datebounds(1,:),datebounds(2,:),'d');
    
%     prices_equities_open_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
%         names_equities(1));
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
        prices_equities_ts = merge(prices_equities_ts, equities_tmp_ts,...
            'DateSetMethod', 'intersection');
        
        %creating a temporary vector for the yahoo download for OPEN
%         equities_tmp = fetch(connect, tickers_equities(ii),'Open',...
%             datebounds(1,:),datebounds(2,:),'d');        
        %converting to FTS object
%         equities_tmp_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
%             names_equities(ii));
%         prices_equities_open_ts = merge(prices_equities_ts, equities_tmp_ts,...
%             'DateSetMethod', 'intersection');                
    end
close(connect)

clear ii equities_tmp_ts

%% GETTING ACTUAL DATA FOR TESTING ACCURACY OF LOGIT PREDICTIONS
%% STEP 1: Specify start- and end-hour here

%Specify start and end hours in 24-hour format 
%(e.g., 5 for 5am, 17 for 5pm)
starthour = 10
startminute = 00

endhour = 11
endminute = 00

fixhour = 11+(00/60);

starthr = starthour+(startminute/60);
endhr = endhour+(endminute/60);

time = [endhour endminute 00];

%Extract start and end dates for equities extraction
datebounds = ftsbound(prices_all_ts,2);

%% STEP 2: Data extraction for FX

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

%Doing it for the fix
dv = cellstr(datestr(fixhour/24:1/60/24:fixhour/24));
prices_11amfix_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
    datebounds(2,:), [], 1, 'd',dv);

%Extract start and end dates
datebounds = ftsbound(prices_all_ts,2);

%Get the OHLC cell array from FX data
fx_names = fieldnames(prices_all_ts,1);

% prices_endhour_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
%     datebounds(2,:), [], 1, 'd',{'11:00'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATING THE TWAP HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create vector of times and extract all prices from
%starthour to starthour + 59mins daily
dv = cellstr(datestr(starthr/24:1/60/24:endhr/24-1/60/24));
prices_hourly_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
    datebounds(2,:), [], 1, 'd',dv);

%% STEP 3: Calculating the TWAP between start-hour and end-hour

%%FIRST: Convert to non-FTS object to be able to use the 
%accumarray() function
prices_hourly = fts2mat(prices_hourly_ts,1);


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


%clearing up variables for mem reasons
clear prices_all prices_TWAP_open prices_TWAP_high...
    prices_TWAP_low prices_TWAP_close subs_hourly dv dates_hourly TWAPtime


%%FIFTH: Merge all data series based on dates where necessary

%Create a vector of datenums for the start- and end-hour prices
datenums_endhour = getfield(prices_endhour_ts,'dates');
datenums_starthour = getfield(prices_starthour_ts,'dates');
datenums_11amfix = getfield(prices_11amfix_ts,'dates');

%Create a vector of datenums for the TWAP prices
datenums_TWAP = datenum(prices_TWAP(:,1:3));

%Merge all the datenums
datenums_fx = intersect(datenums_11amfix, datenums_TWAP, 'rows');    

%Converting the hourly TWAPs to FTS object
prices_TWAP_ts = fints(datenum(prices_TWAP(:,1:6)), ...
    prices_TWAP(:,7:10),fx_names);



%% STEP 4: Removing entries for which common dates do not exist

%get the series of all of the datenums and matching up equities with the 
%FX dates
datenums_equities_prices = getfield(prices_equities_ts,'dates');

datenums_endhour = getfield(prices_endhour_ts,'dates');
datenums_starthour = getfield(prices_starthour_ts,'dates');
datenums_fx = intersect(datenums_endhour, datenums_starthour, 'rows');

common_datenums_all = intersect(datenums_fx, datenums_equities_prices,...
    'rows');

datevec_prices = datevec(common_datenums_all);

%%Now extracting from the end-hour prices only the days for which the datenums 
%%series is the smallest
prices_endhour_ts = prices_endhour_ts(datestr(common_datenums_all));
prices_starthour_ts = prices_starthour_ts(datestr(common_datenums_all));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATING TWAP RESPONSE VARIABLE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prices_11amfix_ts = prices_11amfix_ts(datestr(common_datenums_all));
prices_TWAP_ts = prices_TWAP_ts(datestr(common_datenums_all));
% prices_diff_ts = prices_TWAP_ts - prices_11amfix_ts;
prices_diff_ts = prices_11amfix_ts - prices_starthour_ts;

prices_diff = fts2mat(prices_diff_ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Chopping off the first date in the response variable to make sure that the
%lagged equities data lines up with the FX data
prices_lagged_diff = prices_diff(lagperiod+1:size(prices_diff,1),:);

%Creating a nominal variable here to categorize when the TWAP is greater
%than the fix or vice versa 
%(Recall that prices_diff_ts = prices_TWAP_ts - prices_11amfix_ts;)
%%%
for ii = 1:size(prices_lagged_diff,1)
    if prices_lagged_diff(ii,4) > 0
        %TWAP greater than 11am
        prices_lagged_diff_logical(ii) = 2;
    else
        %11am greater than TWAP
        prices_lagged_diff_logical(ii) = 1;
    end
end
clear ii


%% LOGIT PREDICTIONS HAPPENING HERE

% FIRST we need to log and lag equities data

%adding a lag
prices_equities_ts_lagged = ...
    lagts(prices_equities_ts,lagperiod,NaN);


%convert to non-FTS object
prices_equities_lagged = fts2mat(prices_equities_ts_lagged);


%take logs of equities
prices_equities_log_lagged = log(prices_equities_lagged);

%removing the first 'NaN' value in the lagged equities data
prices_equities_log_lagged = ...
    prices_equities_log_lagged(lagperiod+1:size(prices_equities_log_lagged,1),:);


% SECOND we need to get our probabilities
pihat = mnrval(B, prices_equities_log_lagged);
pihat = [pihat prices_lagged_diff(:,4)];

% THIRD we want to put dates on the probabilities

%Knock the first date entry off (since it's lagged)
datenums_lagged_equities_prices = ... 
    datenums_equities_prices(lagperiod+1:size(datenums_equities_prices,1),:);
% datevec_lagged_equities_prices = datevec(datenums_lagged_equities_prices);

names = {'LessThanZero', 'GreaterThanZero', 'Actual'};
pihat_ts = fints(datenums_lagged_equities_prices, pihat, names, 'DAILY');
clear names

%% CHECK ACCURACY (FINALLY!)

for ii = 1:size(pihat,1)
    if pihat(ii,1) > .5
        check(ii) = 1;
    else check(ii) = 2;
    end
end
clear ii

yes = 0;
for ii = 1:size(pihat,1)
    if check(ii) == prices_lagged_diff_logical(ii)
        yes = yes + 1;
    end
end

% Percent correct
yes/size(pihat,1)

names = {'Prediction', 'Actual'};
check_ts = fints(datenums_lagged_equities_prices, ...
    [check' prices_lagged_diff_logical'], names, 'DAILY');
clear names
