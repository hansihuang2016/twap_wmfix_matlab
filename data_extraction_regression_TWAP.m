clear all
close all

%This script will extract the CLEANED data downloaded from the Histdata.com
%database, check for any missing dates or extra dates and adjust as
%necessary.

%Note that:
%(1) The data MUST be in the format specified by one of the csv files. That
%is, the csv file must have the year, month, date, hour, minute and second
%in a separate cell. Also, the volume information must be deleted.
%
%(2) For missing data 'NaN's are added in. We can perform linear
%interpolation (or any other operation) on missing data if necessary.
%
%(3) The data from Yahoo! should already exist in the workspace. Make sure
%you run yahoo_dat BEFORE you run this code.



%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';

%% Pull and extract FX data

%Get Downloaded Histdata.com FX Data
prices_all = ...
    csvread(strcat(workingdir,'DAT_ASCII_EURUSD_M1_2013_cleaned.csv'));

fx_names = {'Open'; 'High'; 'Low'; 'Close'};

%Create a time-series object using said data
prices_all_ts = fints(datenum(prices_all(:,1:6)),prices_all(:,7:10),fx_names);

save prices_all_ts.mat prices_all_ts

%Extract start and end dates for equities extraction
datebounds = ftsbound(prices_all_ts,2);

%% THIS IS ALL THE EQUITIES STUFF NOW
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

returns_equities_adjclose_ts = tick2ret(prices_equities_adjclose_ts);
returns_equities_open_ts = tick2ret(prices_equities_open_ts);

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
clear connect equities_tmp equities_tmp_ts ii equitiestime...
    datevec_equities prices_equities_adjclose prices_equities_close...
    prices_equities_open prices_equities_diff returns_equities_adjclose...
    returns_equities_close returns_equities_open



%% MOVING BACK TO THE FX STUFF NOW: CALCULATING TWAP
%% First, extract relevant data

%Create vector of times and extract all prices from
%starthour to starthour + 59mins daily
dv = cellstr(datestr(starthour/24:1/60/24:endhour/24-1/60/24));
prices_hourly_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
    datebounds(2,:), [], 1, 'd',dv);



%% Next, calculating TWAP for the hour

%%First: Convert to non-FTS object to be able to use the 
%accumarray() function
prices_hourly = fts2mat(prices_hourly_ts,1);
prices_endhour = fts2mat(prices_endhour_ts,1);



%%Second: Get a matrix of datevectors in the form [y, m, d, h, m, s] to use
%%to use the unique() function and get a unique date for each time period.
%%In other words, there will be a unique value for each set of values
%%from 10:00 to 10:59am, 9:00 to 9:59am, and so on...
dates_hourly = datevec(prices_hourly(:,1));

%Now get the unique dates using the unique() function so that we can 
%eventually get the TWAP for each of those dates
[uniquedates_hourly, ~, subs_hourly]    = ...
    unique(dates_hourly(:,1:3),'rows');



%%Third: Calculate the TWAP for each date using the accumarray() function
%Doing it by column (open, high, low close) 
prices_TWAP_open = [uniquedates_hourly accumarray(subs_hourly, ...
    prices_hourly(:,2), [], @mean)];
prices_TWAP_high = [uniquedates_hourly accumarray(subs_hourly, ...
    prices_hourly(:,3), [], @mean)];
prices_TWAP_low = [uniquedates_hourly accumarray(subs_hourly, ...
    prices_hourly(:,4), [], @mean)];
prices_TWAP_close = [uniquedates_hourly accumarray(subs_hourly, ....
    prices_hourly(:,5), [], @mean)];



%%Fourth: We have to add a time column as we have an average per date
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


%% Third, merge all data series based on dates where necessary

%Create a vector of datenums for the end-hour prices
datenums_endhour = getfield(prices_endhour_ts,'dates');

%Create a vector of datenums for the TWAP prices
datenums_TWAP = datenum(prices_TWAP(:,1:3));

%Converting the hourly TWAPs to FTS object
prices_TWAP_ts = fints(datenum(prices_TWAP(:,1:6)),...
    prices_TWAP(:,7:10),fx_names);


%Keep only the intersection of the end hour and the TWAP dates
common_datenums_fx = intersect(datenums_endhour,datenums_TWAP,'rows');
%Keep the intersection of the above, along with the equities dates
common_datenums_all = intersect(common_datenums_fx, datenums_equities_rtns,...
    'rows');

%Now that we have the common dates, let's make everything match up to those
%dates
prices_TWAP_ts = prices_TWAP_ts(datestr(common_datenums_all));
prices_endhour_ts = prices_endhour_ts(datestr(common_datenums_all));
prices_starthour_ts = prices_starthour_ts(datestr(common_datenums_all));

returns_equities_adjclose_ts = returns_equities_adjclose_ts(datestr(common_datenums_all));
returns_equities_open_ts = returns_equities_open_ts(datestr(common_datenums_all));

prices_equities_adjclose_ts = prices_equities_adjclose_ts(datestr(common_datenums_all));
prices_equities_open_ts = prices_equities_open_ts(datestr(common_datenums_all));

%clearing up variables for mem reasons
clear missingdates datenums_hourly datenums_endhour rowidx_nan...
    check_missing nan_entry nan_entry_ts common_datenums

%% Finally, taking the difference between the TWAP and the end-hour price

prices_diff_fx_ts = prices_TWAP_ts - prices_endhour_ts;
prices_diff_fx = fts2mat(prices_diff_fx_ts);

%% Converting FTS objects back to non-FTS objects to use in regression

prices_diff_fx = fts2mat(prices_diff_fx_ts);
prices_endhour = fts2mat(prices_endhour_ts);
prices_starthour = fts2mat(prices_starthour_ts);
prices_TWAP = fts2mat(prices_TWAP_ts);

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
prices_lagged_TWAP = prices_TWAP(1:size(prices_TWAP,1)-1,:);
prices_lagged_diff_fx = prices_diff_fx(1:size(prices_diff_fx,1)-1,:);
prices_lagged_starthour = prices_starthour(1:size(prices_starthour,1)-1,:);


prices_equities_adjclose = fts2mat(prices_equities_adjclose_ts);
prices_equities_open = fts2mat(prices_equities_open_ts);

%% Stepwise regression on Close of End-Hour FX and Adj Close of equities

if lagreg == 1
    %running the regression on lagged equity/unlagged FX
    reg1 = stepwiselm(prices_equities_adjclose_lagged1, ...
        prices_lagged_endhour(:,4), ...
        'VarNames',[names_equities 'EndHourPrices'])
%     reg2 = stepwiselm(prices_equities_adjclose_lagged1, ...
%        prices_lagged_TWAP(:,4))
%     reg3 = stepwiselm(prices_equities_adjclose_lagged1, ...
%             prices_lagged_diff_fx(:,4), ...
%             'VarNames',[names_equities '11amMinusTWAP'])
%     reg4 = stepwiselm(prices_equities_adjclose_lagged1, ...
%         prices_lagged_starthour(:,4), ...
%         'VarNames',[names_equities 'StartHourPrices'])

else
    %running regression on unlagged equity/unlagged FX
    reg1 = stepwiselm(prices_equities_adjclose, prices_endhour(:,4))
%     reg2 = stepwiselm(prices_equities_adjclose, prices_TWAP(:,4))
end    
    
%% Linear multivariate regression on Close of 11am FX and Adj Close 
% of equities

% fitlm(prices_equities, prices_11am(:,4))

