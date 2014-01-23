
clear all
close all

load prices_all_ts.mat

%% Do you want to take logs of equity prices?

takelogs = 0 %1 = yes; 0 = no

%% Do you want to run a lagged regression?

lagreg = 1 %1 = yes; 0 = no

%% Specify start-hour here

%Specify start and end hours in 24-hour format 
%(e.g., 5 for 5am, 17 for 5pm)
starthour = 11
startminute = 00
endhour = 12
endminute = 00

starthour = starthour+(startminute/60);
endhour = endhour+(endminute/60);

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




%% Converting FTS objects back to non-FTS objects to use in regression

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

if lagreg == 1
    %running the regression on lagged equity/unlagged FX
    reg1 = stepwiselm(prices_equities_adjclose_lagged1, ...
        prices_lagged_endhour(:,4), ...
        'VarNames',[names_equities 'EndHourPrices'])
%     reg4 = stepwiselm(prices_equities_adjclose_lagged1, ...
%         prices_lagged_starthour(:,4), ...
%         'VarNames',[names_equities 'StartHourPrices'])

else
    %running regression on unlagged equity/unlagged FX
    reg1 = stepwiselm(prices_equities_adjclose, prices_endhour(:,4))
end    
    
%% Linear multivariate regression on Close of 11am FX and Adj Close 
% of equities

% fitlm(prices_equities, prices_11am(:,4))

