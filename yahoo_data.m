
clear all
close all


%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';


%% Specify start-hour here

%Specify start and end hours in 24-hour format 
%(e.g., 5 for 5am, 17 for 5pm)
starthour = 10
startminute = 00
endhour = 11
endminute = 00

starthour = starthour+(startminute/60);
endhour = endhour+(endminute/60);

datebounds(1,:) = ['01/01/13'];
datebounds(2,:) = ['12/31/13'];

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

save prices_equities_open_ts.mat prices_equities_open_ts
save prices_equities_adjclose_ts.mat prices_equities_adjclose_ts
save returns_equities_adjclose_ts.mat returns_equities_adjclose_ts
save returns_equities_open_ts

