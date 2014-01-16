
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

%% Download equities data to get trading dates for equities. We will use 
%these to get our dateseries for which we get our FX data

tickers_americas = {'^MERV', '^MXX', '^BVSP', '^GSPTSE', '^GSPC', 'DOW',...
    '^RUT'};
names_americas = {'BuenosAires', 'Mexico', 'SaoPaulo', 'Toronto',...
    'SandP500', 'DowJones', 'Russell2000'};

% tickers_americas = {'^MXX', '^BVSP', '^GSPTSE', '^GSPC', 'DOW',...
%     '^RUT'};
% names_americas = {'Mexico', 'SaoPaulo', 'Toronto', 'SandP500',...
%     'DowJones', 'Russell2000'};

% tickers_americas = {'^BVSP', '^GSPTSE', '^GSPC', 'DOW'};
% names_americas = {'SaoPaulo', 'Toronto', 'SandP500', 'DowJones'};


tickers_europe = {'^ATX', '^BFX', '^FCHI', '^GDAXI', '^SSMI', '^FTSE'};
names_europe = {'Vienna', 'Brussels', 'Paris', 'Frankfurt', 'Swiss',...
    'London'};

tickers_asiapac = {'^AORD', '^HSI', '^BSESN', '^JKSE', '^KLSE',...
    '^N225', '^NZ50', '^STI', '^KS11', '^TWII'};
names_asiapac = {'Australia', 'HangSeng', 'Bombay', 'Jakarta', 'Malaysia',...
    'Nikkei225', 'NewZealand', 'Singapore', 'Korea', 'Taiwan'};

% tickers_europe = {'^ATX', '^BFX', '^SSMI', '^FTSE'};
% names_europe = {'Vienna', 'Brussels','Swiss','London'};


tickers_equities = [tickers_americas tickers_europe tickers_asiapac];
names_equities = [names_americas names_europe names_asiapac];

%connecting to yahoo! once to get basic info on one of the tickers
connect = yahoo
    equities_tmp = fetch(connect, tickers_equities(1),'Close',...
            datebounds(1,:),datebounds(2,:),'d'); 

    prices_equities_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
        names_equities(1));
close(connect)

%connecting to yahoo! again...
connect = yahoo
    %creating a loop (*sigh*) -- we have to do this because matlab does not
    %allow us to download more than one ticker at a time from yahoo...:P
    for ii = 1:size(tickers_equities,2)
        %creating a temporary vector for the yahoo download
        equities_tmp = fetch(connect, tickers_equities(ii),'Adj Close',...
            datebounds(1,:),datebounds(2,:),'d');
        
%         %sorting from earliest to latest (yahoo default is latest to
%         %earliest
%         equities_tmp = sortrows(equities_tmp,1);
        
        %converting to FTS object
        equities_tmp_ts = fints(equities_tmp(:,1),equities_tmp(:,2),...
            names_equities(ii));
        prices_equities_ts = merge(prices_equities_ts, equities_tmp_ts,...
            'DateSetMethod', 'intersection');
    end
close(connect)

returns_equities_ts = tick2ret(prices_equities_ts);

clear connect equities_tmp equities_tmp_ts ii equitiestime...
    datevec_equities

% save returns_equities_ts.mat returns_equities_ts
% save prices_equities_ts.mat prices_equities_ts


