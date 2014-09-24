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



%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';

%% Pull and extract FX data

currency = 'EURUSD' %ALLCAPS please

startyear = 2007
endyear = 2009

yearvector = startyear:endyear;


prices_all = ...
    csvread(strcat(workingdir,...
        'DAT_ASCII_',currency,'_M1_',num2str(startyear),'_cleaned.csv'));
for jj = 2:size(yearvector,2)
    %% Pull and extract initial data
    
    year = num2str(yearvector(jj));

    %Get Downloaded Histdata.com FX Data
    tmp = ...
        csvread(strcat(workingdir,...
        'DAT_ASCII_',currency,'_M1_',year,'_cleaned.csv'));
    
    prices_all = vertcat(prices_all, tmp);  
    clear tmp
end
clear tmp

fx_names = {'Open'; 'High'; 'Low'; 'Close'};

%Create a time-series object using said data
prices_all_ts = fints(datenum(prices_all(:,1:6)),prices_all(:,7:10),fx_names);

save prices_all_ts.mat prices_all_ts

%Extract start and end dates for equities extraction
datebounds = ftsbound(prices_all_ts,2);

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

save prices_equities_open_ts.mat prices_equities_open_ts
save prices_equities_adjclose_ts.mat prices_equities_adjclose_ts





