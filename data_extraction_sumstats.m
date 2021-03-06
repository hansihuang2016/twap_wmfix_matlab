clear all
close all


%% Set working directory
workingdir = ...
    '/Users/arnavsheth/Documents/Research/AdamDuncan/Matlab/Data/';

fx_names = {'Open'; 'High'; 'Low'; 'Close'};

%% Specify start/end years, hours, increments and currency
    
currency = 'EURUSD' %ALLCAPS please

startyear = 2014
endyear = 2014

yearvector = startyear:endyear;

%Specify start and end hours in 24-hour format
    %(e.g., 5 for 5am, 17 for 5pm)
starthour = 7
startminute = 00
endhour = 17
endminute = 00

starthour = starthour+(startminute/60);
endhour = endhour+(endminute/60);

increments = (5/60)

hourvector = starthour:increments:endhour;

sumstats = zeros(4,size(hourvector,2)-1,size(yearvector,2));
large_mat_idx = zeros(size(yearvector,2)+1,1);
large_mat_idx(1) = 0;

for j = 1:size(yearvector,2)
    clear prices_TWAP

    %% Pull and extract initial data
    
    year = num2str(yearvector(j));

    %Get Downloaded Histdata.com FX Data
    prices_all = ...
        csvread(strcat(workingdir,...
        'DAT_ASCII_',currency,'_M1_',year,'_cleaned.csv'));
    
    %Create a time-series object using said data
    prices_all_ts = fints(datenum(prices_all(:,1:6)),...
        prices_all(:,7:10),fx_names);
    
    %Extract start and end dates
    datebounds = ftsbound(prices_all_ts,2);
    
    %Extracting the 11am price for each day in the year
    %This is our proxy for the 11am WM Fix
    prices_11am_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
        datebounds(2,:), [], 1, 'd',{'11:00'});
    
    %Convert to non-FTS object for later use
    prices_11am = fts2mat(prices_11am_ts,1);
    
    %Adding date+time columns for the 11am price matrix for completeness
    prices_11am = [datevec(prices_11am(:,1)) prices_11am(:,2:5)];
    
   
    %% Loop to get summary stats
    
    for i = 1:size(hourvector,2)-1
        %% Extract relevant data
        
        %Create vector of times and extract all prices from
        %starthour to starthour + 59mins daily
        dv = cellstr(datestr(...
            hourvector(i)/24:1/60/24:(hourvector(i+1)/24)-(1/60/24)));
        prices_hourly_ts = fetch(prices_all_ts, datebounds(1,:), [], ...
            datebounds(2,:), [], 1, 'd',dv);
        
        
        %% Calculating TWAP for the hour
        
        %%First: Convert to non-FTS object to be able to use the
        %accumarray() function
        prices_hourly = fts2mat(prices_hourly_ts,1);
        
        
        
        %%Second: Get a matrix of datevectors in the form [y, m, d, h, m, s] to
        %%use to use the unique() function and get a unique date for each time
        %%period. In other words, there will be a unique date for each set of
        %%data values from 10:00 to 10:59am, 9:00 to 9:59am, and so on...
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
        prices_TWAP_close = [uniquedates_hourly accumarray(subs_hourly, ...
            prices_hourly(:,5), [], @mean)];
        
        
        
        %%Fourth: We have to add a time column as we have an average per date
        %We are doing this because once we convert it back to a FTS object we will
        %need an identifier for each time
        if endminute == 0
            time = [(endhour-1) 59 00];
        else
            time = [(endhour-1) (endminute-1) 00];
        end
        time = repmat(time, size(prices_TWAP_close,1), 1);
        prices_TWAP = [prices_TWAP_open(:,1:3) time ...
            prices_TWAP_open(:,4)...
            prices_TWAP_high(:,4) prices_TWAP_low(:,4) ...
            prices_TWAP_close(:,4)];
        
        
        %Keeping the closing prices in a larger matrix
        large_mat_idx(j+1) = size(prices_TWAP_close(:,4),1);
        prices_TWAP_all(sum(large_mat_idx(1:j))+1:...
            sum(large_mat_idx(1:j+1)),i) = ...
            prices_TWAP_close(:,4);
        prices_TWAP_all(prices_TWAP_all == 0) = NaN;
        
        
        %Converting the hourly TWAPs to FTS object
        prices_TWAP_ts = fints(datenum(prices_TWAP(:,1:6)),...
            prices_TWAP(:,7:10),fx_names);
        
                
        %clearing up variables for mem reasons
        clear prices_all time prices_TWAP_open prices_TWAP_high...
            prices_TWAP_low prices_TWAP_close subs_hourly dv dates_hourly
        

        %Histogram of the hourly prices
        %     figure
        %     histfit(prices_TWAP(:,10))
        %     TWAP = [num2str(hourvector(i)) ' to' num2str(hourvector(i+1)) ' TWAP'];
        %     legend(TWAP)
        
        avg = nanmean(prices_TWAP(:,10));
        stdev = nanstd(prices_TWAP(:,10));
        skew = skewness(prices_TWAP(:,10));
        kurt = kurtosis(prices_TWAP(:,10));
        
        sumstats(:,i,j) = [avg;stdev;skew;kurt];
%         sumstats((j*4)-3:j*4,:) = [avg; stdev; skew; kurt];
        
        clear avg stdev skew kurt
    end
end

avg_all = nanmean(prices_TWAP_all);
stdev_all = nanstd(prices_TWAP_all);
skew_all = skewness(prices_TWAP_all);
kurt_all = kurtosis(prices_TWAP_all);
sumstats_all = [avg_all;stdev_all;skew_all;kurt_all];


% figure
% hold on
%     grid on
%     for jj = 1:size(yearvector,2)
%         plot3(sumstats(1,:,jj),repmat(yearvector(jj),1,size(hourvector,2)-1)1:size(hourvector,2)-1)
%     end   
% hold off


