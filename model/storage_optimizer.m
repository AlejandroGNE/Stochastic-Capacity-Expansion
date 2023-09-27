function [sOut] = storage_optimizer(sIn)

in.battery_capacity_MW = sIn.cap;
Case= sIn.storageCase;

switch Case
    case 0 
        %% default, no time clustering
            in.Clearing_price= sIn.ClrPrice;
        [out] = optimizer_local(in);
            Power=out.Power;
    case 1 
        %% use compression, run whole year in one go
            in.Clearing_price= sIn.ClrPrice(sIn.StorHours);
            in.horizon= sum(sIn.StorHours);
            in.interval_size= in.horizon;
        [out] = optimizer_local(in);
            Power= zeros(1,8760);
            Power(sIn.StorHours)= out.Power;
    case 2
        %% use compression, run one week at a time
%         k= 52;
%         df= sIn.ClrPrice;
%         a= 24*7; 
%         b= floor(length(df)/a);
%         df= df(1:a*b);
%         df_matrix= reshape(df,a,b);
%         Power= zeros(1,365*24);
        k= 52; a= 168; b= 52; 
        df_matrix= reshape(sIn.ClrPrice(1,1:168*52),a,b);
        Power= zeros(1,365*24);
         for i=1:k
            StorHours= logical(sIn.StorHours(:,i));
            in.Clearing_price= df_matrix(StorHours,i)';
            in.horizon= length(in.Clearing_price);
            in.interval_size= in.horizon;
            if in.horizon
        [out] = optimizer_local(in);
            df=zeros(1,a);
            df(StorHours)= out.Power; % allocate truncated into full week
            Power(1,i*a-a+1:i*a)= df;
            end
        end
        Power(1,end-24+1:end) = Power(1,a*b-24+1:a*b);
    case 3
        %% cluster weeks in groups, without compression
            k = length(sIn.repWeekIDs); % number of groups/rep. weeks
            df= sIn.ClrPrice;
            a = 24*7; % equal to a weeks' length in hours
            b = floor(size(df,2)/a); % number of weeks in a year
            df= df(1:a*b); % get rid of hours past the 52 weeks
            df_matrix= reshape(df,a,b);
            in.Clearing_price= zeros(1,a*k);
            for i=1:k
            in.Clearing_price(1,i*a-a+1:i*a)=...
            df_matrix(:,sIn.repWeekIDs(i))';
            end
            in.horizon = k*a;
            in.interval_size = a;
        [out] = optimizer_local(in);
            Power= out.Power;
            df= zeros(a,b);
            for i=1:k
            ids= sIn.kmeansMap(:,2) == i;
            df(:,ids') = repmat(Power(1,i*a-a+1:i*a)',1,sum(ids));
            end
            Power= zeros(1,365*24);
            Power(1,1:a*b)= reshape(df,1,a*b);
            Power(1,end-24+1:end) = Power(1,a*b-24+1:a*b);
    case 4
        %% cluster weeks into groups, with compression

        k = length(sIn.repWeekIDs); % number of groups/rep. weeks
        df= sIn.ClrPrice;
        a = 24*7; % weeks' length in hours
        b = floor(size(df,2)/a); % number of weeks in a year
        df= df(1:a*b); % get rid of hours past the 52 weeks
        df_matrix= reshape(df,a,b); % data in matrix form; columns are weeks

        % photw
        df1=zeros(1,a*k); % pre-allocate df for k repWeeks
        df2= zeros(a,b); % pre-allocate df for all 52 weeks
        for i=1:k
            StorHours= logical(sIn.StorHours(:,i));
            in.Clearing_price= df_matrix(StorHours,sIn.repWeekIDs(i))';
            in.horizon= length(in.Clearing_price);
            in.interval_size= in.horizon;
            if in.horizon
        [out] = optimizer_local(in);
            % photw to repWeeks allocation
            df=zeros(1,a);
            df(StorHours)= out.Power; % allocate truncated into full week
            df1(1,i*a-a+1:i*a) = df;
            ids= sIn.kmeansMap(:,2) == i;
            % repeat repWeeks to their groups to complete the 52
            df2(:,ids') = repmat(df1(1,i*a-a+1:i*a)',1,sum(ids));
            else
            end
        end
        Power= zeros(1,365*24);
        Power(1,1:a*b)= reshape(df2,1,a*b);
        Power(1,end-24+1:end) = Power(1,a*b-24+1:a*b);

%     case 2 % <- take average per hour-of-day for all weeks in season 
%         orig_prices = Clearing_price; % 1x8760 format
%         
%         % separate whole dataset in per-season df
%         cutoffs_spring = 78*24+1:170*24;            % 92 days
%         cutoffs_summer = 170*24+1:264*24;           % 94 days
%         cutoffs_fall = 264*24+1:354*24;             % 90 days
%         cutoffs_winter = [1:78*24,354*24+1:8760];   % 89 days
% %         seasonal_cutoffs = [78*24 170*24 264*24 354*24];
% 
%         % separate yearly data into each season so I can determine high and low
%         % weeks per season
%         prices_spring = orig_prices(cutoffs_spring);
%         prices_summer = orig_prices(cutoffs_summer);
%         prices_fall = orig_prices(cutoffs_fall);
%         prices_winter = orig_prices(cutoffs_winter);
% 
%         % for each season, convert to matrix form where 
%         % one_Row=oneWeek
%         % oneColumn=oneHour
%         df = prices_spring;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_avg = mean(df_matrix);
%         week1_avg = df_avg;
% 
%         df = prices_summer;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_avg = mean(df_matrix);
%         week2_avg = df_avg;
% 
%         df = prices_fall;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_avg = mean(df_matrix);
%         week3_avg = df_avg;
% 
%         df = prices_winter;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_avg = mean(df_matrix);
%         week4_avg = df_avg;
%         horizon = 4*168;
%         interval_size = 168;
%         Clearing_price = [week1_avg week2_avg week3_avg week4_avg];
%     case 3 % <- take median per hour-of-day for all weeks in season 
%         orig_prices = Clearing_price; % 1x8760 format
%         
%         % separate whole dataset in per-season df
%         cutoffs_spring = 78*24+1:170*24;            % 92 days
%         cutoffs_summer = 170*24+1:264*24;           % 94 days
%         cutoffs_fall = 264*24+1:354*24;             % 90 days
%         cutoffs_winter = [1:78*24,354*24+1:8760];   % 89 days
% 
%         % separate yearly data into each season so I can determine high and low
%         % weeks per season
%         prices_spring = orig_prices(cutoffs_spring);
%         prices_summer = orig_prices(cutoffs_summer);
%         prices_fall = orig_prices(cutoffs_fall);
%         prices_winter = orig_prices(cutoffs_winter);
% 
%         % for each season, convert to matrix form where 
%         % one_Row=oneWeek
%         % oneColumn=oneHour
%         df = prices_spring;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_med = median(df_matrix);
%         week1_med = df_med;
% 
%         df = prices_summer;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_med = median(df_matrix);
%         week2_med = df_med;
% 
%         df = prices_fall;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_med = median(df_matrix);
%         week3_med = df_med;
% 
%         df = prices_winter;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         df_med = median(df_matrix);
%         week4_med = df_med;
%         horizon = 4*168;
%         interval_size = 168;
%         Clearing_price = [week1_med week2_med week3_med week4_med];
%     case 4 % take peaks and valleys, interpolate the middle (INCOMPLETE)
%         orig_prices = Clearing_price; % 1x8760 format
%         
%         % separate whole dataset in per-season df
%         cutoffs_spring = 78*24+1:170*24;            % 92 days
%         cutoffs_summer = 170*24+1:264*24;           % 94 days
%         cutoffs_fall = 264*24+1:354*24;             % 90 days
%         cutoffs_winter = [1:78*24,354*24+1:8760];   % 89 days
% 
%         % separate yearly data into each season so I can determine high and low
%         % weeks per season
%         prices_spring = orig_prices(cutoffs_spring);
%         prices_summer = orig_prices(cutoffs_summer);
%         prices_fall = orig_prices(cutoffs_fall);
%         prices_winter = orig_prices(cutoffs_winter);
% 
%         % for summer & winter:
%         % 1) convert to matrix form where 
%         %   one_Row=oneWeek
%         %   oneColumn=oneHour
%         % 2) define a typical hour of the week with highest demand
%         % 3) choose the week with the highest value in 2)
% 
%         % summer
%         df = prices_summer;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         [~,I] = max(df_matrix,[],2);
%         topDemandHour= median(I);
%         slicedWeeks = df_matrix(:,topDemandHour);
%         [~,I] = max(slicedWeeks);
%         repWeekSummer = df_matrix(I,:);
%        % PENDING: DEFINE PEAKED AND NON-PEAKED SUMMER WEEKS TO MODEL
%        % STORAGE (maybe? maybe compare their storage outputs first; I'm
%        % thinking they would be the same)
% 
%         % winter
%         df = prices_winter;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
%         [~,I] = max(df_matrix,[],2);
%         topDemandHour= median(I);
%         slicedWeeks = df_matrix(:,topDemandHour);
%         [~,I] = max(slicedWeeks);
%         repWeekWinter = df_matrix(I,:);
% 
%         % for spring and fall:
%         % 1) convert to matrix form where 
%         %   one_Row=oneWeek
%         %   oneColumn=oneHour
%         % 2) define a typical hour of the week with highest demand
%         % 3) choose the week with the *median* value in 2)        
%         
%         % spring
%         df = prices_spring;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
% %         [~,I] = max(df_matrix,[],2);
% %         topDemandHour= median(I);
% %         slicedWeeks= df_matrix(:,topDemandHour);
% %         slicedWeeks= sort(slicedWeeks);
% %         middle = length(slicedWeeks)/2;
% %         add = rem(middle,2);
% %         I = middle+add-1;
% %         repWeekSpring = df_matrix(I,:);
%         repWeekSpring = df_matrix(6,:);
%         
% 
%         % fall
%         df = prices_fall;
%         columnWidth = 24*7; % equal to a weeks' length in hours
%         columnHeight = floor(sizee(df,2)/columnWidth); % equal to the number of weeks in a year
%         df_matrix = reshape(...
%             df(1:columnHeight*columnWidth),...
%             columnHeight, columnWidth);
% %         [~,I] = max(df_matrix,[],2);
% %         topDemandHour= median(I);
% %         slicedWeeks= df_matrix(:,topDemandHour);
% %         slicedWeeks= sort(slicedWeeks);
% %         middle = length(slicedWeeks)/2;
% %         add = rem(middle,2);
% %         I = middle+add-1;
%         repWeekFall = df_matrix(6,:);
% 
% %         % plot all weeks overlaid
% %         figure
% %         hold on
% %         for row=1:size(df_matrix,1)
% % %             figure
% %             plot(df_matrix(row,:))
% %         end
% %         hold off
% % 
% %         % plot all hours in a season as time series
% %         plot(df)
% %         xline(168:168:168*13)
% % 
% %         shoulderAll = [repWeekSpring repWeekFall];
% %         df= shoulderAll;
% 
%         horizon = 4*168;
%         interval_size = 168;
%         Clearing_price =...
%             [repWeekSpring repWeekSummer repWeekFall repWeekWinter];
%     case 5 % two repWeeks, high and low storage output 
% %         orig_prices = Clearing_price; % 1x8760 format
%         orig_prices = Power; % 1x8760 format
% 
%         % separate whole dataset in the two blocks
%         cutoffsHighStor= 119*24+1:287*24; % 168 days (24 weeks)
%         cutoffsLowStor= [1:119*24 287*24+1:8760]; % 197 days (28 weeks)
%         % 8760 hours have 52 weeks and one day
% 
%         pricesLowStor= orig_prices(cutoffsLowStor);
%         pricesHighStor= orig_prices(cutoffsHighStor);
% 
%         % for high storage period:
%         % 1) reshape data, oneColumn=oneWeek and oneRow=OneHour
%         % 2) produce the expected forecast that would result of choosing
%         % each week as representative
%         % 3) get the error of each forecast: e = f - actual
%         % 4) get the std for each vector of errors
%         % 5) choose week that minimizes std
% 
%         % NOTE: THIS SHOULD BE USED USING A FIRST ANNUAL ESTIMATION OF
%         % STORAGE OUTPUT, NOT USING PRICE.  THE FOLLOWING IS JUST AN
%         % EXAMPLE MADE TO DEVELOP THE APPROACH USING PRICES DATA AS INPUT.
% 
%         df = pricesHighStor;
%         columnHeight = 24*7; % equal to a weeks' length in hours
%         columnWidth = floor(sizee(df,2)/columnHeight); % equal to the number of weeks in the period
%         df = df(1:columnHeight*columnWidth);
%         df_matrix = reshape(df,columnHeight,columnWidth);
%         df_forecast = zeros(length(df),columnWidth);
%         df_e = zeros(1,columnWidth);
%         df_std = zeros(1,columnWidth);
%         df_total = sum(df_matrix,"all");
%         for i=1:columnWidth
%             df_forecast(:,i) = repmat(df_matrix(:,i),columnWidth,1);
%             df_forecast_total= sum(df_forecast(:,i));
%             df_e(1,i) = ((df_forecast_total - df_total)/df_total)*100;
% %             df_e(:,i)= ((df_forecast(:,i)'-df)./df)*100;
% %             df_std(1,i) = sum(df_e(:,i));%std(df_e(:,i));
%             df_std(1,i) = df_e(1,i);
%         end
%         df_std= abs(df_std);
%         [~,I]=min(df_std);
%         repWeekHigh= df_matrix(:,I)';
% 
% 
%         df = pricesLowStor;
%         columnHeight = 24*7; % equal to a weeks' length in hours
%         columnWidth = floor(sizee(df,2)/columnHeight); % equal to the number of weeks in the period
%         df = df(1:columnHeight*columnWidth);
%         df_matrix = reshape(df,columnHeight,columnWidth);
%         df_forecast = zeros(length(df),columnWidth);
%         df_e = zeros(1,columnWidth);
%         df_std = zeros(1,columnWidth);
%         df_total = sum(df_matrix,"all");
%         for i=1:columnWidth
%             df_forecast(:,i) = repmat(df_matrix(:,i),columnWidth,1);
%             df_forecast_total= sum(df_forecast(:,i));
%             df_e(1,i) = ((df_forecast_total - df_total)/df_total)*100;
% %             df_e(:,i)= ((df_forecast(:,i)'-df)./df)*100;
% %             df_std(1,i) = sum(df_e(:,i));%std(df_e(:,i));
%             df_std(1,i) = df_e(1,i);
%         end
%         df_std= abs(df_std);
%         [~,I]=min(df_std);
%         repWeekLow= df_matrix(:,I)';
% 
%         horizon = 2*168;
%         interval_size = 168;
%         Clearing_price =...
%             [repWeekLow repWeekHigh]; 
%     case 6 % systematically choose best weeks
%         k = 4; % number of groups/rep. weeks
% 
%         df = Power;
%         columnHeight = 24*7; % equal to a weeks' length in hours
%         columnWidth = floor(sizee(df,2)/columnHeight); % equal to the number of weeks in the period
%         df= df(1:columnHeight*columnWidth); % get rid of hours past the 52 weeks
%         df_matrix= reshape(df,columnHeight,columnWidth);
%         df_weekly_totals= sum(df_matrix)';
%         df_weekly_totals(:,2)= 1:columnWidth; % add week of the year identifier
% %         df_sorted= sortrows(df_weekly_totals,'descend');
%         idx= kmeans(df_weekly_totals(:,1),k);
%         df_weekly_totals(:,3)= idx;
%         
%         % now that we know how weeks group together, pick one rep. week per
%         % group; before, we did this by forecasting each group's throughput
%         % using each week as the chosen one, then going for the week with
%         % the least error; that seemed to lead to the week with the most
%         % average throughput (or the most "middle" week), so we could just
%         % do that this time
%         
%         % separate data into groups
%         repWeekIds= zeros(k,1);
%         vals= zeros(k,1);
%         repWeeks= zeros(1,columnHeight*k);
%         for i=1:k
%             id= df_weekly_totals(:,3) == i;
%             df= df_weekly_totals(id,1:2);
%             avg= repmat(mean(df(:,1)),sizee(df,1),1);
%             dev= abs(df(:,1)-avg);
%             [~,I]= min(dev);
%             repWeekIds(i)= df(I,2);
%             vals(i)= df(I,1);
%             repWeeks(1,i*columnHeight-columnHeight+1:i*columnHeight)=...
%                 df_matrix(:,repWeekIds(i))';
%         end
% 
%         % done! now, here we determined what weeks fit best based on Power,
%         % because that's running on 8760... how do we run 8760 a limited
%         % number of times but enough so we can systematically determine
%         % those representative weeks? when should we run the first 8760 
%         % storage optimizer? and once we have that data for the rep. weeks,
%         % what code will be here, in its place, to be run with the regular
%         % dispatch every period for the price maker?
% 
%         % option A: run storage optimizer in initializer, retrieve indices
%         % for rep. weeks (for a given number of rep. weeks); for the
%         % accounting model part always run the clustered version
% 
%         % option B: calibrate the rep. weeks indices for each grid before
%         % running any scenario, assuming these stay constant
% 
%         % let's go with option B
%         
%         % not really an option if we want the flexibility of changing how
%         % many rep weeks we want to use
end

NetLoad_after_storage = (sIn.Load2015)'+Power;
profit = -(Power*(sIn.ClrPrice)');
sOut.Loadstor= NetLoad_after_storage;
sOut.Power = Power;
sOut.profit = profit;
end




function [out] = optimizer_local(in)
% This function finds the optimal operation of the
% battery. We cannot optimize for a single year at a
% time, so the year must be partitioned into smaller
% intervals. Subsequent intervals are dependent on the
% previous one, so they are run serially. Because
% optimizing an interval is not identical to optimizing
% the full year, we add enough extra time to each
% interval, such that it emulates the operation of a
% full year for the portion of th interval. This extra
% time is called overlap, and the joining of an
% interval with its corresponding overlap is called
% compartment. The overlap is discarded at the end of
% each compartment / start of the next compartment.
% This is know in the literature as "rolling" the
% optimizations.
% 2021/10/11 update:
% We found that using a smaller interval size leads to
% huge gains in runtime. We also found that using no
% overlap is virtually the same to using many small
% overlaps; all results tested were within 1% of the
% best-case scenario, which is obtained by using 8760
% as the interval size.

%% Import inputs and set parameters 
% essential inputs
battery_capacity_MW= in.battery_capacity_MW;
Clearing_price= in.Clearing_price;

% parameters
try 
    horizon= in.horizon;
catch
    horizon=8760;
end

try
    overlap_size= in.overlap_size;
catch
    overlap_size= 0;
end

try
    interval_size= in.interval_size;
catch
    interval_size= 168;
end

%% Define the size of the compartments
% tic

% compartment size in hours
compartment_size = interval_size + overlap_size;

% position where intervals are stitched together
starting_hour = 1; 

% number of times the optimization will be run
number_of_compartments = ceil(horizon/interval_size); 
total_hours_run =...
    number_of_compartments*interval_size...
    +overlap_size;
hours_past_8760 = total_hours_run-horizon;

%% Properties of the battery

% Roundtrip efficiency
RTE = 0.80;
round_trip_efficiency = sqrt(RTE);

% Number of hours it takes to attain full charge when 
% charged at max Energy capacity
time_to_full_charge = 4; 

% Max charge or Energy capacity of the battery
battery_energy_capacity_MWh =...
    battery_capacity_MW*time_to_full_charge; 

initial_SOC = 0;

% Initial state of charge (how much charge batteries 
% had "brand new", before operations begin)
remaining_SOC_previous_compartment_MWh = initial_SOC; 

%% Extend variables to the size of all compartments
% electricity_demand(horizon+1:total_hours_run) =...
%     electricity_demand(1:hours_past_8760); 

% rename clearing_price to preserve original data
clearing_price = Clearing_price; 

clearing_price(horizon+1:total_hours_run) =...
    Clearing_price(1:hours_past_8760); 

%% Matrix to set constraints

charge = ones(compartment_size,1)...
    .* round_trip_efficiency; 

discharge = -ones(compartment_size,1)...
    ./ round_trip_efficiency;

SOC = zeros(2*compartment_size,2*compartment_size);

% for each hour of the year included in the compartment
for counter = 1:compartment_size 

    % copy "charge" into each cell of the lower 
    % triangle of the diagonal of cutting the left half
    % of SOC matrix from left-top to bottom-right
    SOC(counter,1:counter) = charge(counter); 
    
    % copy "discharge" into each cell of the lower 
    % triangle of the diagonal of cutting the right 
    % half of SOC matrix from left-top to bottom-right
    SOC(counter,(compartment_size+1):...
        counter+compartment_size) = discharge(counter); 
end

% extend current SOC matrix by putting a negative SOC 
% below it
SOC(compartment_size+1:...
    2*compartment_size,1:2*compartment_size) = ...
    -SOC(1:compartment_size,1:(2*compartment_size)); 

% this matrix is helpful to set the constraints on 
% state of charge for each hour based on the amount of 
% charge/discharge of storage in the previous hour

%% Pre-allocate variables
Maxcharge = zeros(2*compartment_size,1);
compartment_clearing_prices =...
    zeros(2*compartment_size,1);
Energy = zeros(1,total_hours_run);
Power = zeros(1,total_hours_run);

% uncomment to store charge and discharge variables
%         X1 = zeros(1,total_hours_run); 
%         X2 = zeros(1,total_hours_run);

% NetLoad_after_storage = zeros(1,total_hours_run);
optimal_state_of_charge = zeros(1,total_hours_run);

% lower bound = zero i.e. do nothing
lb = zeros(2*compartment_size,1); 

% upper bound is maximum charging/discharging per hour
ub = repmat(battery_capacity_MW,2*compartment_size,1); 

%% Run serial optimizations of compartments

for iterations = 1:number_of_compartments
    
    % interval: number of hours for which storage 
        % operation is optimized at a time (hours)

    % overlap: number of hours added to interval for 
        % optimization but discarded later

    % compartment: interval + overlap
    
    current_compartment = ...
        starting_hour...
        :starting_hour+compartment_size-1;

    % initial charge is predefined for 1st iteration 
    % and estimated for 2nd iteration onwards, in MWh
    remaining_SOC_previous_compartment_MWh =...
        repmat(...
        remaining_SOC_previous_compartment_MWh,...
        compartment_size,1); 
    
    % How much energy can the battery take in? 
    % Can't buy/charge more energy than the total 
    % capacity of the battery, in MWh
    % This is the RHS constraint of the St contraint, 
    % Smax
    Maxcharge(1:compartment_size,1) =...
        repmat(...
        battery_energy_capacity_MWh,...
        compartment_size,1)...
        - remaining_SOC_previous_compartment_MWh; 
    
    % How much energy can the battery give out? 
    % Can't sell/discharge more energy than what's 
    % already in the battery, in MWh
    % concatenating constraints for charge and 
    % discharge limits for the SOC matrix
    Maxcharge(...
        compartment_size+1:2*compartment_size,1) =...
        remaining_SOC_previous_compartment_MWh; 
    
    % clearing prices are positive for buying/charging
    compartment_clearing_prices(...
        1:compartment_size,1) =...
        clearing_price(current_compartment)'; 

    % and negative for selling/discharging
    compartment_clearing_prices(...
        compartment_size+1:2*compartment_size,1) =...
        -compartment_clearing_prices(...
        1:compartment_size,1); 

    % linprog internally reshapes f to column vector
    f = compartment_clearing_prices; 
    opts.Algorithm='dual-simplex';
    opts.Display='off';
    
    % Linprog problem is structured in two parts, one
    % for charging (upper half), one for discharging
    % (lower half).
    % x: decision variable, energy in and out, MWh
    % f: clearing prices, $/MWh
    % SOC: constraints
    
    [x] = linprog(f,SOC,Maxcharge,[],[],lb,ub,[],opts);

    if isempty(x)
        x = zeros(2*compartment_size,1);
    end

    % for charge variables
    x1 = x(1:compartment_size); 

    % for discharge variables
    x2 = x((compartment_size+1):(2*compartment_size)); 

    % Energy in/out each hour, units are MWh after 
    % considering the losses
    Energy(current_compartment) =...
        charge.*x1+(discharge.*x2); 

    % Hourly energy flows before affecting efficiency
    Power(current_compartment) = x1-x2; 

    % the SOC of the current hour is the sum of Energy 
    % up to the current hour plus the initial SOC
    optimal_state_of_charge(current_compartment) =...
        cumsum(Energy(current_compartment)) ...
        +remaining_SOC_previous_compartment_MWh(1); 
    
    % prepare for next iteration:
    % going from the current to the next interval, 
    % subtracting the overlap
    starting_hour = starting_hour + interval_size; 

    remaining_SOC_previous_compartment_MWh =...
        optimal_state_of_charge(starting_hour-1);
end

%% Export variables
% crop variables in case last compartment size > horizon
if hours_past_8760>0
    out.Power = Power(1:horizon);
else
    out.Power=Power;
end

end