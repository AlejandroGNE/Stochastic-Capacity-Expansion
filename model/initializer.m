function [s] = initializer(s)  
% Takes modeling choices and assumptions and prepares 
% them to be used by:
% 1) the accounting module (cashflows.m), if estimating
% the system costs for a given set of additions and 
% retirements, or
% 2) the solver, if searching for the optimal set of 
% additions and retirements.
%% Unpack inputs structure

% The general structure for this section is:
% try: read in modeling choices and assumptions
% catch: where no choices were made, use a default

try
    Renewable_standard=s.Renewable_standard;  
catch
    Renewable_standard=0;
    s.Renewable_standard=Renewable_standard;
end

try
    Storage_standard=s.Storage_standard;  
catch
    Storage_standard=0;
    s.Storage_standard=Storage_standard;
end

try 
    Renewable_subsidies=s.Renewable_subsidies;
catch
    Renewable_subsidies=0; 
    s.Renewable_subsidies=Renewable_subsidies; 
end

try 
    Storage_subsidies=s.Storage_subsidies;
catch 
    Storage_subsidies=0; 
    s.Storage_subsidies=Storage_subsidies; 
end

try 
    region=s.region; 
catch 
    region = 1; 
    s.region=region; 
end

try 
    grid_data=s.grid_data;
catch
    grid_data=2022; 
    s.grid_data=grid_data; 
end

try 
    new_plants_data=s.new_plants_data;
catch
    new_plants_data=0; 
    s.new_plants_data=new_plants_data; 
end

try
    ramping=s.ramping;
catch
    s.ramping = 0;
end

try 
    storage=s.storage; 
catch
    storage=0; s.storage=storage; 
end

try 
    monte=s.monte; %#ok<*NASGU> 
catch
    monte=1; s.monte=monte; 
end

try 
    mutation_shrink=s.mutation_shrink;
catch
    mutation_shrink=0.95; 
    s.mutation_shrink=mutation_shrink;
end

try 
    interest=s.interest; 
catch
    interest=0.05; 
    s.interest=interest; 
end

try 
    Periods=s.Periods; 
catch
    Periods=7; 
    s.Periods=Periods; 
end

% PeriodSize = Number of years per period
% To reduce dimensionality, this model estimates
% additions and retirements every period, not every
% year, and each period can be any size, for example, 
% two, three, five years long. To run for 2030-2050, if
% a PeriodSize of 5 years is defined, the model will
% only run for 2030, 2035, 2040, and 2045, four times,
% as opposed to running twenty one times.
try 
    PeriodSize=s.PeriodSize; 
catch
    PeriodSize=5; 
    s.PeriodSize=PeriodSize; 
end

try 
    fixed_costs_options=s.fixed_costs_options;  
catch
    fixed_costs_options=1; 
    s.fixed_costs_options=fixed_costs_options; 
end

try 
    mutation_scale=s.mutation_scale;
catch
    mutation_scale=0.05; 
    s.mutation_scale=mutation_scale; 
end

try 
    Mesh_Tolerance=s.Mesh_Tolerance;
catch
    Mesh_Tolerance=1e-6; 
    s.Mesh_Tolerance=Mesh_Tolerance; 
end

try 
    FunctionTolerance_GA=s.FunctionTolerance_GA;
catch 
    FunctionTolerance_GA=1e-7;
    s.FunctionTolerance_GA=FunctionTolerance_GA;
end

try 
    PopulationSize=s.PopulationSize;
catch 
    PopulationSize=300;
    s.PopulationSize=PopulationSize;
end

try
    FunctionTolerance_PS=s.FunctionTolerance_PS;
catch
    FunctionTolerance_PS=1e-6;
    s.FunctionTolerance_PS=FunctionTolerance_PS;
end

try 
    MaxGenerations=s.MaxGenerations;
catch
    MaxGenerations=600; 
    s.MaxGenerations=MaxGenerations; 
end

try 
    learning=s.learning;
catch
    learning=2; 
        s.learning=learning; 
end

try 
    plot_solvers=s.plot_solvers;
catch
    plot_solvers=0; 
    s.plot_solvers=plot_solvers; 
end

try 
    increase_highest_demand_hour=...
        s.increase_highest_demand_hour;
catch
    increase_highest_demand_hour=1; 
    s.increase_highest_demand_hour=...
        increase_highest_demand_hour; 
end

try
    unmet_demand_penalization=...
        s.unmet_demand_penalization;
catch
    unmet_demand_penalization=1e6;
    s.unmet_demand_penalization=...
        unmet_demand_penalization;
end

try
    report = s.report;
catch
    s.report = 0;
end

try
    penalty_approach = s.penalty_approach;
catch
    s.penalty_approach = 2;
end

try
    square_pop = s.square_pop;
catch
    s.square_pop = 1;
end

try
    alpha = s.alpha;
catch
    if region == 1
        s.alpha = 21;
    elseif region == 2
        s.alpha = 32;
    end
end

try
    cluster = s.time_cluster;
catch
    s.time_cluster = 0;
    cluster= 0;
end

try
    compression = s.compression;
catch
    s.compression = 0;
    compression= 0;
end

try
    devDebug= s.devDebug;
catch
    s.devDebug=0;
    devDebug=0;
end

try
    slurmid = s.slurmid;
    if s.cpuspertask > 36
        s.slurmid = s.slurmid+1;
    end
    rng(s.slurmid);
catch
    rng shuffle
%     s.slurmid = randi(1e9);
end

% Make sure random number engines are always different
% rng(s.slurmid);

%% Load region data

% $65 for each MWh out of RPS compliance
s.RPS_non_compliance_penalty = 65;

% $65 for each MWh out of SPS compliance
s.SPS_non_compliance_penalty = 65;

% Over-periods for proper cost accounting
s.BPeriods = Periods + 4; 

s.Year = 2016; % Starting year

% Load power plants data, renewable resource data, and
% demand data for a given region and year
[plant,Var_Energy,Load2] =...
    Region_data(region,grid_data);
s.Var_Energy = Var_Energy;

% multiply highest-demand hour of the year by 1.15 to
% see if model stops retiring storage in the first
% period
if increase_highest_demand_hour == 1
    [~,I] = maxk(Load2',1);
    for i = 1:Periods
        Load2(i,I(i)) = Load2(i,I(i)) * 1.15;
    end
    s.I = I;
end
s.Load2 = Load2;

plant.LCOE_dolMWh(1) = 0; % pre-allocate new column for LCOE data

if learning == 0
    % Use exogenous learning model with old data
    load('Global_pop','P0'); 
elseif learning == 1
    % Use  exogenous learning model with new data
    load('Global_pop',...
        'icap_additions_world_IEO_2021_v0'); 
    P0 = icap_additions_world_IEO_2021_v0;
elseif learning > 1
    % Use endogenous learning model; we load P0 just so
    % we have something to pass on to cashflows
    load('Global_pop','P0'); 
end
s.P0 = P0;

% Cost and performance assumptions for new power plants
[Cap_small,Cap_small_storage] =...
    New_power_plants_data(new_plants_data);

% If not modeling storage
if storage == 0
    % Delete storage plants from database and
    plant(plant.Category == 'STORAGE',:) = [];
    % Import cost data without storage
    Cap_cost = Cap_small;
else
    % Otherwise, import cost data with storage
    Cap_cost = Cap_small_storage;
end
s.plant = plant;

% Load EPA FOM data
load("EPA_v6_table_4_9_FOM.mat","EPA_v6_table_4_9")
s.EPA_FOM = EPA_v6_table_4_9;

Cap_cost.Tech = categorical(Cap_cost.Tech);
if fixed_costs_options == 0
    plant.OnM_dolkW(...
        plant.Category == "NUCLEAR") =...
        Cap_cost.OnM_dolkW(...
        Cap_cost.Tech == "NUCLEAR");
end
Cap_cost.Tech = string(Cap_cost.Tech);

plant.Emissions_control_all =...
    zeros(size(plant.PlantType));
% Append row for gas_CC retirements
Cap_cost(end+1,:) = Cap_cost(2,:);

% Number of new technologies that can be added
num_tech = size(Cap_cost,1); s.num_tech = num_tech;

try
    smooth = s.smooth;
catch
    s.smooth = ones(num_tech,1);
end

% % Add 1 to account for split in gas_CC R&As
% num_tech = num_tech + 1; s.num_tech = num_tech;

% Names of new power plants
s.name1 = strcat(Cap_cost.Tech,'1'); 

% Number of decisions variables for solver
% Decision variables are additions or retirements of
% each technology in each period
nvars = Periods*num_tech; s.nvars = nvars;

%% Subsidies

% The user has the option to add subsidies for solar,
% wind or storage technologies.

subsidies = zeros(size(Cap_cost.Capital_cost_dolKW));
% subsidies(end) = [];

if Renewable_subsidies ~= 0
    subsidies([5 6],1) = repelem(Renewable_subsidies,2);
end
if Storage_subsidies ~= 0
    subsidies(8,1) = Storage_subsidies;
end
s.subsidies = subsidies;

%% Set maximum additions and retirements  

% Get installed capacity to set limit on retirements
installed_capacity = get_installed_capacity(plant);
% Cap_cost.Fuel = categorical(Cap_cost.Fuel);
Cap_cost.Tech = categorical(Cap_cost.Tech);
for fuel = 1:length(Cap_cost.Tech)
    try
        Cap_cost.Retire_cap(fuel) =...
            -installed_capacity.installed_capacity_MW(...
            installed_capacity.technologies ==...
            Cap_cost.Tech(fuel));
    catch
        Cap_cost.Retire_cap(fuel) = 0;
    end
end
Cap_cost.Tech = cellstr(Cap_cost.Tech);

% Let model install up to 20 GW of something per year
% (100 GW per 5-year period)
Cap_cost.max_cap(:) = 100000;

% For renewables, let model install up to 30 GW/year
% 150 GW/period
if storage == 0
    Cap_cost.max_cap([5 6]) =   150000; 
else
    Cap_cost.max_cap([5 6]) = 150000;
    % For storage, same as renewables
    Cap_cost.max_cap(8) = 150000;
end

% Add 1 MW to existing capacity to allow model to
% remove all plants without splitting retirements
% between periods
Cap_cost.Retire_cap = Cap_cost.Retire_cap -1;

% first gas_CC can only add plants
Cap_cost.Retire_cap(2) = 0;

% second gas_CC can only retire
Cap_cost.max_cap(end) = 0;

% Ban gas_CC retirements to reduce search space
% Model can still retire gas through gas_CT
% (since retiring gas CC and CT makes no difference 
% within the model; pricier plants can be retired first
% regardless of them being CC or CT)

% Cap_cost.Retire_cap(2) = 0;

% 8/20/2022 update: removed ban to gas_CC retirements,
% since we considered that gas_CC and gas_CT represent
% different types of resources and the model should be
% able to decide which resource to add/keep/retire

% 100% renewable constraints
try
if s.decarbGoal == 100
    Cap_cost.max_cap([1 2 3]) = 0;
end
catch
end

s.Cap_cost = Cap_cost;

%% Do not retire more capacity than exists
% A = zeros(num_tech,num_tech*Periods);
% 
% % each technology = one row in A
% % for each technology,
% % the sum of its capacity additions and retirements 
% % throughout all periods...
% for tech_period = 1:num_tech
%     A(tech_period,tech_period:num_tech:num_tech*Periods) =...
%         -ones(...
%         size(A(tech_period,tech_period:num_tech:num_tech*Periods))); 
% end
% 
% % ... cannot exceed its existing capacity.
% % in other words, the solver cannot pick retirement
% % values that retire more capacity than exists
% b = zeros(num_tech,1);
% b(1:num_tech) = -Cap_cost.Retire_cap;
% s.A = A;
% s.b = b;

%% Implement SPS
if s.Storage_standard > 0 && s.storage > 0

    storage_power = sum(plant.Capacity_MW(...
        plant.Category == 'STORAGE'));
    storage_power = repmat(storage_power,1,Periods);

    SPS_target = zeros(1,Periods);

    % storage standard applied to last period
    SPS_target(1,Periods) = max(Load2(7,:))*...
        s.Storage_standard;

    % linearly interpolate standard to intermediate
    % periods
    SPS_target(1,2:Periods-1) = interp1(...
        [1 Periods],...
        [SPS_target(1) SPS_target(end)],...
        2:Periods-1);

    % add enough storage to comply with mandate
    add_storage = zeros(1,Periods);
    for i=2:Periods
        gap = SPS_target(i)-storage_power(i-1);
        storage_power(i) = storage_power(i-1) + gap;
        add_storage(i) = gap;
    end

%     % check if system complies every period
%     SPS_non_compliance_MW =...
%         SPS_target - storage_power;
%     SPS_non_compliance_MW(SPS_non_compliance_MW<0) = 0;
% 
%     storage_row = Cap_cost.Tech == "Storage";
% 
%     % Penalty is the capital cost of storage itself
%     s.SPS_non_compliance_penalty =... % dol/MW
%         Cap_cost.Capital_cost_dolKW(storage_row)...
%         * 1 ... % <-- scale up penalty here
%         * 1000; % dol/kW * 1000 kW/MW = dol/MW
%     
%     Cost_of_SPS_non_compliance=...
%         SPS_non_compliance_MW.* ...
%         repelem(s.SPS_non_compliance_penalty,Periods);
end

% apply bounds defined above to decision variables in
% a solver-ready format
max_cap=repmat(Cap_cost.max_cap,1,Periods);
low_cap=repmat(Cap_cost.Retire_cap,1,Periods);
if s.Storage_standard > 0 && s.storage > 0
    low_cap(8,:) = add_storage;
end
lb=reshape(low_cap,1,Periods*num_tech); s.lb = lb;
ub=reshape(max_cap,1,Periods*num_tech); s.ub = ub;

% force nuclear by 2030
try
    if s.forceNuclear
        nucTech= 7;
        gridSize= max(Load2,[],'all');
        reICAP= plant.Capacity_MW(plant.Category == 'BIOMASS')...
            + plant.Capacity_MW(plant.Category == 'HYDRO');
        nucSize= gridSize - reICAP;
        lb1= reshape(lb,[],Periods);
        ub1= reshape(ub,[],Periods);
        lb1(nucTech,s.decarbYear) = nucSize;
        ub1(nucTech,s.decarbYear) = nucSize;
        lb1(5:end-1,:)= 0;
        s.lb= reshape(lb1,1,[]);
        s.ub= reshape(ub1,1,[]);
    end
catch
end

% hard-code extra retirements
try
    if s.hardCodeFossilRets
        theFiles = dir('**/*.mat');
        for k= 1:length(theFiles)
            if startsWith(theFiles(k).name,s.baseOptimal)
                baseFileName= theFiles(k).name;
                fullFileName=...
                    fullfile(theFiles(k).folder,baseFileName);
                baseOutputs= load(fullFileName,'outputs');
                baseS= load(fullFileName,'s');
                baseOptimal= baseOutputs.outputs.Optimal_DeltaMW;
            end
        end
        % for fossils, copy all retirements from
        % optimal buildout, hard code them
        baseOptimal1= reshape(baseOptimal,[],Periods);
        
        if s.storage
        reshaped= baseOptimal1;
        reshaped2= reshaped; % make it a matrix
        reshaped2(9,:)=reshaped(8,:); % adapt to gas_ret format
        reshaped2(8,:)=zeros(1,size(reshaped2,2)); % add row for storage
%         reshaped3= reshape(reshaped2,1,[]); % make it vector again
        baseOptimal1= reshaped2;
        end
        
        % set goals for target year
        targetPeriod= s.decarbYear;
        goal= Cap_cost.Retire_cap;
        goal(4:end-1)=0;
        fossils= goal<0;
        remainingFossils= cumsum(baseOptimal1(fossils,1:targetPeriod)')';
        remainingFossils= remainingFossils(:,end);
        difference= goal(fossils) - remainingFossils;
            % <0 means we need to hard code that retirement
            % 0 means we can just hard code zeros after this
        difference(difference>=0) = 0;
        baseOptimal1(fossils,targetPeriod) =...
            baseOptimal1(fossils,targetPeriod) + difference;
        baseOptimal1(fossils,targetPeriod+1:end)= 0;
        try
        if s.forbidBiomass
            lb1(4,:) = zeros(1,Periods);
            ub1(4,:) = zeros(1,Periods);
        end
        catch
        end
        
        lb1= reshape(lb,[],Periods);
        ub1= reshape(ub,[],Periods);
        lb1(fossils,:) = baseOptimal1(fossils,:);
        ub1(fossils,:) = baseOptimal1(fossils,:);
        lb1(5:end-1,:)= 0;
        
        try
            if s.onlyWindSolar
                nonWindSolar= [4 7];
                ub1(nonWindSolar,:)= 0;
            end
        catch
        end
        
        s.lb= reshape(lb1,1,[]);
        s.ub= reshape(ub1,1,[]);
    end
catch
end

%% Apply portfolio standards policy

% Use peak load in last period as a reference to 
% constrain portfolio standards applied to capacity 
% For example: "Storage must be 50% of the total grid
% capacity in 2050"
peak_load_in_last_period = max(Load2,[],'all');

% Use total generation in last period as a reference to
% constrain portfolio standards applie to electricity
% generation
% For example: "Renewables must supply 90% of the
% electricity by 2050"
total_load_in_2050 = sum(Load2(7,:));

% Existing storage counts towards the storage standard
existing_storage = sum(plant.Capacity_MW(...
    plant.Category == 'STORAGE'));
Minimum_storage_capacity_in_2050 =...
    peak_load_in_last_period*Storage_standard...
    -existing_storage;

% Existing solar and wind capacity
% (ICAP = Installed CAPacity = existing capacity)
icap_solar =...
    installed_capacity.installed_capacity_MW(...
    installed_capacity.technologies == 'SOLAR');
icap_wind =...
    installed_capacity.installed_capacity_MW(...
    installed_capacity.technologies == 'WIND');

% Solar & wind annual average capacity factor (hr/yr)
CapFactor_solar_annual_average =...
    mean(Var_Energy.Solar) * 8760;    
CapFactor_wind_annual_average =...
    mean(Var_Energy.Wind) * 8760;

% We can only constrain additions and retirements of
% *capacity*, not generation. Since capacity factors
% are fixed, we can move between capacity and
% generation equivalents, but they are different for
% wind and solar, so we take their ratio and get a
% solar-to-wind equivalent capacity factor. This way,
% we can scale up the *capacity* addition/retirement
% constraint to its *generation* equivalent in A.
Solar_to_wind_equivalent_capacity =...
    CapFactor_solar_annual_average...
    /CapFactor_wind_annual_average;

% Existing solar and wind count towards the standard
policy_complied_by_icap =...
    icap_solar*CapFactor_solar_annual_average...
    + icap_wind*CapFactor_wind_annual_average;
Minimum_renewable_generation_in_2050 =...
    total_load_in_2050*Renewable_standard;
Minimum_renewable_capacity_in_2050 =...
    (Minimum_renewable_generation_in_2050...
    -policy_complied_by_icap)...
    /CapFactor_wind_annual_average;

%% Pre-dispatch to set penalties

% get fuel prices
[NG,Coal,Nuc,Oil,Bio] =...
    fuel_prices(2020,s.fuel_data);

% estimate fuel VOM, add to non-fuel VOM to get marginal cost
% plant = table2struct(plant);
[s] = read_data_local(plant,NG,Coal,Nuc,Oil,Bio,zeros(Periods,num_tech),...
            s.name1,s.Cap_cost,num_tech,2020,s.fixed_costs_options,s);

if s.penalty_approach == 0
plant = s.plant;
% run dispatch, get real capacity factors, estimate real LCOEs for
% retirement script (and remove this from read_data)
dispatch_inputs.plant_struct = table2struct(plant);
dispatch_inputs.Load_d = Load2(1,:);
dispatch_inputs.Var_Energy = Var_Energy;
dispatch_inputs.year = 2016;
dispatch_inputs.RPS_target = 0;

[dispatch_outputs] = dispatch_model(dispatch_inputs,s);

x = plant.Category == 'OWN';
plant.MC(x) = dispatch_outputs.highest_LCOE;
s.plant = plant;
end
%% Calibrate storage representative weeks
% k = s.time_cluster; % number of groups/rep. weeks
% 
% df = Power;
% a = 24*7; % equal to a weeks' length in hours
% b = floor(size(df,2)/a); % equal to the number of weeks in the period
% df= df(1:a*b); % get rid of hours past the 52 weeks
% df_matrix= reshape(df,a,b);
% df_weekly_totals= sum(df_matrix)';
% df_weekly_totals(:,2)= 1:b; % add week of the year identifier
% idx= kmeans(df_weekly_totals(:,1),k);
% df_weekly_totals(:,3)= idx;
% 
% % separate data into groups
% repWeekIDs= zeros(k,1);
% % vals= zeros(k,1);
% % repWeeks= zeros(1,columnHeight*k);
% for i=1:k
%     id= df_weekly_totals(:,3) == i;
%     df= df_weekly_totals(id,1:2);
%     avg= repmat(mean(df(:,1)),size(df,1),1);
%     dev= abs(df(:,1)-avg);
%     [~,I]= min(dev);
%     repWeekIDs(i)= df(I,2);
% end
% s.sIn.repWeekIDs= repWeekIDs;
% s.sIn.kmeansMap= df_weekly_totals(:,2:3);

% photw (per-hour-of-the-week) clustering
% StorHours=zeros(a,k);
% for i=1:k
%     df_k= df_matrix(:,s.sIn.kmeansMap((s.sIn.kmeansMap(:,2)==i),1));
%     noStorhours= ~sum(df_k,2);
%     yesStorHours= sum(df_k,2)~=0;
%     StorHours(noStorhours,i)= 0;
%     StorHours(yesStorHours,i)= 1;
% end
% s.sIn.StorHours=StorHours;

%% Export storage case
% Case=0: default: run year with 8760 hours, no compression/clustering
% Case=1: use compression (ignore hours with zero storage output)
% Case=2: run week by week
% Case=3: cluster weeks in groups without compression (pass hours with zero storage output, too)
% Case=4: cluster weeks in groups, and compress them

if cluster
    if compression
        Case=4;    % cluster weeks in groups, and compress them
    else
        Case=3;    % cluster weeks in groups without compression
        % (pass hours with zero storage output, too)
    end
else
    if compression
        Case=1;    % use compression (ignore hours with zero storage output)
        if compression>1
            Case=2;    % run week by week
        end
    else
        Case=0 ;    % default: run year with 8760 hours, no compression/clustering
    end
end

s.sIn.storageCase= Case;

%% Run dispatch and storage if needed
switch Case
    case {1,2,3,4}
        % need to run dispatch to get prices for storage
        plant = s.plant;
        dispatch_inputs.plant_struct = table2struct(plant);
        dispatch_inputs.Load_d = Load2(1,:);
        dispatch_inputs.Var_Energy = Var_Energy;
        dispatch_inputs.year = 2016;
        dispatch_inputs.RPS_target = 0;

        [dispatch_outputs] = dispatch_model(dispatch_inputs,s);

        storage_power = sum(plant.Capacity_MW(...
            plant.Category == 'STORAGE'));

        so_in.ClrPrice = dispatch_outputs.clearing_price;
        so_in.Load2015 = Load2(1,:);
        so_in.cap = storage_power;
        so_in.time_cluster = 0;
        so_in.storageCase= 0;
        [so_out] = storage_optimizer(so_in); % run storage 8760
        Power= so_out.Power; % use 8760 outputs to choose best rep. weeks
end
%% Prepare data for storage cases
mustHours= [1 7 8 9 14 15 16 17 23 24];
inTheFence= [2 3 10 13];
% safelyDelete= [4 5 6 11 12 18 19 20 21 22];
switch Case
    case 1        
        day= zeros(1,24);
        day(1,mustHours)= mustHours;
        day(1,inTheFence)= inTheFence;
        day= logical(day);
        s.sIn.StorHours= repmat(day,1,365);
%         s.sIn.StorHours= logical(Power);
    case 2
        day= zeros(1,24);
        day(1,mustHours)= mustHours;
        day(1,inTheFence)= inTheFence;
        day= logical(day)';
        week= repmat(day,7,1);
        s.sIn.StorHours= repmat(week,1,52);
%         s.sIn.StorHours= reshape(df(1,1:168*52),168,52);
    case 3
        k = cluster; % number of groups/rep. weeks

        df = Power;
        a = 24*7; % equal to a weeks' length in hours
        b = floor(size(df,2)/a); % equal to the number of weeks in the period
        df= df(1:a*b); % get rid of hours past the 52 weeks
        df_matrix= reshape(df,a,b);
        df_weekly_totals= sum(df_matrix)';
        df_weekly_totals(:,2)= 1:b; % add week of the year identifier
        idx= kmeans(df_weekly_totals(:,1),k);
        df_weekly_totals(:,3)= idx;

        % separate data into groups
        repWeekIDs= zeros(k,1);
        for i=1:k
            id= df_weekly_totals(:,3) == i;
            df= df_weekly_totals(id,1:2);
            avg= repmat(mean(df(:,1)),size(df,1),1);
            dev= abs(df(:,1)-avg);
            [~,I]= min(dev);
            repWeekIDs(i)= df(I,2);
        end
        s.sIn.repWeekIDs= repWeekIDs;
        s.sIn.kmeansMap= df_weekly_totals(:,2:3);
    case 4
        k = cluster; % number of groups/rep. weeks
        df = Power;
        a = 24*7; % equal to a weeks' length in hours
        b = floor(size(df,2)/a); % equal to the number of weeks in the period
        df= df(1:a*b); % get rid of hours past the 52 weeks
        df_matrix= reshape(df,a,b);
        df_weekly_totals= sum(df_matrix)';
        df_weekly_totals(:,2)= 1:b; % add week of the year identifier
        idx= kmeans(df_weekly_totals(:,1),k);
        df_weekly_totals(:,3)= idx;

        % separate data into groups
        repWeekIDs= zeros(k,1);
        for i=1:k
            id= df_weekly_totals(:,3) == i;
            df= df_weekly_totals(id,1:2);
            avg= repmat(mean(df(:,1)),size(df,1),1);
            dev= abs(df(:,1)-avg);
            [~,I]= min(dev);
            repWeekIDs(i)= df(I,2);
        end
        s.sIn.repWeekIDs= repWeekIDs;
        s.sIn.kmeansMap= df_weekly_totals(:,2:3);

        % compressing
        day= zeros(1,24);
        day(1,mustHours)= mustHours;
        day(1,inTheFence)= inTheFence;
        day= logical(day)';
        week= repmat(day,7,1);
        s.sIn.StorHours= repmat(week,1,k);
%         s.sIn.StorHours= logical(df_matrix(:,repWeekIDs));

%         s.sIn.StorHours=zeros(a,k);
%         for i=1:k
%             df_k= df_matrix(:,s.sIn.kmeansMap((s.sIn.kmeansMap(:,2)==i),1));            
%             s.sIn.StorHours(:,i)= logical(sum(df_k,2));
%         end
end

%% Suggest feasible starting populations

% While the matlab genetic algorithm documentation
% states that the population search uses a uniform
% distribution, in practice, I have seen the solver
% exploring the bounds rather than values around the
% center of those bounds. So in this section I decided 
% to suggest an initial population that truly followed 
% a uniform distribution.

init_population = ...
    zeros(PopulationSize,Periods*num_tech);

% Use shorter bounds for the initial population
% suggested, otherwise the solver wastes a lot of time
% wandering in the high cost space

% Let model install up to 2 GW of something per year
% (10 GW per 5-year period) by default
Cap_cost.max_cap(:) = 10000;
Cap_cost.max_cap(end) = 0;

if storage == 0
    Cap_cost.max_cap([5 6]) =   150000; 
else
    Cap_cost.max_cap([5 6]) = 150000;
    % For storage, same as renewables
    Cap_cost.max_cap(8) = 150000;
end

max_cap=repmat(Cap_cost.max_cap,1,Periods);
low_cap=repmat(Cap_cost.Retire_cap,1,Periods);
if s.Storage_standard > 0 && s.storage > 0
    low_cap(8,:) = add_storage;
end
lb=reshape(low_cap,1,Periods*num_tech);
ub=reshape(max_cap,1,Periods*num_tech);

% Generate population with individuals within bounds
for tech_period=1:Periods*num_tech
    init_population(:,tech_period) =...
        lb(tech_period) +...
        (ub(tech_period)-lb(tech_period)).*rand(PopulationSize,1);
end

% Plot the initial population suggested, for
% verification purposes during code development
if plot_solvers == 2
    subplot_rows = 1; 
    subplot_columns = 3;
    if Renewable_standard > 0
        subplot_columns = subplot_columns + 1;
    end
    if Storage_standard > 0
        subplot_columns = subplot_columns + 1;
    end
    
    hAxes(1) = subplot(subplot_rows,subplot_columns,1);
    check_totals = zeros(PopulationSize,num_tech);
    for indiv=1:PopulationSize
        for tech_period=1:num_tech
            check_totals(indiv,tech_period) =...
            sum(init_population(...
            indiv,tech_period:num_tech:num_tech*Periods));
        end
    end
    plot(check_totals(:,1),'LineWidth',2)
    hold on
    for tech_period=2:num_tech
        plot(check_totals(:,tech_period),'LineWidth',2)
    end
    set(gca,'FontSize',12)
    set(gcf,'Position',...
        [100 100 600*subplot_columns 500*subplot_rows])
    hold off
    legend(Cap_cost.Tech)

end

% Generate population of individuals that do not 
% retire more capacity than the total existing capacity
% for each technology category
iterations = 0;
iterations_per_tech = zeros(1,num_tech);
for indiv=1:PopulationSize
for tech_period=1:(num_tech-1)
% for period = 1:Periods
is_individual_feasible = 1;
while is_individual_feasible > 0

    init_population(indiv,...
        tech_period:num_tech:num_tech*Periods) =...
        lb(tech_period) +...
        (ub(tech_period)-lb(tech_period)).*rand(1,Periods);

    % cumulative additions & retirements
    xx = init_population(indiv,...
        tech_period:num_tech:num_tech*Periods);
    xx = cumsum(xx')';

    % installed capacity
    icap = -repmat(Cap_cost.Retire_cap(tech_period)...
        ,1,Periods);

    new_icap = icap + xx;

    % if we reach gas_CC retirements, sum R&As
    % to check for this constraint
    if tech_period == num_tech
        yy = init_population(indiv,...
            2:num_tech:num_tech*Periods);
        xx = xx + yy;
        xx = cumsum(xx')';

        % installed capacity
        icap = -repmat(Cap_cost.Retire_cap(tech_period)...
            ,1,Periods);

        new_icap = icap + xx;
    end

    % if the individual retires more than
    % exists in any period, new_icap will have
    % negative values, so we check that
    % new_icap is always above zero; if it is
    % not, then the individual isn't feasible
    is_individual_feasible =...
        new_icap <= 0;
    is_individual_feasible = double(...
        is_individual_feasible);
    is_individual_feasible =...
        sum(is_individual_feasible);

    %             is_individual_feasible =...
    %                 sum(init_population(indiv,...
    %                 tech:num_tech:num_tech*Periods))...
    %                 > Cap_cost.Retire_cap(tech);
    iterations = iterations+1;
end
iterations_per_tech(tech_period) = iterations;
iterations = 0;
% end
end
% disp(iterations_per_tech)
end

%% Suggest initial retirements of gas CC
% Generate population of individuals that do not 
% retire more capacity than the total existing capacity
% for gas CC
gas_CC_add = 2;
gas_CC_remove = 8 + s.storage ;
icap = -lb(gas_CC_remove); % existing capacity
for indiv=1:PopulationSize
    % For the first period, we can only retire the
    % existing capacity, which was already implemented
    % in the previous for loops.
    % For the subsequent periods, we can only retire
    % whatever is left after adding additions and
    % subtracting retirements from all periods before
    % the current period
for period = 2:Periods

    cumulative_additions =...
        sum(init_population(indiv,...
        gas_CC_add:num_tech:num_tech*(period-1)));

    cumulative_retirements =...
        sum(init_population(indiv,...
        gas_CC_remove:num_tech:num_tech*(period-1)));
    
    tech_period = gas_CC_remove*period;
    new_gas_CC_lb = icap +...
        cumulative_additions +...
        cumulative_retirements;
    lb(tech_period) = -new_gas_CC_lb;

    init_population(indiv,tech_period) =...
        lb(tech_period) +...
        (ub(tech_period)-lb(tech_period)).*rand(1);
end
end

% Plot the initial population suggested, for
% verification purposes during code development
if plot_solvers == 2
    hAxes(2) = subplot(subplot_rows,subplot_columns,2);
    check_totals = zeros(PopulationSize,num_tech);
    for indiv=1:PopulationSize
        for tech_period=1:num_tech
            check_totals(indiv,tech_period) =...
                sum(init_population(...
                indiv,tech_period:num_tech:num_tech*Periods));
        end
    end
    plot(check_totals(:,1),'LineWidth',2)
    hold on
    for tech_period=2:num_tech
        plot(check_totals(:,tech_period),'LineWidth',2)
    end
    set(gca,'FontSize',12)
    set(gcf,'Position',...
        [100 100 600*subplot_columns 500*subplot_rows])
    hold off
    linkaxes(hAxes, 'y' )
end

lb_mat = init_population;


%% Generate individuals that ensure ICAP > Load

% derating multipliers
CF = ones(size(Cap_cost.Tech))';
CF(5) = 0.224; CF(6) = 0.348;

iterations = 0;
iterations_per_period = zeros(1,Periods);
for indiv=1:PopulationSize
% for period=1:(num_tech-1)
for period = 1:Periods
is_unfeasible = 1;
while is_unfeasible > 0
    
    lb = lb_mat(indiv,:);

    index =...
        num_tech*period-(num_tech-1):num_tech*period;
    
    init_population(indiv,index)=...
        lb(index) +...
        (ub(index)-lb(index)).*rand(1,num_tech);

    delta_icap = init_population(indiv,index);

    icap = -Cap_cost.Retire_cap';
    if period > 1
        for i = period-1:-1:1
            j = num_tech*i-(num_tech-1):num_tech*i;
            icap = icap...
                + init_population(indiv,j);
        end
    end
    new_icap = icap + delta_icap;
    new_icap = new_icap.*CF;
    all_icap = sum(new_icap);
    is_unfeasible=...
        all_icap <= Load2(period,s.I(period));
    is_unfeasible = double(...
        is_unfeasible);
    
    iterations = iterations+1;
end
iterations_per_period(period) = iterations;
iterations = 0;
% end
end
% disp(iterations_per_period);
end

% Plot the initial population suggested, for
% verification purposes during code development
if plot_solvers == 2
    hAxes(3) = subplot(subplot_rows,subplot_columns,3);
    check_totals = zeros(PopulationSize,num_tech);
    for indiv=1:PopulationSize
        for tech_period=1:num_tech
            check_totals(indiv,tech_period) =...
                sum(init_population(...
                indiv,tech_period:num_tech:num_tech*Periods));
        end
    end
    plot(check_totals(:,1),'LineWidth',2)
    hold on
    for tech_period=2:num_tech
        plot(check_totals(:,tech_period),'LineWidth',2)
    end
    set(gca,'FontSize',12)
    set(gcf,'Position',...
        [100 100 600*subplot_columns 500*subplot_rows])
    hold off
    linkaxes(hAxes, 'y' )
end

%% Adjust initial population to comply with RPS
% Generate population with individuals that comply with
% renewable portfolio standards
iterations = 0;
if Renewable_standard > 0
    for indiv=1:PopulationSize
        tech_period = [5 6];
        columns = [tech_period(1):num_tech:num_tech*Periods...
            tech_period(2):num_tech:num_tech*Periods];
        is_individual_feasible = 0;
        while is_individual_feasible == 0
            init_population(indiv,columns) =...
                zeros(1,2*Periods)...
                +(repelem(100000,2*Periods)-...
                zeros(1,2*Periods)).*rand(1,2*Periods);
            is_individual_feasible =...
                sum(init_population(indiv,columns))...
                > Minimum_renewable_capacity_in_2050;
            iterations = iterations+1;
        end
        iterations_per_tech(tech_period) =...
            iterations_per_tech(tech_period)+iterations;
        iterations = 0;
    end
end

if plot_solvers == 2 && Renewable_standard > 0
    subplot(subplot_rows,subplot_columns,3)
    check_totals = zeros(PopulationSize,num_tech);
    for indiv=1:PopulationSize
        for tech_period=1:num_tech
            check_totals(indiv,tech_period) = ...
                sum(init_population(...
                indiv,tech_period:num_tech:num_tech*Periods));
        end
    end
    plot(check_totals(:,1),'LineWidth',2)
    hold on
    for tech_period=2:num_tech
        plot(check_totals(:,tech_period),'LineWidth',2)
    end
    hold off
end

%% Adjust initial population to comply with SPS
% Generate population with individuals that comply with
% storage portfolio standards
iterations = 0;
if Storage_standard > 0
    for indiv=1:PopulationSize
        tech_period = 8;
        columns = tech_period(1):num_tech:num_tech*Periods;
        is_individual_feasible = 0;
        while is_individual_feasible == 0
            init_population(indiv,columns) =...
                repmat(lb(tech_period),1,Periods)...
                +(repelem(50000,Periods)-repmat(...
                lb(tech_period),1,Periods)).*rand(1,Periods);
            is_individual_feasible =...
                sum(init_population(indiv,columns))...
                > Minimum_storage_capacity_in_2050;
            iterations = iterations+1;
        end
        iterations_per_tech(tech_period) =...
            iterations_per_tech(tech_period)+iterations;
        iterations = 0;
    end
end

% Generate population with individuals that comply with
% renewable portfolio standards
if plot_solvers == 2 && Storage_standard > 0
    subplot(subplot_rows,subplot_columns,...
        2+double(Renewable_standard>0)+1)
    check_totals = zeros(PopulationSize,num_tech);
    for indiv=1:PopulationSize
        for tech_period=1:num_tech
            check_totals(indiv,tech_period) =...
                sum(init_population(...
                indiv,tech_period:num_tech:num_tech*Periods));
        end
    end
    plot(check_totals(:,1),'LineWidth',2)
    hold on
    for tech_period=2:num_tech
        plot(check_totals(:,tech_period),'LineWidth',2)
    end
    hold off
end

%% start from optimal
try
    if s.runFromOptimal
        theFiles = dir('**/*.mat');
        for k= 1:length(theFiles)
            if startsWith(theFiles(k).name,s.runFromOptimal)
                baseFileName= theFiles(k).name;
                fullFileName=...
                    fullfile(theFiles(k).folder,baseFileName);
                baseOutputs= load(fullFileName,'outputs');
                baseS= load(fullFileName,'s');
                baseOptimal= baseOutputs.outputs.Optimal_DeltaMW;
                reshaped= reshape(baseOptimal,[],7);
                reshaped2= reshaped; % make it a matrix
                reshaped2(9,:)=reshaped(8,:); % adapt to gas_ret format
                reshaped2(8,:)=zeros(1,size(reshaped2,2)); % add row for storage
                reshaped3= reshape(reshaped2,1,[]); % make it vector again

                % change initial population so they are
                % all small variations of the optimal
                % and storage is not always zero
                init_population= repmat(reshaped3,size(init_population,1),1);
                lb2= zeros(1,PopulationSize);
                ub2= repelem(20000,PopulationSize);
                for period=1:Periods
                    steps= period*num_tech-1;
                    init_population(:,steps) =...
                        lb2' +...
                        (ub2-lb2)'.*rand(PopulationSize,1);
                end
                % but keep the original population
                init_population(1,:)= reshaped3;
%                 lb2= reshaped3-reshaped3*.001;
%                 ub2= reshaped3+reshaped3*.0010;
%                 for indiv=2:PopulationSize
%                     init_population(indiv,:)=...
%                         lb2 +...
%                         (ub2-lb2).*rand(length(reshaped3),1)';
%                 end

            end
        end
    end
catch
end

%% stage 1 params
try
    use_two_stage = s.use_two_stage;
catch
    s.use_two_stage = 0;
end

if s.use_two_stage == 1
s.stage1.nvars = num_tech;
s.stage1.ub = Cap_cost.max_cap;
s.stage1.lb = Cap_cost.Retire_cap;
s.stage1.ga_options = optimoptions( ...
    'ga'...
    ,'UseParallel',true...
    ,'FunctionTolerance',FunctionTolerance_GA...
    ,'MaxGenerations',MaxGenerations...
    ...,'HybridFcn',{@patternsearch,ps_options}...
    ,'Display','iter'...
    ,'MutationFcn',{...
        @mutationgaussian...
        mutation_scale mutation_shrink}...
    ...,'InitialPopulationMatrix',init_population...
    ,'PopulationSize',PopulationSize...
    ...,'OutputFcn',@checkpoint...
    );
end

%% Find least-cost retirements and additions

% Use solver to find least-cost capacity additions and
% retirements. We use a hybrid of genetic algorithm and
% pattern search matlab solvers.
ps_options = optimoptions( ...
    'patternsearch'...
    ,'UseParallel',true...
    ,'MeshTolerance',Mesh_Tolerance...
    ,'FunctionTolerance',FunctionTolerance_PS...
    ,'Display','iter'...
    ...,'OutputFcn',@checkpoint...
    );
ga_options = optimoptions( ...
    'ga'...
    ,'UseParallel',true...
    ,'FunctionTolerance',FunctionTolerance_GA...
    ,'MaxGenerations',MaxGenerations...
    ...,'HybridFcn',{@patternsearch,ps_options}...
    ,'Display','iter'...
    ,'MutationFcn',{...
        @mutationgaussian...
        mutation_scale mutation_shrink}...
    ,'InitialPopulationMatrix',init_population...
    ,'PopulationSize',PopulationSize...
    ...,'OutputFcn',@checkpoint...
    );

% If you want to see plots with the progress of the
% algorithm, plot_solvers should be equal to 1. This
% should not be used if running in the cluster or
% remotely without a graphic user interface (GUI).
if plot_solvers == 1
    ps_options = optimoptions( ...
        ps_options,...
        'PlotFcn',...
            {@psplotbestf, ...
            @psplotfuncount,...
            @psplotmeshsize, ...
            @psplotbestx} ...
        );
    ga_options = optimoptions( ...
        ga_options,...
        'PlotFcn', ...
            {@gaplotbestf, ...
            @gaplotbestindiv} ...
        );
end

% if running a test, just run a few polls
if s.devDebug == 1
    ps_options = optimoptions( ...
        ps_options,...
        'MaxFunctionEvaluations',800, ...
        'MaxTime',60*5, ...
        'MaxIterations',10);
end

if s.devDebug == 1
    ps_options = optimoptions(ps_options,'UseParallel',false);
    ga_options = optimoptions(ga_options,'UseParallel',false);
end

% make sure new ps_options are chained to ga_options
ga_options = optimoptions(ga_options, ...
    'HybridFcn',{@patternsearch,ps_options}...
    );
s.ga_options = ga_options;
s.pop = init_population;
%% Rest of the code for solver

% ObjectiveFunction = @(x) parameterized_objective(s,x);

% Run the solver
% [OPT_x,system_cost,~,output] = ga(ObjectiveFunction,...
%     nvars,A,b,[],[],lb,ub,[],ga_options);

% Store solver outputs
% Optimal_DeltaMW = OPT_x;
% outputs.Optimal_DeltaMW = Optimal_DeltaMW;
% outputs.system_cost = system_cost;
% outputs.generations = output.generations;

%% Export subsidy cost (government spending)
% Popbuild = zeros(num_tech,BPeriods*5);
% OPT_x = reshape(OPT_x,num_tech,Periods);
% peryear = OPT_x/5;
% Popbuild(:,1:Periods*5) = repelem(peryear,1,5);
% subsidy_cost = repmat(subsidies,1,BPeriods*5);
% subsidy_cost(Popbuild<0) = 0;
% subsidy_cost = subsidy_cost.*Popbuild*1000;
% subsidy_cost = sum(subsidy_cost);
% real = (1/(1+interest)).^(0:BPeriods*5-1);
% subsidy_cost =...
%     sum(subsidy_cost.*repmat(real,1,monte))/10^9;
% outputs.subsidy_cost = subsidy_cost;

%% Export runtimes
% runtime = toc;
% runtime = seconds(runtime);
% runtime.Format = 'hh:mm';
% outputs.runtime = runtime;

end

%% Local read_data function to get highest_LCOE
function [s] =...
        read_data_local(plant,NG, Coal, Nuc,Oil,Bio,...
        cap1,name1,Cap_cost2,num_tech,year,...
        fixed_costs_options,s) %#ok<INUSL> 
%% start
warning off
% plant=struct2table(plant);
EPA_FOM = s.EPA_FOM;

% convert to categorical data types
Cap_cost2.Tech = categorical(Cap_cost2.Tech);
Cap_cost2.Fuel = categorical(Cap_cost2.Fuel);
Cap_cost2.Prime_Mover =...
    categorical(Cap_cost2.Prime_Mover);

% % add emissions of new power plants
% Cap_cost2.CO2eq_lbMWh(Cap_cost2.Fuel == 'COAL') =...
%     0.09552 * 1/1000000 *...
%     Cap_cost2.HR_btuKWh(Cap_cost2.Fuel == 'COAL') *...
%     1000 * 2204.62;
% Cap_cost2.CO2eq_lbMWh(Cap_cost2.Fuel == 'GAS') =...
%     0.05306 * 1/1000000 *...
%     Cap_cost2.HR_btuKWh(Cap_cost2.Fuel == 'GAS')...
%     * 1000 * 2204.62;
% Cap_cost2.CO2eq_lbMWh(Cap_cost2.Fuel == 'OIL') =...
%     0.07315 * 1/1000000 *...
%     Cap_cost2.HR_btuKWh(Cap_cost2.Fuel == 'OIL')...
%     * 1000 * 2204.62;
% 
% %% add new plants
% for cc = 1:num_tech
%     plant.Plant_Name(end+1)= strcat('new_',string(...
%         Cap_cost2.Tech(cc)),'_',string(year)); ...name1(cc);
%     plant.Category(end)=Cap_cost2.Tech(cc);
%     plant.Heatrate_BtukWh(end)=Cap_cost2.HR_btuKWh(cc);
%     plant.Prime_Mover(end)=Cap_cost2.Prime_Mover(cc);
%     plant.UB(end)=Cap_cost2.UB(cc);
%     plant.LB(end)=Cap_cost2.LB(cc);
%     plant.OnM_dolkW(end)=Cap_cost2.OnM_dolkW(cc);
%     plant.Var_dolMWh(end)=Cap_cost2.Var_dolMWh(cc);
%     plant.CO2eq_lbMWh(end)=Cap_cost2.CO2eq_lbMWh(cc);
%     plant.OnLineYear(end)=year;
% end
% 
% plant.Capacity_MW(end-num_tech+1:end)=...
%     cap1(1:num_tech,1);

%% Fixed costs

% use Sargent & Lundy FOM
if fixed_costs_options == 1

%% fixed costs for single-value technologies
% from Sargent & Lundy (2017$/kW-Yr)
gas_oil_steam_turbine_plants = plant.PlantType ==...
    "O/G Steam" & plant.Capacity_MW <= 500;
plant.OnM_dolkW(gas_oil_steam_turbine_plants) = 50.43;
gas_oil_steam_turbine_plants = plant.PlantType ==...
    "O/G Steam" & plant.Capacity_MW > 500 &...
    plant.Capacity_MW <= 1000;
plant.OnM_dolkW(gas_oil_steam_turbine_plants) = 31.39;
gas_oil_steam_turbine_plants = plant.PlantType ==...
    "O/G Steam" & plant.Capacity_MW > 1000;
plant.OnM_dolkW(gas_oil_steam_turbine_plants) = 27.17;

gas_oil_combined_cycle = plant.PlantType ==...
    "Combined Cycle" & plant.Capacity_MW <= 500;
plant.OnM_dolkW(gas_oil_combined_cycle) = 33;
gas_oil_combined_cycle = plant.PlantType ==...
    "Combined Cycle" & plant.Capacity_MW > 500 &...
    plant.Capacity_MW <= 1000;
plant.OnM_dolkW(gas_oil_combined_cycle) = 23.05;
gas_oil_combined_cycle = plant.PlantType ==...
    "Combined Cycle" & plant.Capacity_MW > 1000;
plant.OnM_dolkW(gas_oil_combined_cycle) = 25.25;

gas_oil_combustion_turbines = plant.PlantType ==...
    "Combustion Turbine" & plant.Capacity_MW <= 100;
plant.OnM_dolkW(gas_oil_combustion_turbines) = 14.96;
gas_oil_combustion_turbines = plant.PlantType ==...
    "Combustion Turbine" & plant.Capacity_MW > 100 &...
    plant.Capacity_MW <= 300;
plant.OnM_dolkW(gas_oil_combustion_turbines) = 12.61;
gas_oil_combustion_turbines = plant.PlantType ==...
    "Combustion Turbine" & plant.Capacity_MW > 300;
plant.OnM_dolkW(gas_oil_combustion_turbines) = 10.94;

solar_plants = plant.PlantType == "Solar PV" &...
    plant.Capacity_MW <= 5;
plant.OnM_dolkW(solar_plants) = 67;
solar_plants = plant.PlantType == "Solar PV" &...
    plant.Capacity_MW > 5;
plant.OnM_dolkW(solar_plants) = 36;

pumped_storage = plant.PlantType == "Pumped Storage";
plant.OnM_dolkW(pumped_storage) = 41.46;

% from EPA's IPM Documentation Chapter 4 Table 4-9 
% FOM Assumptions in v6 (2019$/kW-Yr)
biomass_plants = plant.PlantType == "Biomass";
plant.OnM_dolkW(biomass_plants) = 149.3;

landfill_gas_MSW_and_fossil_waste_plants =...
    plant.PlantType == "Fossil Waste" |...
    plant.PlantType == "Landfill Gas" |...
    plant.PlantType == "Municipal Solid Waste" |...
    plant.PlantType == "Non-Fossil Waste";
plant.OnM_dolkW(...
    landfill_gas_MSW_and_fossil_waste_plants) = 259.23;

igcc_plants = plant.PlantType == "IGCC";
plant.OnM_dolkW(igcc_plants) = 108.71;

%% fixed costs for age-dependent non-coal plants
plant.Age = year - plant.OnLineYear;
plant.Age(plant.Age<0) = 0;

conventional_hydroelectric_plants =...
    plant.PlantType == "Hydro";
plant.OnM_dolkW(conventional_hydroelectric_plants) =...
    29.629 + (0.369*plant.Age(...
    conventional_hydroelectric_plants));

wind_plants = plant.PlantType == "Onshore Wind" |...
    plant.PlantType == "Offshore Wind";
wind_plants_2 = wind_plants & plant.Capacity_MW <= 100;
plant.OnM_dolkW(wind_plants_2) = 59.56 +...
    (1.22*plant.Age(wind_plants_2));
wind_plants_2 = wind_plants &...
    plant.Capacity_MW > 100 & plant.Capacity_MW <= 200;
plant.OnM_dolkW(wind_plants_2) = 39.83 +...
    (1.17*plant.Age(wind_plants_2));
wind_plants_2 = wind_plants & plant.Capacity_MW > 200;
plant.OnM_dolkW(wind_plants_2) = 40.26 +...
    (0.92*plant.Age(wind_plants_2));

%% fixed costs for age-dependent coal plants
plant.FGD = repelem(0,length(plant.SO2Control))'; 
plant.FGD(plant.SO2Control ~= "No_SO2_control") =...
    repelem(1,length(plant.FGD(plant.SO2Control ~=...
    "No_SO2_control")))';
coal_plants = plant.PlantType ==...
    "Coal Steam" & plant.Capacity_MW <= 500;
plant.OnM_dolkW(coal_plants) = 69.94 +...
(0.126*plant.Age(coal_plants)) +...
(5.68*plant.FGD(coal_plants));
coal_plants = plant.PlantType ==...
    "Coal Steam" &...
    plant.Capacity_MW > 500 &...
    plant.Capacity_MW <= 1000;
plant.OnM_dolkW(coal_plants) = 59.75 +...
    (0.126*plant.Age(coal_plants)) +...
    (5.68*plant.FGD(coal_plants));
coal_plants = plant.PlantType ==...
    "Coal Steam" &...
    plant.Capacity_MW > 1000 &...
    plant.Capacity_MW <= 2000;
plant.OnM_dolkW(coal_plants) = 54.25 +...
    (0.126*plant.Age(coal_plants)) +...
    (5.68*plant.FGD(coal_plants));
coal_plants = plant.PlantType ==...
    "Coal Steam" & plant.Capacity_MW > 2000;
plant.OnM_dolkW(coal_plants) = 59 +...
    (0.126*plant.Age(coal_plants)) +...
    (5.68*plant.FGD(coal_plants));
plant.Age = []; plant.FGD = [];
end
% use EIA FOM (like above but only coal varies with age)
if fixed_costs_options == 2

%% fixed costs for single-value technologies
% from Sargent & Lundy (2017$/kW-Yr)
gas_oil_steam_turbine_plants = plant.PlantType ==...
    "O/G Steam" & plant.Capacity_MW <= 500;
plant.OnM_dolkW(gas_oil_steam_turbine_plants) = 50.43;
gas_oil_steam_turbine_plants = plant.PlantType ==...
    "O/G Steam" & plant.Capacity_MW > 500 &...
    plant.Capacity_MW <= 1000;
plant.OnM_dolkW(gas_oil_steam_turbine_plants) = 31.39;
gas_oil_steam_turbine_plants = plant.PlantType ==...
    "O/G Steam" & plant.Capacity_MW > 1000;
plant.OnM_dolkW(gas_oil_steam_turbine_plants) = 27.17;

gas_oil_combined_cycle = plant.PlantType ==...
    "Combined Cycle" & plant.Capacity_MW <= 500;
plant.OnM_dolkW(gas_oil_combined_cycle) = 33;
gas_oil_combined_cycle = plant.PlantType ==...
    "Combined Cycle" & plant.Capacity_MW > 500 &...
    plant.Capacity_MW <= 1000;
plant.OnM_dolkW(gas_oil_combined_cycle) = 23.05;
gas_oil_combined_cycle = plant.PlantType ==...
    "Combined Cycle" & plant.Capacity_MW > 1000;
plant.OnM_dolkW(gas_oil_combined_cycle) = 25.25;

gas_oil_combustion_turbines = plant.PlantType ==...
    "Combustion Turbine" & plant.Capacity_MW <= 100;
plant.OnM_dolkW(gas_oil_combustion_turbines) = 14.96;
gas_oil_combustion_turbines = plant.PlantType ==...
    "Combustion Turbine" & plant.Capacity_MW > 100 &...
    plant.Capacity_MW <= 300;
plant.OnM_dolkW(gas_oil_combustion_turbines) = 12.61;
gas_oil_combustion_turbines = plant.PlantType ==...
    "Combustion Turbine" & plant.Capacity_MW > 300;
plant.OnM_dolkW(gas_oil_combustion_turbines) = 10.94;

solar_plants = plant.PlantType == "Solar PV";
plant.OnM_dolkW(solar_plants) = 18;

pumped_storage = plant.PlantType == "Pumped Storage";
plant.OnM_dolkW(pumped_storage) = 38.46;

% from EPA's IPM Documentation Chapter 4 Table 4-9 
% FOM Assumptions in v6 (2019$/kW-Yr)
biomass_plants = plant.PlantType == "Biomass";
plant.OnM_dolkW(biomass_plants) = 149.3;

landfill_gas_MSW_and_fossil_waste_plants =...
    plant.PlantType == "Fossil Waste" |...
    plant.PlantType == "Landfill Gas" |...
    plant.PlantType == "Municipal Solid Waste" |...
    plant.PlantType == "Non-Fossil Waste";
plant.OnM_dolkW(...
    landfill_gas_MSW_and_fossil_waste_plants) = 259.23;

igcc_plants = plant.PlantType == "IGCC";
plant.OnM_dolkW(igcc_plants) = 108.71;

%% fixed costs for age-dependent non-coal plants
plant.Age = year - plant.OnLineYear;
plant.Age(plant.Age<0) = 0;

conventional_hydroelectric_plants =...
    plant.PlantType == "Hydro";
plant.OnM_dolkW(conventional_hydroelectric_plants) =...
    44.5;

wind_plants = plant.PlantType == "Onshore Wind" |...
    plant.PlantType == "Offshore Wind";
plant.OnM_dolkW(wind_plants) = 39;

%% fixed costs for age-dependent coal plants
plant.FGD = repelem(0,length(plant.SO2Control))'; 
plant.FGD(plant.SO2Control ~= "No_SO2_control") =...
    repelem(1,length(plant.FGD(plant.SO2Control ~=...
    "No_SO2_control")))';
coal_plants = plant.PlantType ==...
    "Coal Steam" & plant.Capacity_MW <= 500;
plant.OnM_dolkW(coal_plants) = 69.94 +...
(0.126*plant.Age(coal_plants)) +...
(5.68*plant.FGD(coal_plants));
coal_plants = plant.PlantType ==...
    "Coal Steam" &...
    plant.Capacity_MW > 500 &...
    plant.Capacity_MW <= 1000;
plant.OnM_dolkW(coal_plants) = 59.75 +...
    (0.126*plant.Age(coal_plants)) +...
    (5.68*plant.FGD(coal_plants));
coal_plants = plant.PlantType ==...
    "Coal Steam" &...
    plant.Capacity_MW > 1000 &...
    plant.Capacity_MW <= 2000;
plant.OnM_dolkW(coal_plants) = 54.25 +...
    (0.126*plant.Age(coal_plants)) +...
    (5.68*plant.FGD(coal_plants));
coal_plants = plant.PlantType ==...
    "Coal Steam" & plant.Capacity_MW > 2000;
plant.OnM_dolkW(coal_plants) = 59 +...
    (0.126*plant.Age(coal_plants)) +...
    (5.68*plant.FGD(coal_plants));
plant.Age = []; plant.FGD = [];
end
% use EPA FOM
if fixed_costs_options == 3
plant.Age = year - plant.OnLineYear;
plant.Age(plant.Age<0) = 0;

plant.AgeCategory =...
    repmat("All_Years",size(plant.Age));
plant.AgeCategory = categorical(plant.AgeCategory);
aging_plants = plant.PlantType == "Coal Steam" |...
    plant.PlantType == "O/G Steam";
plant.AgeCategory(aging_plants & plant.Age < 40) =...
    "0_to_40_Years";
plant.AgeCategory(aging_plants & plant.Age > 50) =...
    "Greater_than_50_Years";
plant.AgeCategory(aging_plants & plant.Age <= 50 ...
    & plant.Age >= 40) =...
    "40_to_50_Years";

plant.Emissions_control_all =...
    join([...
    string(plant.Emissions_control_all) ...
    string(plant.AgeCategory) ...
    ]);

% determine rows of NEEDS that match with EPA_v6
[idxA,idxB] =...
    ismember(...
    plant.Emissions_control_all,...
    EPA_FOM.Emissions_control_all); 
plant(idxA,'OnM_dolkW') =...
    EPA_FOM(idxB(idxA),'FOM_2019__kW_Yr_'); 

plant.Age = []; plant.AgeCategory = [];
end
% override fixed costs for storage
storage = plant.Category == "STORAGE";... &...
%     plant.OnLineYear < 2020;
plant.OnM_dolkW(storage) = 24.93;
% used value suggested in Table 10 for small PHS plants
% https://www.pnnl.gov/sites/default/files/media/file/PSH_Methodology_0.pdf
% 2023_03_03 update: FOM varies with plant size but is around 24 which
% matches the FOM of new storage plants in AEO, so I'll just use that

%% estimate marginal costs for new plants
% add carbon tax
try
    if s.carbonTax>0 
        plant.carbonTax = s.carbonTax * plant.CO2eq_lbMWh;
    end
catch
end

x = plant.Category == 'GAS_CC';
plant.MC(x) = NG * plant.Heatrate_BtukWh(x)* 10^-3;

x = plant.Category == 'GAS_CT';
plant.MC(x) = NG * plant.Heatrate_BtukWh(x)* 10^-3;

x = plant.Category == 'COAL';
plant.MC(x) = Coal * plant.Heatrate_BtukWh(x)* 10^-3;

x = plant.Category == 'NUCLEAR';
plant.MC(x) = Nuc * plant.Heatrate_BtukWh(x)* 10^-3;

x = plant.Category == 'OIL';
plant.MC(x) = Oil * plant.Heatrate_BtukWh(x)* 10^-3;

x = plant.Category == 'BIOMASS';
plant.MC(x) = Bio * plant.Heatrate_BtukWh(x)* 10^-3;

% plant.MC here in reality means fuel costs, not
% marginal costs, since marginal costs consider the
% variables costs too; this summation is wrongly 
% represented here with plant.TC; I left it as is for
% now but it should be changed when time permits
try
    plant.TC = plant.MC + plant.Var_dolMWh + plant.carbonTax;
catch
    plant.TC = plant.MC + plant.Var_dolMWh;
end

% systematically set penalty TC (marginal cost) to be higher than the
% costliest marginal cost in the existing system by a factor given by the
% penalty slack
x = plant.Category == 'OWN';
plant.OnLineYear(x) = 2016;
costliest_VOM = max(plant.TC(~x));
costliest_FOM = max(plant.OnM_dolkW(~x));
if s.penalty_approach == 1 % use FOM for VOM
    plant.MC(x) = costliest_FOM * 1000;
end
if s.penalty_approach == 2 % let FOM dilute by penalty annual MWh
    plant.MC(x) = costliest_VOM * s.penalty_slack;
end
if s.penalty_approach == 3
    plant.MC(x) = costliest_VOM * s.penalty_slack;
end
s.penalty_FOM = costliest_FOM * s.penalty_slack;
% plant.TC(x)=s.unmet_demand_penalization;
plant.TC(x) = plant.MC(x) + plant.Var_dolMWh(x);
%% roughly estimate LCOE for retirement criteria
% these capacity factors were estimated as the weighted
% average capacity factor for each technology category
% from running the dispatch model in the first period
% GAS_CC_CF   = 0.714571;
% GAS_CT_CF   = 0.015963;
% COAL_CF     = 0.38116;
% NUCLEAR_CF  = 1;
% OIL_CF      = 0.005; % <- assumed 0.5%, but was zero
% BIOMASS_CF  = 0.5916;
% 
% x = plant.Category == 'GAS_CC';
% plant.LCOE_dolMWh(x) =...
%     ((plant.OnM_dolkW(x)...
%     .* plant.Capacity_MW(x) ...
%     * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
%     * 8760 * GAS_CC_CF)...
%     ./ (plant.Capacity_MW(x) * 8760 * GAS_CC_CF);
% 
% x = plant.Category == 'GAS_CT';
% plant.LCOE_dolMWh(x) =...
%     ((plant.OnM_dolkW(x)...
%     .* plant.Capacity_MW(x) ...
%     * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
%     * 8760 * GAS_CT_CF)...
%     ./ (plant.Capacity_MW(x) * 8760 * GAS_CT_CF);
% 
% x = plant.Category == 'COAL';
% plant.LCOE_dolMWh(x) =...
%     ((plant.OnM_dolkW(x)...
%     .* plant.Capacity_MW(x) ...
%     * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
%     * 8760 * COAL_CF)...
%     ./ (plant.Capacity_MW(x) * 8760 * COAL_CF);
% 
% x = plant.Category == 'NUCLEAR';
% plant.LCOE_dolMWh(x) =...
%     ((plant.OnM_dolkW(x)...
%     .* plant.Capacity_MW(x) ...
%     * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
%     * 8760 * NUCLEAR_CF)...
%     ./ (plant.Capacity_MW(x) * 8760 * NUCLEAR_CF);
% 
% x = plant.Category == 'OIL';
% plant.LCOE_dolMWh(x) =...
%     ((plant.OnM_dolkW(x)...
%     .* plant.Capacity_MW(x) ...
%     * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
%     * 8760 * OIL_CF)...
%     ./ (plant.Capacity_MW(x) * 8760 * OIL_CF);
% 
% x = plant.Category == 'BIOMASS';
% plant.LCOE_dolMWh(x) =...
%     ((plant.OnM_dolkW(x)...
%     .* plant.Capacity_MW(x) ...
%     * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
%     * 8760 * BIOMASS_CF)...
%     ./ (plant.Capacity_MW(x) * 8760 * BIOMASS_CF);

plant = sortrows(plant,'TC');
s.plant = plant;
%         x = isnan(plant.TC);
%         plant(x,:)=[];
% plant = table2struct(plant,'ToScalar',true);
end

function [Load2,Var_Energy] = compress_time(s,Load2,Var_Energy) %#ok<INUSL> 
    i=1;

    spring = 78*24+1:170*24;            % 92 days
    summer = 170*24+1:264*24;           % 94 days
    fall = 264*24+1:354*24;             % 90 days
    winter = [1:78*24,354*24+1:8760];   % 89 days

    %               R           G       B           alpha
    spring_color =  [1,         0,      1,          0.2];
    summer_color =  [0.929,     0.694,  0.125,      0.2];
    fall_color =    [0.,        1,      0,          0.2];
    winter_color =  [0,         0,      1,          0.2];
    linewidth = 1;

    MEF1_marker = '-o';

    var = Var_Energy.Solar(:);
    var_spring = reshape(var(spring),24,[]);

%     plot(1:24,MEF1_spring,'Color',spring_color,'LineWidth',linewidth);
%     hold on
%     plot(1:24,MEF1_summer,'Color',summer_color,'LineWidth',linewidth);
%     plot(1:24,MEF1_fall,'Color',fall_color,'LineWidth',linewidth);
%     plot(1:24,MEF1_winter,'Color',winter_color,'LineWidth',linewidth);
%     hold off
%     box off; ylabel('tCO_2e/MWh')
%     set(gca,'FontSize',15)
%     %title('Incremental MEF 24SA (tCO2e/MWh)')
%     legend({'spring','summer','fall','winter'},'Location','best'); legend boxoff
%     saveas(gcf,strcat(Plotpath,region,'_','MEF1_24SA.png'));

    plot(1:24,var_spring(:,1),MEF1_marker,'Color',spring_color);
    hold on
    for day = size(var_spring,2):-1:2
        plot(1:24,var_spring(:,day),MEF1_marker,'Color',spring_color);
    end
    hold off
    box off; ...ylabel('tCO_2e/MWh')
%     set(gca,'FontSize',15)


% load('Load2.mat','Load2'); % MISO
% % load('Load_NY_deter.mat','Load_NY_deter'); % NYISO
% 
% l = Load2(1,:);
% l = reshape(l,24,365);
% plot(1:24,l(:,1),'Color',[1,0.647,0,0.2]);
% hold on
% for day = 365:-1:2
%     disp(day)
%     if day>309
%     plot(1:24,l(:,day),'Color',[1,0.647,0,0.2]);
%     else
%     plot(1:24,l(:,day),'Color',[1,0,0,0.2]);
%     end
% end
% hold off
% box off
% set(gca,'FontSize',15)
% % set(gcf,'Position',[100 100 300 300])
% title('MISO load')
% axis([1 24 ylim])
% % saveas(gcf,strcat(Plotpath,'loadprofile_PV.png'));




















end