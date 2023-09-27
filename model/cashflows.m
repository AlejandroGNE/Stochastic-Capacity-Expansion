function [cashflows_outputs]=cashflows(s,x)
%% Unpack inputs
price_maker = s.price_maker;
monte = s.monte;
CVaR = s.CVaR;
Periods = s.Periods;
interest = s.interest;
learning = s.learning;
fuel_data = s.fuel_data;
BPeriods = s.BPeriods;
Year = s.Year;
Load2 = s.Load2;
Var_Energy = s.Var_Energy;
P0 = s.P0;
num_tech = s.num_tech;
name1 = s.name1;
plant = s.plant;
Cap_cost = s.Cap_cost;
subsidies = s.subsidies;
fixed_costs_options = s.fixed_costs_options;

%% switch opt stages
try
    if s.opt_stage == 1
        x = repmat((x')/7,1,Periods);
        x = reshape(x,1,num_tech*Periods);
    end
catch
end

%% trim over-retirements

% remove gas_cc_ret
Cap_cost2 = Cap_cost;
Cap_cost2(2,:) = Cap_cost2(end,:);
Cap_cost2(end,:) = [];
icap = -repmat(Cap_cost2.Retire_cap,1,Periods);

% work with temporary R&A for this section
x2 = reshape(x,num_tech,Periods);
x2(2,:) = x2(2,:) + x2(num_tech,:); % merge gas CC R&A
x2 = x2(1:end-1,:); % remove gas cc ret
xx = x2; % copy temp R&A to get cumulative R&A
xx = cumsum(xx')';

% new icap is the current icap plus R&A each period
% i.e., coal capacity in 2025 is coal capacity in 2020
% plus coal additions in 2020 minus coal retirements in
% 2020; repeat for all periods and technologies
new_icap = icap + xx;

% find when (which period) was the max retirements reached; values for
% later periods are zero
excess = zeros(size(new_icap));
excess(new_icap<0) = new_icap(new_icap < 0);
make_zero = zeros(size(excess));
num_techs = size(excess,1);
for i=1:num_techs
    n = find(excess(i,:),1);
    % if n is empty, the rest is just zero, always
    excess(i,n+1:end) = 0;
    make_zero(i,n+1:end) = ones(1,Periods-n);
end
make_zero = find(make_zero);
x2_trimmed = x2 - excess;
x2_trimmed(make_zero) = 0; %#ok<FNDSB> 

% excess_retirements = zeros(size(new_icap));
% excess_retirements(new_icap<0) = new_icap(new_icap < 0);
% x2_trimmed = x2 - excess_retirements; %   <-- new pop2
% xx_trim = cumsum(x2_trimmed')';
% new_icap_trim = icap + xx_trim; % no more negatives here (no ort)
cashflows_outputs.o.pop2 = x2_trimmed;

% if new_icap shows negative values, the model is
% over-retiring for those technologies/periods; we need
% to trim these, but keep the MW over-retired so
% we can penalize these and nudge model away from them

% penalize over-retirements by squaring their capacity
% ort_mw = new_icap(new_icap < 0);
% divide by 1e9 to bring penalty down to fitness
% function levels
% ortp = sum(ort_mw.*ort_mw)/1e8;
ortp = 0;

% % Check individual is within upper and lower bounds
% individual_is_out_of_bounds = double(...
%     ... % if there is at least one flag, population is
%     ... % not feasible
%     ~isempty(find(...
%     ... % flag if individual is not within bounds
%     ~(...
%     x>=lb & x<=ub... % is individual within bounds?
%     ), 1)));
% retire_more_than_existing_icap = zeros(1,Periods);
% 
% % Check we don't retire more than exists
% % I removed this later on, when I implemented the
% % smoothing, since retiring more than exists does not
% % really matter for cost purpose and we can take care
% % of it post optimization by literally trimming values
% % to realistic values
% 
% % Temporary variables that mix gas_CC R&As
% Cap_cost2 = Cap_cost;
% Cap_cost2(2,:) = Cap_cost2(end,:);
% Cap_cost2(end,:) = [];
% icap = -repmat(Cap_cost2.Retire_cap,1,Periods);
% 
% x2 = reshape(x,num_tech,Periods);
% x2(2,:) = x2(2,:) + x2(num_tech,:);
% x2 = x2(1:end-1,:);
% 
% % x2(1,2:num_tech:num_tech*Periods) =...
% %     x2(1,2:num_tech:num_tech*Periods) +...
% %     x2(1,num_tech:num_tech:num_tech*Periods);
% % x2(1,num_tech:num_tech:num_tech*Periods) = [];
% 
% % num_tech2 = num_tech - 1;
% 
% % xx = reshape(x2,num_tech2,Periods);
% xx = x2;
% xx = cumsum(xx')';
% new_icap = icap + xx;
% retire_more_than_existing_icap = new_icap >= 0;
% retire_more_than_existing_icap =...
%     double(~retire_more_than_existing_icap);
% retire_more_than_existing_icap =...
%     sum(retire_more_than_existing_icap);
% 
% % If individual is not feasible
% if sum(individual_is_out_of_bounds...
%         +retire_more_than_existing_icap) ~= 0
% %     % penalize with the typical lowest cost, plus the
% %     % sum of the square of the capacities, times alpha
% %     % this should make individuals different and nudge
% %     % solver to smaller R&A values, while making
% %     % unfeasible individuals always more expensive than
% %     % feasible ones
% %     capacity_leverager = sum(x.*x,"all")*(1/1e9);
% %     Total = 800 + s.alpha*capacity_leverager;
%     % penalize individual with a 2 trillion system cost
%     Total = 2000;
%     cashflows_outputs.y = Total;
%     return % end calculation
% end % otherwise, proceed with system cost calculation

%% 1 - Setup major variables
pop2 = reshape(x,num_tech,Periods);

gas_cc_retirements = pop2(end,:);
gas_cc_retirements_row = num_tech;
pop2(end,:) = [];
num_tech = num_tech - 1;

% replace over-retirements with zeros
pop2([1 3:num_tech],:) = x2_trimmed([1 3:num_tech],:);

capital_cost = zeros(num_tech,BPeriods*5);
Popbuild = zeros(num_tech,BPeriods*5);
peryear = pop2/5;
Popbuild(:,1:Periods*5) = repelem(peryear,1,5);
%% 2 - Learning Model

if learning ~= 2 % If not using endogenous learning

    if learning == 3 % no learning
        learning_multipliers_per_year = ones(num_tech,Periods*5);
    else
    global_capacity_forecast = P0(1:num_tech,2:end);
    % sources cited in Srujana's Thesis for P0:
    % [23] IEA: WEO Model,
    % https://www.iea.org/weo/weomodel/(2018)
    % [24] EIA: World Energy Projection System Plus:
    % Overview. 35 (2017)

    global_capacity_additions =...
        diff(P0(1:num_tech,:),[],2);

    learning_multipliers_per_year =...
        learning_model_exogenous(...
        Cap_cost.LR(1:num_tech),...
        global_capacity_forecast,...
        global_capacity_additions,...
        Periods,Year,num_tech);  
    end
elseif learning == 2 % If using endogenous learning

    existing_capacity =...
        -Cap_cost.Retire_cap(1:num_tech);

    % copy gas_CT because gas_CC is empty
    existing_capacity(2) =...
        -Cap_cost.Retire_cap(gas_cc_retirements_row); 

    capacity_additions_cumulative = pop2;
    capacity_additions_cumulative(...
        capacity_additions_cumulative<0) = 0;
    capacity_additions_cumulative =...
        cumsum(capacity_additions_cumulative,2);

    capacity_forecast = existing_capacity +...
        capacity_additions_cumulative;

    learning_multipliers_per_year =...
        learning_model_endogenous(...
        Cap_cost.LR(1:num_tech),...
        capacity_forecast,...
        existing_capacity,...
        Periods,Year,num_tech);
end
%% 3 - Capital costs with learning; subsidies
capex = Cap_cost.Capital_cost_dolKW(1:num_tech); %$/kW
subsidy_matrix = repmat(...
    subsidies(1:num_tech),1,Periods*5); %$/kW
subsidy_matrix(Popbuild(:,1:Periods*5)<0)=0;
cap_cost = repmat(capex,1,Periods*5);
cap_cost(Popbuild(:,1:Periods*5)<0) = 0; % make capex=0 for retirements
capital_cost(:,1:Periods*5) =... % apply learning multipliers
    cap_cost.* learning_multipliers_per_year;
gov_expenditures = capital_cost; % keep cap costs w/ learning
% if subsidy>capex (because we are using full or
% near-full subsidy and learning decreases capex), then
% we only pay for the actual capex subsidized
subsidy_limit = capital_cost(:,1:Periods*5)<subsidy_matrix; % instances where subsidy > capex
capital_cost(:,1:Periods*5) =...
    capital_cost(:,1:Periods*5)-subsidy_matrix; % subtract subsidies from capital costs
capital_cost(capital_cost<0) = 0;
capital_cost = capital_cost.*Popbuild*1000;
subsidy_cost=subsidy_matrix;
subsidy_cost(Popbuild(:,1:Periods*5)<0)=0;
% adjust excessive subsidies to match capex
subsidy_cost(subsidy_limit) =...
gov_expenditures(subsidy_limit);
% multiply by additions to get total gov expenditures
subsidy_cost=subsidy_cost.*...
Popbuild(:,1:Periods*5)*1000;
subsidy_cost = sum(subsidy_cost);
cashflows_outputs.subsidy_cost = subsidy_cost;
% disp('Capital costs of RW & storage')
% disp(capital_cost([5 6 8],1:35))

%% 4 - Retirements & additions (R&A)
pop2(end+1,:) = gas_cc_retirements;
Ind = 5:5:Periods*5;
year = Year+4:5:(Periods)*5+(Year);
y_chr = int2str(year');
y_chr = cellstr(y_chr);
y_chr = strcat('y',y_chr);
% s.last_period = year(end);
plant_struct(1:monte,1) =repmat(...
    table2struct(plant,'ToScalar',true),monte,1);
storage_power = zeros(monte,Periods);
for counter1=1:monte
    for counter2=1:Periods
        if (counter2~=1) % copy plant data from previous to current period
            plant_struct(counter1,counter2)=...
                plant_struct(counter1,counter2-1);
        end

        [NG,Coal,Nuc,Oil,Bio] =...
            fuel_prices(year(counter2),fuel_data);

        A=find(pop2(:,counter2)>0);

        plant_struct(counter1,counter2) =...
            read_data(... % <-- add new plants in R&A
            plant_struct(counter1,counter2),NG,...
            Coal,Nuc,Oil,Bio, pop2(A,counter2),...
            name1(A),Cap_cost(A,:),size(A,1),...
            year(counter2),fixed_costs_options,...
            s);

        B=find(pop2(:,counter2)<0);

        [plant_struct(counter1,counter2),...
            storage_power(counter1,counter2)] =...
            retire_data(... % <-- retire plants in R&A
            plant_struct(counter1,counter2),...
            pop2(B,counter2),...
            Cap_cost(B,:));
    end
end
plant_struct = plant_struct';
Load_d =reshape(permute(...
    Load2, [2 1 3]), size(Load2, 2), [])';
plant_struct=reshape(plant_struct,Periods*monte,1);
%% penalize non-compliance of storage standards
if s.Storage_standard > 0 && s.storage > 0

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
    
    % check if system complies every period
    SPS_non_compliance_MW =...
        SPS_target - storage_power;
    SPS_non_compliance_MW(SPS_non_compliance_MW<0) = 0;

    storage_row = Cap_cost.Tech == "Storage";

    % Penalty is the capital cost of storage itself
    s.SPS_non_compliance_penalty =... % dol/MW
        Cap_cost.Capital_cost_dolKW(storage_row)...
        * 1 ... % <-- scale up penalty here
        * 1000; % dol/kW * 1000 kW/MW = dol/MW
    
    Cost_of_SPS_non_compliance=...
        SPS_non_compliance_MW.* ...
        repelem(s.SPS_non_compliance_penalty,Periods);
end

RPS_target = zeros(1,Periods);
if s.Renewable_standard > 0

    % renewable standard applied to last period
    RPS_target(1,Periods) = s.Renewable_standard;

    % linearly interpolate standard to intermediate
    % periods
    RPS_target(1,2:Periods-1) = interp1(...
        [1 Periods],...
        [RPS_target(1) RPS_target(end)],...
        2:Periods-1);

end
%% 5 - Variable costs and emissions
Variable_cost = zeros(1,Periods*monte);
% GHG_period = zeros(1,Periods*monte);
curtailment = zeros(1,Periods*monte);
unserved_energy = zeros(1,Periods*monte);

% VOM_costs = zeros(1,Periods*monte);
% fuel_costs = zeros(1,Periods*monte);
% fixed_cost = zeros(1,Periods*monte);
% variable_cost_renewables = zeros(1,Periods*monte);

Cost_of_RPS_non_compliance = zeros(1,Periods*monte);
% dispatch_inputs.last_period = last_period;

% dd_in = dispatch model inputs
% dd_out = dispatch model outputs
% ss_in = storage model inputs
% ss_out = storage model outputs
% p = periods

for p = 1:Periods*monte
    dd_in.plant_struct = plant_struct(p,1);
    dd_in.Load_d = Load_d(p,:);
    dd_in.Var_Energy = Var_Energy;
    dd_in.year = year(p);
    dd_in.RPS_target = RPS_target(p);
%     if s.plot_solvers >0
%         disp(strcat('dispatch before storage, year ',string(dd_in.year)))
%     end
    [dd_out] = dispatch_model(dd_in,s);
    if s.report>0
        d.price_before_stor.(y_chr{p})= dd_out.clearing_price;
    end
%     if s.plot_solvers >0
%         disp(strcat('dispatch after storage, year ',string(dd_in.year)))
%     end
    if s.storage == 1
        ss_in = dd_in;
        ss_in.Clearing_price = dd_out.clearing_price;
        ss_in.storage_power = storage_power(1,p);
        ss_in.iter = 20;
        ss_in.s = s;
        try
        ss_in.dispatchable = dd_out.dispatchable_capacity;
        catch
        end
        ss_in.residual_load = dd_out.residual_load;
        if price_maker == 0
            [ss_out] = Price_taker(ss_in);
            dd_in.Load_d = ss_out.Netload_after_storage;
            [dd_out] = dispatch_model(dd_in,s);
        elseif price_maker == 1
            [ss_out] = Price_maker(ss_in);
            dd_out = ss_out.dispatch_outputs;
        end
    end
    if s.report > 0
        d.dispatch.(y_chr{p})= dd_out.dispatch;
        d.VOM_cost.(y_chr{p})= dd_out.VOM_cost;
        d.fuel_cost.(y_chr{p})= dd_out.fuel_cost;
        d.fixed_cost.(y_chr{p})= dd_out.fixed_cost;
        d.dispatched_plants.(y_chr{p})= dd_out.plant_tb;
        d.all_plants.(y_chr{p})= dd_out.plant_main;
        d.emissions.(y_chr{p})= dd_out.GHG;
        d.clearing_price.(y_chr{p})= dd_out.clearing_price;
        d.hourly_unmet_demand.(y_chr{p})= dd_out.hourly_unmet_demand;
        d.MWh_year.(y_chr{p})= dd_out.MWh_year;
        try
            d.storage_power.(y_chr{p}) = ss_out.power;
            d.storage_profit.(y_chr{p}) = ss_out.profit;
        catch
        end
    end

    Variable_cost(p) = dd_out.total_operation_cost;
    curtailment(p) = dd_out.curtailment;
    unserved_energy(p) = dd_out.unmet_demand;
%     GHG_period(counter) = dispatch_outputs.GHG;
%     variable_cost_renewables(counter) = dispatch_outputs.variable_cost_renewables;

    if RPS_target(p) > 0
        Cost_of_RPS_non_compliance(p) = dd_out.Cost_of_RPS_non_compliance;
    end
end

%% 6 - Annualize variable cost
Variable_cost = reshape(Variable_cost,Periods,monte);
year2 = Year:(Year)+(Periods*5);
year2([Ind,end]) = [];
Variable2 = zeros(BPeriods*5,monte);
Variable2(Ind,1:monte) = Variable_cost;
Variable2((Periods*5+1):end,1:monte) =...
    repmat(Variable2(Periods*5,:),20,1);

% Interpolate years not estimated from periodical 
% estimated data
Variable2(Variable2==0) =...
    interp1(year(1:Periods)',...
    Variable_cost,year2','linear');
Variable2(1:4,:) =...
    interp1(year(1:Periods)',...
    Variable_cost,year2(1:4),'linear','extrap');
% disp('Variable Costs ='); disp(sum(Variable2))
Variable3 = zeros(BPeriods*5,monte);

if s.Storage_standard > 0 && s.storage > 0

% Interpolate non-compliance costs from periods to
% years
all_years = (Year:Year+BPeriods*5-1)';
years_with_data = zeros(BPeriods*5,monte);
years_with_data(Ind) = all_years(Ind);
years_with_data_ind = years_with_data == all_years;
% years_without_data_ind = ~years_with_data_ind;
% years_without_data = all_years(years_without_data_ind);

SPS_withdata = SPS_non_compliance_MW ~= 0;
if sum(SPS_withdata) > 0
SPS_Ind = Ind(SPS_withdata);
years_with_data_SPS  = zeros(BPeriods*5,monte);
years_with_data_SPS(SPS_Ind) = all_years(SPS_Ind);
% years_with_data_SPS_ind =...
%     years_with_data_SPS == all_years;

years_without_data_SPS = all_years;
years_without_data_SPS(SPS_Ind) = 0; 
years_without_data_SPS((Periods*5+1):end) = 0; 
first_year_with_data_SPS = find(years_with_data_SPS,1);
years_without_data_SPS(1:first_year_with_data_SPS) = 0; 
years_without_data_SPS_ind =...
    years_without_data_SPS == all_years;

Variable3(years_with_data_ind,1:monte) =...
    Cost_of_SPS_non_compliance;

% copy value of last year to all subsequent years
Variable3((Periods*5+1):end,1:monte) =...
    repmat(Variable3(Periods*5,:),20,1);

% interpolate intermediate years
Variable3(years_without_data_SPS_ind) =...
    interp1(year(SPS_withdata)',...
    Cost_of_SPS_non_compliance(SPS_withdata)',...
    all_years(years_without_data_SPS_ind),'linear');
cashflows_outputs.Cost_of_SPS_non_compliance =...
    Cost_of_SPS_non_compliance;
end
end
%% 7 - Sum and discount capital and variable costs
Cost = repmat(sum(sum(capital_cost,3),1)',1,monte)...
    + Variable2 + Variable3;
real = (1/(1+interest)).^(0:BPeriods*5-1);
% disp('Capital costs =')
% disp(sum(repmat(sum(sum(capital_cost,3),1)',1,monte)))
% disp('Variable costs ='); disp(sum(Variable2))
% disp('SPS non-compliance cost ='); disp(sum(Variable3))
% Convert to US billions
Total = sum(Cost.*repmat(real',1,monte))/10^9;
%% 8 - Find Conditional Value at Risk (CVaR)
if monte>1
    Total = sort(Total,'asc');
    ind_cvar = round(CVaR*(length(Total)));
    Total = mean(Total(ind_cvar:end));
end

% store discounted total system cost separately
cashflows_outputs.DTSC = Total;

% add over-retirement penalty
Total = Total + ortp;
cashflows_outputs.ortp = Total;

% leverage capacity by multiplying some alpha factor
% times the sum of the square of all additions and
% retirements; this penalizes solutions that add too
% much in a single period and nudges model to solutions
% that spread out additions and retirements; since
% there are many local optima with very different
% buildouts, we are basically choosing a region of
% local optimas that spread out additions and
% retirements rather than lumping them
% if s.square_pop == 1
%     capacity_leverager = sum(pop2.*pop2,"all")*(1/1e9);
%     Total = Total + s.alpha*capacity_leverager;
% end

% 2023_03_09 update: determine penalty beforehand
% 2023_03_10 update: leave default option as we are
% likely not using penalty for dissertation (if fossil
% to renewable gap is revived, check code in lines
% 240-260 in smoothing_sweet_spot.m)
square_ras = (pop2.*pop2).*repmat(s.smooth,1,Periods);
capacity_leverager = sum(square_ras,"all")*(1/1e9);
Total = Total + s.alpha*capacity_leverager;

% % 2023_03_06 update: 
% % penalize energy-equivalent additions/retirements 
% solar_avg_CF = mean(Var_Energy.Solar);    
% wind_avg_CF = mean(Var_Energy.Wind);    
% square_ras = (pop2.*pop2);
% wind_row = string(Cap_cost.Tech) == "WIND";
% solar_row = string(Cap_cost.Tech) == "SOLAR";
% square_ras(wind_row,:) =...
%     square_ras(wind_row,:)*wind_avg_CF;
% square_ras(solar_row,:) =...
%     square_ras(solar_row,:)*solar_avg_CF;
% capacity_leverager = sum(square_ras,"all")*(1/1e9);
% Total = Total + s.alpha*capacity_leverager;
%% 9 - Export outputs
cashflows_outputs.y = Total;
cashflows_outputs.curtailment = sum(curtailment,'all');
try
cashflows_outputs.Cost_of_RPS_non_compliance =...
    Cost_of_RPS_non_compliance;
catch
end
cashflows_outputs.unmet_demand =...
    ...sum(unserved_energy);
    unserved_energy;
%% 10 - Report detailed outputs for deep analysis
if s.report > 0
    cashflows_outputs.o.dispatch = d.dispatch;
    cashflows_outputs.o.VOM_cost = d.VOM_cost;
    cashflows_outputs.o.fuel_cost = d.fuel_cost;
    cashflows_outputs.o.fixed_cost = d.fixed_cost;
    capcost = sum(capital_cost);
    capcost = capcost(1,1:35);
    capcost2 = zeros(1,35);
    for i = 5:5:35
        capcost2(i) = sum(capcost(i-4:i));
    end
    capcost = capcost2(Ind);
    cashflows_outputs.o.cap_cost = capcost;
    cashflows_outputs.o.dispatched_plants=...
        d.dispatched_plants;
    cashflows_outputs.o.all_plants = d.all_plants;
    cashflows_outputs.o.emissions = d.emissions;
    cashflows_outputs.o.clearing_price=...
        d.clearing_price;
    cashflows_outputs.o.hourly_unmet_demand=...
        d.hourly_unmet_demand;
    cashflows_outputs.o.MWh_year=...
        d.MWh_year;
%     cashflows_outputs.o.storage_energy=...
%         d.storage_energy;
%     cashflows_outputs.o.storage_soc=...
%         d.storage_soc;
try
    cashflows_outputs.o.storage_power=...
        d.storage_power;
    cashflows_outputs.o.price_before_stor=d.price_before_stor;
    cashflows_outputs.o.storage_profit=d.storage_profit;
catch % if no storage, ignore it
end
%     cashflows_outputs.o.storage_clearing_price=...
%         d.storage_clearing_price;
end
% cashflows_outputs.VOM_costs = VOM_costs;
% cashflows_outputs.fuel_costs = fuel_costs;
% cashflows_outputs.fixed_cost = fixed_cost;
% cashflows_outputs.variable_cost_renewables =...
%     variable_cost_renewables;
% disp('VOM_costs='); disp(VOM_costs)
% disp('fuel_costs='); disp(fuel_costs)
% disp('fixed_cost='); disp(fixed_cost)
% disp('RW_VOM_costs='); disp(variable_cost_renewables)
% try
% disp('Cost_of_RPS_non_compliance='); 
% disp(sum(Cost_of_RPS_non_compliance))
% catch
% end
% cashflows_outputs.emission_rate =...
%     dispatch_outputs.emission_rate;
% cashflows_outputs.dispatch =...
%     dispatch_outputs.dispatch;
% cashflows_outputs.hourlyGHG =...
%     dispatch_outputs.hourlyGHG;
% cashflows_outputs.plant_tb =...
%     dispatch_outputs.plant_tb;