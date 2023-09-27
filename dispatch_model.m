function [dispatch_outputs] =...
        dispatch_model(dispatch_inputs,s)
%% Unpack inputs
plant_tb = dispatch_inputs.plant_struct;
demand = dispatch_inputs.Load_d;
Var_Energy = dispatch_inputs.Var_Energy;
year = dispatch_inputs.year;
RPS_target = dispatch_inputs.RPS_target;
%% load data
plant_tb=struct2table(plant_tb);
% plant_tb(plant_tb.OnLineYear>year,:) = [];
plant_tb=sortrows(plant_tb,'TC');
technology = plant_tb.Category;

% keep all plants apart
plant_main = plant_tb; 
main_capacity = plant_main.Capacity_MW;

H = technology == 'HYDRO';
S = technology == 'SOLAR';
W = technology == 'WIND';
B = technology == 'STORAGE';

hydro_generation = ...
    sum(main_capacity(H))*Var_Energy.Hydro;
solar_generation = ...
    sum(main_capacity(S))*Var_Energy.Solar;
wind_generation = ...
    sum(main_capacity(W))*Var_Energy.Wind;

var_energy =...
    hydro_generation + ...
    solar_generation + ...
    wind_generation;

% residual load = total load - renewable generation
residual_load = demand'-var_energy;
dispatch_outputs.residual_load = residual_load;
oversupply = residual_load<0;

curtailment = -sum(residual_load(oversupply));
residual_load(oversupply) = 0;

H = find(H); 
S = find(S); 
W = find(W); 
B = find(B); 

% dispatcheable plants = 
% all plants - renewable plants - storage
plant_tb([H;S;W;B],:)=[]; 

% plant_renewables = plant_main([H;S;W],:);
% annual_renewable_generation =...
% table(unique(plant_renewables.Category),...
% [sum(hydro_generation);sum(solar_generation);...
% sum(wind_generation)],'VariableNames',...
% {'Category','Annual_generation_MWh'});

capacity = plant_tb.Capacity_MW; 
marginal_cost = plant_tb.TC;
technology = plant_tb.Category;
clearing_price = zeros(8760,1);
dispatch = zeros(size(capacity,1),8760);

% --- ramp vs no ramp differ from here ---

if s.ramping == 0
    for hour = 1:8760
        resid_load = residual_load(hour);
        stack = cumsum(capacity,1);
        Ind = find(stack>=resid_load,1);
        clearing_price(hour) = marginal_cost(Ind);
        %     emission_rate(k) = plant_tb.CO2eq_lbMWh(Ind);
        dispatch(1:Ind,hour) = capacity(1:Ind);
        if (Ind>2)
            dispatch(Ind,hour) = resid_load-stack(Ind-1);
        elseif resid_load == 0
            dispatch(Ind,1) = 0;
        elseif Ind == 1
            dispatch(Ind,1) = stack(Ind)-resid_load;
        end
    end
    clearing_price = clearing_price';
end

if s.ramping == 1
    %% Set up dispatch variables; run for hour zero

    % upper bound or maximum capacity dispatcheable
    ramp_up_rate = plant_tb.UB.*capacity;

    % lower bound or minimum capacity shed
    ramp_down_rate = plant_tb.LB.*capacity;

    % marginal_emission_rate = zeros(8760,1);
    cap = plant_tb.Capacity_MW;
    capacity = cumsum(cap,1);
    resid_load = residual_load(8760);
    last_dispatched =...
        find(capacity>=resid_load,1);
    dispatch(1:last_dispatched,1) = cap(1:last_dispatched);

    % dispatch only the remaining capacity (all residual
    % load less all plants except for the last plant
    % dispatched)
    if last_dispatched>2
        dispatch(last_dispatched,1) =...
            resid_load...
            -capacity(last_dispatched-1);
    elseif resid_load == 0
        dispatch(last_dispatched,1)=0;
    elseif last_dispatched==1
        dispatch(last_dispatched,1)=...
            capacity(last_dispatched)-resid_load;
    end

    dispatch_in_previous_hour=dispatch(:,1);

    % find plants affected by ramp rates
    ramp = find(ramp_down_rate);

    % ramp down rates for plants affected by ramping only
    ramp_down_rates = ramp_down_rate(ramp);

    last_dispatched = zeros(8760,1);
    capacity = plant_tb.Capacity_MW;

    fleet_size = size(ramp_down_rate);
    min_MWh = zeros(fleet_size);

    %% Run dispatch for all hours of the year

    for hour = 1:8760

        % upper bound is the plant's full capacity, or its
        % capacity in previous hour + its ramp up rate,
        % whatever is less
        max_MWh = min(capacity,...
            (dispatch_in_previous_hour+ramp_up_rate));

        % lower bound is zero (turn plant off), or its
        % capacity in previous hour - its ramp down rate
        % (i.e. must-dispatch capacity), whatever is more
        MWh_ramped = dispatch_in_previous_hour(ramp);
        must_generate = max(0,MWh_ramped-ramp_down_rates);
        min_MWh(ramp) = must_generate;

        % Capacity that must be dispatched because plants
        % would not be able to ramp down further even if we
        % wanted to
        min_MWh_total = sum(min_MWh);

        % maximum capacity dispatchable per plant
        plants_available = max_MWh-min_MWh;
        resid_load = residual_load(hour);

        if resid_load > min_MWh_total
            new_residual_load = resid_load...
                -min_MWh_total;
            cases = 1;
        elseif resid_load > 0
            % if residual load is smaller than the
            % must-dispatch capacity, then it will be met
            % by the must-dispatch capacity
            new_residual_load = min_MWh_total...
                -resid_load;
            plants_available = min_MWh;
            cases = 2;
        elseif resid_load == 0
            % if there is curtailment, no dispatch
            % needs to be run; clearing price is zero
            new_residual_load = 0;
            cases = 3;
        end

        plants_meeting_new_residual_load =...
            zeros(fleet_size);

        if new_residual_load > 0
            % if there is no curtailment:    

            % stack maximum capacity dispatcheable
            % per plant (they're already sorted
            % by ascending marginal cost)
            cumulative_dispatchable_capacity =...
                cumsum(plants_available,1);

            % find where the stack is greater than
            % the residual load i.e. marginal plant
            last_plant_dispatched =...
                find(cumulative_dispatchable_capacity...
                >= new_residual_load,1);

            if cases == 1
                plants_meeting_new_residual_load(...
                    1:last_plant_dispatched) =...
                    plants_available(1:last_plant_dispatched);

                %if we dispatch 2 or more plants, remove excess
                % capacity from last plant dispatched
                if (last_plant_dispatched>1)
                    plants_meeting_new_residual_load(...
                        last_plant_dispatched) =...
                        new_residual_load-...
                        cumulative_dispatchable_capacity(...
                        last_plant_dispatched-1);
                else % if we only dispatch the first plant
                    plants_meeting_new_residual_load(...
                        last_plant_dispatched,1) =...
                        cumulative_dispatchable_capacity(...
                        last_plant_dispatched)...
                        -new_residual_load;
                end
            end
        else 
            % if there is curtailment:
            last_plant_dispatched = 0;
        end

        % sum must-dispatch generation with
        % meet-residual-load generation
        current_dispatch = min_MWh...
            +plants_meeting_new_residual_load;

        % export dispatch of current hour
        dispatch(:,hour) = current_dispatch;
        dispatch_in_previous_hour = current_dispatch;

        if last_plant_dispatched > 0
            clearing_price(hour) =...
                marginal_cost(last_plant_dispatched);
            %         marginal_emission_rate(hour) =...
            %             emission_rates(last_plant_dispatched);
                    last_dispatched(hour) = last_plant_dispatched;
        end
    end
    clearing_price = clearing_price';
    clearing_price = round(clearing_price,8);
    % marginal_emission_rate = marginal_emission_rate';
end

%% Get total costs of operating the grid; emissions

annual_MWh_per_plant = sum(dispatch,2);

fixed_cost =...
    sum(main_capacity.*...
    plant_main.OnM_dolkW * 1000); 

if s.penalty_approach == 3
    annual_MWh_per_plant(technology == "OWN") =...
        annual_MWh_per_plant(technology == "OWN")^2;
end

variable_cost_fossils =...
    sum((annual_MWh_per_plant.*marginal_cost));

marginal_cost_all = plant_main.TC;

var_cost_solar = mean(marginal_cost_all(S)); 
if isnan(var_cost_solar), var_cost_solar = 0; end
var_cost_wind = mean(marginal_cost_all(W)); 
if isnan(var_cost_wind), var_cost_wind = 0; end
var_cost_hydro = mean(marginal_cost_all(H)); 
if isnan(var_cost_hydro), var_cost_hydro = 0; end

variable_cost_renewables = ...
    sum(Var_Energy.Solar) *...
    sum(main_capacity(S)) *...
    var_cost_solar +...
    sum(Var_Energy.Wind) *...
    sum(main_capacity(W)) *...
    var_cost_wind +...
    sum(Var_Energy.Hydro) *...
    sum(main_capacity(H)) *...
    var_cost_hydro;

penalty_FOM = 0;
penalty_size = max(dispatch((technology == "OWN"),:),[],'all');
if s.penalty_approach == 2
penalty_FOM = penalty_size * s.penalty_FOM * 1000;
end

% this section is to calculate the compensation we would need to apply in
% order for penalty plant to have the costliest LCOE in the system
% we get highest LCOE, and use this value in place of the VOM
% this leads to very high $/MWh but ensures penalty is priciest option

if year == 2016 % if running pre-dispatch
% weighted average CF whole year
plant_tb.CF = annual_MWh_per_plant ./ (capacity*8760);

% get LCOEs

% pre-allocate variable where LCOEs will be stored
plant_tb.LCOE_dolMWh_post(1) = 0;

% get LCOE for penalty plant
x = plant_tb.Category == 'OWN';
plant_tb.LCOE_dolMWh_post(x) =...
    ((plant_tb.OnM_dolkW(x)...
    .* penalty_size ...
    * 1000) + plant_tb.TC(x) .* penalty_size ...
    * 8760 * plant_tb.CF(x))...
    ./ (plant_tb.Capacity_MW(x) * 8760 * plant_tb.CF(x));

% get LCOE for all other plants
x = ~x;
plant_tb.LCOE_dolMWh_post(x) =...
    ((plant_tb.OnM_dolkW(x)...
    .* plant_tb.Capacity_MW(x) ...
    * 1000) + plant_tb.TC(x) .* plant_tb.Capacity_MW(x) ...
    * 8760 .* plant_tb.CF(x))...
    ./ (plant_tb.Capacity_MW(x) * 8760 .* plant_tb.CF(x));

x = plant_tb.LCOE_dolMWh_post == Inf | ...
isnan(plant_tb.LCOE_dolMWh_post);
plant_tb.LCOE_dolMWh_post(x) = 0;


highest_LCOE = max(plant_tb.LCOE_dolMWh_post);
dispatch_outputs.highest_LCOE = highest_LCOE;
% unmet_demand = sum(dispatch((technology == "OWN"),:));
% period_penalty = highest_LCOE * unmet_demand;
end

% total cost of operating the grid = 
% fixed costs + variable costs
total_operation_cost =...
    variable_cost_fossils +...
    fixed_cost +...
    variable_cost_renewables +...
    ...period_penalty +...
    penalty_FOM;
% disp('total_operation_cost');disp(total_operation_cost)
% disp('variable cost, fossils ='); 
% disp(variable_cost_fossils)
% disp('fixed costs ='); 
% disp(fixed_cost)
% disp('variable cost, renewables ='); 
% disp(variable_cost_renewables)

% GHG = sum(annual_MWh_per_plant.*...
%     emission_rates,'omitnan');
% hourlyGHG = sum(dispatch.*...
%     repmat(emission_rates,1,8760),1,'omitnan');

% temporary, for debugging new penalty approaches
% disp(plant_tb.TC(technology == "OWN")) % check VOM corresponds to penalty approach selected
% disp(penalty_FOM)

if RPS_target>0
    fossil_generation = sum(dispatch,'all');
    total_demand = sum(demand);
    percent_fossil_generation =...
        fossil_generation / total_demand;
    percent_renewable_generation =...
        1 - percent_fossil_generation;
    RPS_non_compliance_MWh =...
        ... % RPS target in MWh
        total_demand*RPS_target...
        ... % Actual Renewable generation
        -total_demand*percent_renewable_generation;

    Cost_of_RPS_non_compliance =...
        ... % Non-compliance in MWh
        RPS_non_compliance_MWh...
        ... % Penalty for non-compliance in $/MWh
        *s.RPS_non_compliance_penalty...
        ... % If renewable generation > RPS,
        ... % don't estimate any non-compliance cost
        * double(RPS_non_compliance_MWh>0);
    
    % add non-compliance costs to the total
    total_operation_cost = total_operation_cost+...
        Cost_of_RPS_non_compliance; 

    % export non-compliance costs
    dispatch_outputs.Cost_of_RPS_non_compliance =...
        Cost_of_RPS_non_compliance;
end
% disp('total_operation_cost');
% disp(total_operation_cost)

% set a breakpoint before this line and run code in
% "get_cap_factors.m" to estimate weighted average
% capacity factors for each technology category
%% pack outputs
dispatch_outputs.unmet_demand =...
    sum(dispatch((technology == "OWN"),:));
dispatch_outputs.total_operation_cost =...
    total_operation_cost;
dispatch_outputs.clearing_price = clearing_price;
dispatch_outputs.curtailment = curtailment;
% disp('Highest unmet demand (MW) =')
% disp(max(dispatch((technology == "OWN"),:),[],'all'))
%% only export if needed

% dispatch_outputs.variable_cost_renewables =...
%     variable_cost_renewables;
% dispatch_outputs.plant_tb = plant_tb;
% dispatch_outputs.plant_main = plant_main;
% dispatch_outputs.residual_load = residual_load;
% dispatch_outputs.last_dispatched = last_dispatched;
    
try
    dispatch_outputs.emission_rate =...
        marginal_emission_rate;
catch
    dispatch_outputs.emission_rate = zeros(1,8760);
end

try
    dispatch_outputs.hourlyGHG = hourlyGHG;
catch
    dispatch_outputs.hourlyGHG = zeros(1,8760);
end

if s.storage ==1 && s.price_maker == 0
    try
        % 3/22/2023:
        % export dispatchable capacity to
        % quick fix storage price taker
        % dispatchable capacity comprises idle MW
        % between marginal and penalty plant,
        % affected by their ramp rates
        dispatch_outputs.dispatchable_capacity = zeros(8760,1);
        for h = 1:8760
            dispatch_outputs.dispatchable_capacity(h) =...
                sum(...
                plant_tb.Capacity_MW(last_dispatched(h)+1:end-1).*...
                plant_tb.UB(last_dispatched(h)+1:end-1));
        end
    catch
    end
end

if s.report > 0

    dispatch_outputs.dispatch = dispatch;

    VOM_cost = sum(annual_MWh_per_plant.*...
    plant_tb.Var_dolMWh);
    dispatch_outputs.VOM_cost = VOM_cost;

    fuel_cost = sum(annual_MWh_per_plant.*plant_tb.MC);
    dispatch_outputs.fuel_cost = fuel_cost;

    dispatch_outputs.fixed_cost = fixed_cost;

    dispatch_outputs.plant_tb = plant_tb;
    
    dispatch_outputs.plant_main = plant_main;

    emission_rates = plant_tb.CO2eq_lbMWh;
    GHG = sum(annual_MWh_per_plant.*...
        emission_rates,'omitnan');
    dispatch_outputs.GHG = GHG;

    dispatch_outputs.hourly_unmet_demand =...
    dispatch((technology == "OWN"),:);

end
%% Export summary of generation per technology
pause_here=1;

if s.report > 0
    %% annual MWh sums

    % sum MWh of each plant for all hours of the year
    all_thermal= sum(dispatch,2);

    % group generation by technology group
    tablee= table(all_thermal,technology);
    summary= groupsummary(tablee,"technology","sum");
    
    % repeat for renewables, merge all
    tech_names= Var_Energy.Properties.VariableNames;
    MWh= [sum(hydro_generation) sum(wind_generation) sum(solar_generation)];
    count= [length(H) length(W) length(S)];
    summary_RW= table(tech_names',count',MWh',...
        'VariableNames',summary.Properties.VariableNames);
    all= [summary; summary_RW];
    all= renamevars(all,"sum_all_thermal","MWh");
    all= removevars(all,"GroupCount");
    all.year= repelem(year,height(all),1);
    all= movevars(all,'year','Before','MWh');
    dispatch_outputs.MWh_year= all;
end

% % thermal plants annual generation
% annual_generation_summary =...
%     table(categorical(technology),...
%     sum(dispatch,2),'VariableNames',...
%     {'Category','Annual_generation_MWh'}); 
% 
% % all plants
% annual_generation_summary =...
%     [annual_generation_summary;...
%     annual_renewable_generation]; 
% 
% MWh_per_technology =...
%     table(unique(annual_generation_summary.Category)...
%     ,'VariableNames',{'Technology'});
% for i = 1:length(MWh_per_technology.Technology)
%     MWh_per_technology.Annual_generation_MWh(...
%         MWh_per_technology.Technology == ...
%         MWh_per_technology.Technology(i)) = sum(...
%         annual_generation_summary.Annual_generation_MWh(...
%         annual_generation_summary.Category == ...
%         MWh_per_technology.Technology(i)));
% end
% 
% % summary of installed capacities for all technologies
% MWh_per_technology = sortrows(MWh_per_technology,...
%     "Annual_generation_MWh","descend"); 
% 
% MWh_per_technology.Technology(end+1) = 'TOTAL'; 
% MWh_per_technology.Annual_generation_MWh(end) =...
%     sum(MWh_per_technology.Annual_generation_MWh);
% for i = 1:length(MWh_per_technology.Technology)
%     MWh_per_technology.Percent_Share_of_Generation(i)...
%         = (MWh_per_technology.Annual_generation_MWh(i)...
%         / MWh_per_technology.Annual_generation_MWh(end))*100;
% end
% MWh_per_technology.Annual_generation_MWh =round(...
%     MWh_per_technology.Annual_generation_MWh(:),2);
% 
% data_str = string(...
%     MWh_per_technology.Percent_Share_of_Generation);
% for i = 1:numel(data_str)
%     data_str(i) = sprintf('%.2f',data_str(i));
% end
% MWh_per_technology.Percent_Share_of_Generation =...
%     data_str;
% 
% MWh_per_technology;
