function [plant] =...
        read_data (plant,NG, Coal, Nuc,Oil,Bio,...
        cap1,name1,Cap_cost2,num_tech,year,...
        fixed_costs_options,s)
%% start    
warning off
plant=struct2table(plant);
EPA_FOM = s.EPA_FOM;

if fixed_costs_options > 3
    storage_FOM_options = fixed_costs_options;
    fixed_costs_options = 3;
end

% convert to categorical data types
Cap_cost2.Tech = categorical(Cap_cost2.Tech);
Cap_cost2.Fuel = categorical(Cap_cost2.Fuel);
Cap_cost2.Prime_Mover =...
    categorical(Cap_cost2.Prime_Mover);

% add emissions of new power plants
Cap_cost2.CO2eq_lbMWh(Cap_cost2.Fuel == 'COAL') =...
    0.09552 * 1/1000000 *...
    Cap_cost2.HR_btuKWh(Cap_cost2.Fuel == 'COAL') *...
    1000 * 2204.62;
Cap_cost2.CO2eq_lbMWh(Cap_cost2.Fuel == 'GAS') =...
    0.05306 * 1/1000000 *...
    Cap_cost2.HR_btuKWh(Cap_cost2.Fuel == 'GAS')...
    * 1000 * 2204.62;
Cap_cost2.CO2eq_lbMWh(Cap_cost2.Fuel == 'OIL') =...
    0.07315 * 1/1000000 *...
    Cap_cost2.HR_btuKWh(Cap_cost2.Fuel == 'OIL')...
    * 1000 * 2204.62;

%% add new plants
for cc = 1:num_tech
    plant.Plant_Name(end+1)= strcat('new_',string(...
        Cap_cost2.Tech(cc)),'_',string(year)); ...name1(cc);
    plant.Category(end)=Cap_cost2.Tech(cc);
    plant.Heatrate_BtukWh(end)=Cap_cost2.HR_btuKWh(cc);
    plant.Prime_Mover(end)=Cap_cost2.Prime_Mover(cc);
    plant.UB(end)=Cap_cost2.UB(cc);
    plant.LB(end)=Cap_cost2.LB(cc);
    plant.OnM_dolkW(end)=Cap_cost2.OnM_dolkW(cc);
    plant.Var_dolMWh(end)=Cap_cost2.Var_dolMWh(cc);
    plant.CO2eq_lbMWh(end)=Cap_cost2.CO2eq_lbMWh(cc);
    plant.OnLineYear(end)=year;
end

plant.Capacity_MW(end-num_tech+1:end)=...
    cap1(1:num_tech,1);

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

storage = plant.Category == "STORAGE";
plant.Var_dolMWh(storage) = 0;

try
if storage_FOM_options == 4
    % override fixed costs for storage
    storage = plant.Category == "STORAGE" &...
        plant.OnLineYear < 2020;
    plant.OnM_dolkW(storage) = 30.4;
end

if storage_FOM_options == 5
    % override fixed costs for storage
    storage = plant.Category == "STORAGE";
    plant.OnM_dolkW(storage) = 24.93;
end

if storage_FOM_options == 6
    % override fixed costs for storage
    storage = plant.Category == "STORAGE";
    plant.OnM_dolkW(storage) = 18;
end

if storage_FOM_options == 7
    % override fixed costs for storage
    storage = plant.Category == "STORAGE";
    plant.OnM_dolkW(storage) = 0;
end
catch
end
end
% 
% % override fixed costs for storage
% storage = plant.Category == "STORAGE";... &...
% %     plant.OnLineYear < 2020;
% plant.OnM_dolkW(storage) = 24.93;
% % plant.OnM_dolkW(storage) = 30.4;
% % used value suggested in Table 10 for small PHS plants
% % https://www.pnnl.gov/sites/default/files/media/file/PSH_Methodology_0.pdf
% % 2023_03_03 update: FOM varies with plant size but is around 24 which
% % matches the FOM of new storage plants in AEO, so I'll just use that

%% estimate marginal costs for new plants
% add carbon tax
try
    if s.carbonTax>0 
        years= 2020:5:2050;
        multipliers= zeros(1,7);
        multipliers(7)= 1;
        multipliers(2:6)= interp1([2020 2050],[0 1],2025:5:2045);
        thisYearMultiplier= multipliers(years==year);
        thisYearTax= s.carbonTax * thisYearMultiplier;
        plant.carbonTax = thisYearTax * plant.CO2eq_lbMWh;
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

% % temporary, for debugging new penalty approaches
% x = plant.Category == 'OWN';
% disp('read_data MC =')
% disp(plant.MC(x))

% % systematically set penalty TC (marginal cost) to be higher than the
% % costliest marginal cost in the existing system by a factor given by the
% % penalty slack
% if year == 2020
% x = plant.Category == 'OWN';
% costliest_VOM = max(plant.TC(~x));
% costliest_FOM = max(plant.OnM_dolkW(~x));
% plant.TC(x) = costliest_VOM * s.penalty_slack;
% s.penalty_FOM = costliest_FOM * s.penalty_slack;
% end
% % plant.TC(x)=s.unmet_demand_penalization;

%% roughly estimate LCOE for retirement criteria

% these capacity factors were estimated as the weighted
% average capacity factor for each technology category
% from running the dispatch model in the first period
GAS_CC_CF   = 0.714571;
GAS_CT_CF   = 0.015963;
COAL_CF     = 0.38116;
NUCLEAR_CF  = 1;
OIL_CF      = 0.005; % <- assumed 0.5%, but was zero
BIOMASS_CF  = 0.5916;

x = plant.Category == 'GAS_CC';
plant.LCOE_dolMWh(x) =...
    ((plant.OnM_dolkW(x)...
    .* plant.Capacity_MW(x) ...
    * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
    * 8760 * GAS_CC_CF)...
    ./ (plant.Capacity_MW(x) * 8760 * GAS_CC_CF);

x = plant.Category == 'GAS_CT';
plant.LCOE_dolMWh(x) =...
    ((plant.OnM_dolkW(x)...
    .* plant.Capacity_MW(x) ...
    * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
    * 8760 * GAS_CT_CF)...
    ./ (plant.Capacity_MW(x) * 8760 * GAS_CT_CF);

x = plant.Category == 'COAL';
plant.LCOE_dolMWh(x) =...
    ((plant.OnM_dolkW(x)...
    .* plant.Capacity_MW(x) ...
    * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
    * 8760 * COAL_CF)...
    ./ (plant.Capacity_MW(x) * 8760 * COAL_CF);

x = plant.Category == 'NUCLEAR';
plant.LCOE_dolMWh(x) =...
    ((plant.OnM_dolkW(x)...
    .* plant.Capacity_MW(x) ...
    * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
    * 8760 * NUCLEAR_CF)...
    ./ (plant.Capacity_MW(x) * 8760 * NUCLEAR_CF);

x = plant.Category == 'OIL';
plant.LCOE_dolMWh(x) =...
    ((plant.OnM_dolkW(x)...
    .* plant.Capacity_MW(x) ...
    * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
    * 8760 * OIL_CF)...
    ./ (plant.Capacity_MW(x) * 8760 * OIL_CF);

x = plant.Category == 'BIOMASS';
plant.LCOE_dolMWh(x) =...
    ((plant.OnM_dolkW(x)...
    .* plant.Capacity_MW(x) ...
    * 1000) + plant.TC(x) .* plant.Capacity_MW(x) ...
    * 8760 * BIOMASS_CF)...
    ./ (plant.Capacity_MW(x) * 8760 * BIOMASS_CF);

% overwrite coal ramp rates
% adding as try-catch-end so it works with previous
% models versions that don't define s.coalramp
try
x = plant.Prime_Mover == 'ST' &...
    plant.Category == 'COAL';
plant.UB(x) = s.coalramp; plant.LB(x) = s.coalramp;
catch
end

plant = sortrows(plant,'LCOE_dolMWh');
%         x = isnan(plant.TC);
%         plant(x,:)=[];
plant = table2struct(plant,'ToScalar',true);
end
