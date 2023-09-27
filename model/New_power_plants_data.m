function [Cap_small,Cap_small_storage] = New_power_plants_data(new_plants_data)

if new_plants_data == 0
    load('Cap_small.mat','Cap_small','Cap_small_storage');      % Load input parameters aggregated by technology (e.g. HR, cap cost, fuel cost, etc.)
elseif new_plants_data == 1
    load('Cap_small_old.mat','Cap_small','Cap_small_storage');   
end

