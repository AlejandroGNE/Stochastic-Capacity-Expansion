function [plant_retire_data,storage_power] =...
        retire_data (plant_retire_data,...
        retire_this_much_of,retire_cost)

warning off

% Convert from structure to table
plant_retire_data =...
    struct2table(plant_retire_data); 

% Convert fuel data from string to categorical
retire_cost.Tech = categorical(retire_cost.Tech); 

% Sort plants by descending order of variable cost
plant_retire_data = sortrows(...
    plant_retire_data,'LCOE_dolMWh','descend'); 

% list of categories being retired
categories_being_retired = retire_cost.Tech; 

% for each category being retired
for technology = 1:length(categories_being_retired) 

    % extract all generators of that category
    retireable_generators =...
        plant_retire_data(...
        plant_retire_data.Category ==...
        categories_being_retired(technology),:); 

    % get the cumulative sum of the capacities 
    % (previously sorted from most to least costly)
    cumulative_retireable_capacity =...
        cumsum(retireable_generators.Capacity_MW,1);

    % locate row of cumulative capacities that is 
    % greater than the capacity being retired
    ind = find(cumulative_retireable_capacity>=...
        -retire_this_much_of(technology),1); 


    if ind>1 % if we must retire 2 plants or more...

        % from located plants, make zero the capacity 
        % of all plants but the last one
        retireable_generators.Capacity_MW(1:ind-1) = 0; 

        % retire remaining capacity from last plant
        retireable_generators.Capacity_MW(ind) =...
            retireable_generators.Capacity_MW(ind)+...
            (retire_this_much_of(technology)+...
            cumulative_retireable_capacity(ind-1)); 

        % from located plants, erase all of them but 
        % the last one
        retireable_generators(1:ind-1,:)=[]; 

    elseif ind==1 % if we must retire 1 plant only

        retireable_generators.Capacity_MW(ind) =...
            retireable_generators.Capacity_MW(ind)+...
            retire_this_much_of(technology);

    else % if we must retire all plants
        % this includes over-retirements; if a request
        % to retire more than exists makes it to this
        % point, nothing will really happen
        retireable_generators = [];

    end

    % erase all plants of category retired from 
    % original database
    plant_retire_data(plant_retire_data.Category ==...
        categories_being_retired(technology),:) = [];

    % append new database with retirements applied to
    % the original plants database
    plant_retire_data =...
        [plant_retire_data;retireable_generators];
end

% sum all storage capacity in this period 
% (exporting it for use by storage model)
storage_power = sum(plant_retire_data.Capacity_MW(...
    plant_retire_data.Category == 'STORAGE'));

% convert plants database with retirements applied back 
% to structure
plant_retire_data=table2struct(...
    plant_retire_data,'ToScalar',true);    
