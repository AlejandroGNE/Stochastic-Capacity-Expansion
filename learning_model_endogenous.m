function [learning_multipliers_per_year] = learning_model_endogenous(learning_rates,capacity_forecast,existing_capacity,Periods,Year,num_tech)

learning_coefficient = log(1-learning_rates)/log(2);
% existing_capacity_matrix = repmat(capacity_forecast(:,1)-capacity_additions(:,1),1,Periods);
existing_capacity_matrix = repmat(existing_capacity,1,Periods);
learning_multipliers_per_period = (capacity_forecast(:,1:Periods)./existing_capacity_matrix).^repmat(learning_coefficient,1,Periods); %2016$/Kw

% Extrapolating values
years_interpolated = Year:(Year)+(Periods*5);
years_interpolated([1 5:5:(Periods)*5]) = [];
years_interpolated(end) = [];
learning_multipliers_per_year = zeros(num_tech,Periods*5);
learning_multipliers_per_year(:,5:5:(Periods)*5) = learning_multipliers_per_period;
years_with_data = [Year Year+4:5:(Periods)*5+(Year)]; %to include the starting year
learning_multipliers_per_year(:,1) = 1;
for j=1:num_tech
learning_multipliers_per_year(j,(learning_multipliers_per_year(j,:)==0)) =...
    interp1(years_with_data,[1 learning_multipliers_per_period(j,:)],years_interpolated,'linear');
end