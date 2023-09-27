function [NG,Coal,Nuc,Oil,Bio] = fuel_prices(Year,fuel_data)

if fuel_data == 2020
    load('database_fuel.mat','fuel');
elseif fuel_data == 2016
    load('database_fuel_old.mat','fuel');
elseif fuel_data == 2005
    load('database_fuel_2005.mat','fuel');
elseif fuel_data == 2023
    load('database_fuel_2023.mat','fuel');
end
    ind = find(fuel.Year==Year,1);
    Coal = fuel.Coal(ind);
    Nuc = fuel.Nuc(ind);
    Oil = fuel.Oil(ind);
    Bio = fuel.Biomass(ind);
    NG = fuel.NG(ind);