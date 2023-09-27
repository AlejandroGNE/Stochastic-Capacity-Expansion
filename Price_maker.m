function [price_maker_outputs] = Price_maker(price_maker_inputs)
%% Unpack inputs
Load2015 = price_maker_inputs.Load_d';
ClrPrice = price_maker_inputs.Clearing_price;
cap = price_maker_inputs.storage_power;
iter = price_maker_inputs.iter;
s = price_maker_inputs.s;
dispatch_inputs = price_maker_inputs;

%% Run price maker model
counter=1;

Powerall = zeros(iter,8760);
loadresult = zeros(iter,8760);
clearprice = zeros(iter,8760);

% first three iterations with storage and clearing prices before average is considered
% tic
while counter<=3
        try 
            so_in= s.sIn;
        catch
        end
        so_in.ClrPrice = ClrPrice;
        so_in.Load2015 = Load2015;
        so_in.cap = cap;
        so_in.time_cluster = s.time_cluster;
    [so_out] = storage_optimizer(so_in);
        Loadstor= so_out.Loadstor;
        Power= so_out.Power;
        profit = so_out.profit;

    Powerall(counter,:) = Power;%units are MW in a given hour. This considers RTE of 80% in the storage
    loadresult (counter,:) = Loadstor; %for keeping a record of load after storage
    clearprice(counter,:) = ClrPrice; %for keeping a record of clearing price used for storage
    
    [dispatch_outputs] =  dispatch_model(dispatch_inputs,s);
    ClrPrice = dispatch_outputs.clearing_price;
    counter = counter + 1;
end
% time_to_run_first_three_iterations = toc

% Equal weighted average is calculated below, variables with suffix 2 indicate the variables used for storing average values
Powerall2(counter-3,:) = mean(Powerall(2:3,:),1);
% clrprice2(counter-3,:)=mean(clearprice(2:3,:),1); %mean of clearing price is considered too 
Loadstor = Load2015' + Powerall2(counter-3,:); %resultant load from average energy storage
loadresult2 (counter-3,:) = Loadstor; % keeping a record of resultant load from average energy storage in a sep variable

% average value in the above steps are considered in the loop below in each
% iteration average storage is calculated->clearing prices->resultant storage
% tic
while counter<=iter %can change the number of iterations manually
    dispatch_inputs.Load_d = Loadstor;
    [dispatch_outputs] =  dispatch_model(dispatch_inputs,s);
    ClrPrice = dispatch_outputs.clearing_price;

        try 
            so_in= s.sIn;
        catch
        end
        so_in.ClrPrice = ClrPrice;
        so_in.Load2015 = Load2015;
        so_in.cap = cap;
        so_in.time_cluster = s.time_cluster;
    [so_out] = storage_optimizer(so_in);
        Loadstor= so_out.Loadstor;
        Power= so_out.Power;
        profit = so_out.profit;

%     [Loadstor,~,~,Power,profit] = storage_optimizer(ClrPrice,Load2015,cap);
    Powerall(counter,:) = Power;
    clearprice(counter,:) = ClrPrice; %resultant energy storage is a resultant of clearing price from n-1th iteration
    loadresult (counter,:) = Loadstor;

    %calculating average values
    Powerall2(counter-2,:) = mean(Powerall(2:counter,:),1);
    Loadstor = Load2015' + Powerall2(counter-2,:); %for third iteration values
    loadresult2 (counter-2,:) = Loadstor;
    counter = counter + 1; %entering next iteration
end
% time_to_run_last_20_iterations = toc

% revenue=-Powerall2(1:end-1,:).*clearprice(4:end,:);

% returning value of the function
Loadresult = loadresult2(end,:);
PowerALL = Powerall2(end,:);
% Loadresult = Loadstor;
price_maker_outputs.Netload_after_storage = Loadresult;
price_maker_outputs.Storage_power = PowerALL;
price_maker_outputs.dispatch_outputs = dispatch_outputs;
price_maker_outputs.power = PowerALL;
price_maker_outputs.profit = profit;
end


   
    

    






  
   
    
      
                   
        
        
        
        

