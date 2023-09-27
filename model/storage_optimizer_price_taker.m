function...                         
    [...                            % Outputs
        NetLoad_after_storage,...   % 1
        Energy,...                  % 2
        optimal_state_of_charge,... % 3
        Power,...                   % 4
        profit...                   % 5
    ] = storage_optimizer_price_taker...
    (...                            % Inputs
        Clearing_price,...          % 1
        electricity_demand,...      % 2
        battery_capacity_MW,...     % 3
        dispatchable,...            % 4
        residual_load)              % 5
%         interval_size,...           % 6
%         overlap_size...             % 7       

% This function finds the optimal operation of the
% battery. We cannot optimize for a single year at a
% time, so the year must be partitioned into smaller
% intervals. Subsequent intervals are dependent on the
% previous one, so they are run serially. Because
% optimizing an interval is not identical to optimizing
% the full year, we add enough extra time to each
% interval, such that it emulates the operation of a
% full year for the portion of the interval. This extra
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

%% Define the size of the compartments
% tic

% after studying space of interval & overlap sizes, 
% decided on small interval with no overlap is close 
% enough to optimal (>99%)
interval_size = 120; % 120 hours = 5 days
overlap_size = 0;

% compartment size in hours
p_sz = interval_size + overlap_size; 

% position where intervals are stitched together
starting_hour = 1; 

% number of times the optimization will be run
number_of_compartments = ceil(8760/interval_size); 
total_hours_run =...
    number_of_compartments*interval_size...
    +overlap_size;
hours_past_8760 = total_hours_run-8760;

%% Properties of the battery

% Roundtrip efficiency
RTE = 0.80; 
round_trip_efficiency = sqrt(RTE);

% Number of hours it takes to attain full charge when 
% charged at max Energy capacity
time_to_full_charge = 4; 

% Max charge or Energy capacity of the battery
MWh = battery_capacity_MW*time_to_full_charge;
MW = battery_capacity_MW;

initial_SOC = MWh/4;

% Initial state of charge (how much charge batteries 
% had "brand new", before operations begin)
SOC_previous = initial_SOC; 

%% Extend variables to the size of all compartments
electricity_demand(8761:total_hours_run) =...
    electricity_demand(1:hours_past_8760); 

% rename clearing_price to preserve original data
clearing_price = Clearing_price; 

clearing_price(8761:total_hours_run) =...
    Clearing_price(1:hours_past_8760); 

dispatchable(8761:total_hours_run) =...
    dispatchable(1:hours_past_8760); 

residual_load(8761:total_hours_run) =...
    residual_load(1:hours_past_8760); 

%% Matrix to set constraints

charge = ones(p_sz,1)...
    .* round_trip_efficiency; 

discharge = -ones(p_sz,1)...
    ./ round_trip_efficiency;

SOC = zeros(2*p_sz,2*p_sz);

% for each hour of the year included in the compartment
for counter = 1:p_sz 

    % copy "charge" into each cell of the lower 
    % triangle of the diagonal of cutting the left half
    % of SOC matrix from left-top to bottom-right
    SOC(counter,1:counter) = charge(counter); 
    
    % copy "discharge" into each cell of the lower 
    % triangle of the diagonal of cutting the right 
    % half of SOC matrix from left-top to bottom-right
    SOC(counter,(p_sz+1):...
        counter+p_sz) = discharge(counter); 
end

% extend current SOC matrix by putting a negative SOC 
% below it
SOC(p_sz+1:...
    2*p_sz,1:2*p_sz) = ...
    -SOC(1:p_sz,1:(2*p_sz)); 

% this matrix is helpful to set the constraints on 
% state of charge for each hour based on the amount of 
% charge/discharge of storage in the previous hour
% it chains each hour with all the hours before it

%% Pre-allocate variables
Maxcharge = zeros(2*p_sz,1);
p_clearing_prices =...
    zeros(2*p_sz,1);
Energy = zeros(1,total_hours_run);
Power = zeros(1,total_hours_run);

% uncomment to store charge and discharge variables
%         X1 = zeros(1,total_hours_run); 
%         X2 = zeros(1,total_hours_run);

NetLoad_after_storage = zeros(1,total_hours_run);
optimal_state_of_charge = zeros(1,total_hours_run);

% 3/22/2023:
% bounds consider the grid capabilities to handle storage operation

% lower bound = zero i.e. do nothing
lb = zeros(2*p_sz,1); 

% upper bound is maximum charging/discharging per hour
ub = zeros(2*p_sz,1);
ub2 = zeros(p_sz,1);

%% Run serial optimizations of compartments

for iterations = 1:number_of_compartments
    
    % interval: number of hours for which storage 
        % operation is optimized at a time (hours)

    % overlap: number of hours added to interval for 
        % optimization but discarded later

    % compartment: interval + overlap
    p = starting_hour:starting_hour+p_sz-1;

    % upper bounds are based on the limits of the battery for
    % low-storage systems, and based on the limits of the grid for
    % high-storage systems
    
    % for buying/charging:
    % battery_limits = battery_capacity_MW - SOC
    % grid_limits = dispatchable_capacity_MW

    % if battery_limits > grid_limits, use grid_limits
    
    % dispatchable capacity from the grid i.e., how many MW can we charge
    % the battery based on the limits imposed by the grid
    grid_MW = dispatchable(p');
    
    % battery limit = MW
    batt_MW = repmat(MW,p_sz,1);

    % find hours where grid limits are/aren't reached
    grid_limits = batt_MW > grid_MW; % grid limits reached in these hours
    battery_limits = ~grid_limits; % grid limits not reached in these hours
    ub(grid_limits,1) = grid_MW(grid_limits); %
    ub(battery_limits,1) = batt_MW(battery_limits);

    net_load = residual_load(p');
    grid_limits = batt_MW > net_load;
    battery_limits = ~grid_limits;

    ub2(grid_limits,1) = net_load(grid_limits);
    ub2(battery_limits,1) = batt_MW(battery_limits);
    
    ub(p_sz+1:2*p_sz,1) = ub2;

    % initial charge is predefined for 1st iteration 
    % and estimated for 2nd iteration onwards, in MWh
    SOC_previous = repmat(SOC_previous,p_sz,1); 

    % How much energy can the battery take in? 
    % Can't buy/charge more energy than the total 
    % capacity of the battery, in MWh, minus what's 
    % already in it
    % This is the RHS constraint of the St contraint, 
    % Smax
    Maxcharge(1:p_sz,1) = repmat(MWh,p_sz,1) - SOC_previous; 

    % How much energy can the battery give out? 
    % Can't sell/discharge more energy than what's 
    % already in the battery, in MWh
    % concatenating constraints for charge and 
    % discharge limits for the SOC matrix
    Maxcharge(p_sz+1:2*p_sz,1) = SOC_previous; 
    
    % clearing prices are positive for buying/charging
    p_clearing_prices(1:p_sz,1) = clearing_price(p)'; 

    % and negative for selling/discharging
    p_clearing_prices(p_sz+1:2*p_sz,1) = -p_clearing_prices(1:p_sz,1); 

    % linprog internally reshapes f to column vector
    f = p_clearing_prices; 
    opts.Algorithm='dual-simplex';
    opts.Display='off';
    
    % Linprog problem is structured in two parts, one
    % for charging (upper half), one for discharging
    % (lower half).
    % x: decision variable, energy in and out, MWh
    % f: clearing prices, $/MWh
    % SOC: constraints
    try
        [x] = linprog(f,SOC,Maxcharge,[],[],lb,ub,[],opts);
    catch
        % newer matlab versions omit one of the options
        [x] = linprog(f,SOC,Maxcharge,[],[],lb,ub,opts);
    end

    if isempty(x)
        x = zeros(2*p_sz,1);
    end

    % for charge variables
    x1 = x(1:p_sz); 

    % for discharge variables
    x2 = x((p_sz+1):(2*p_sz)); 

    % Energy in/out each hour, units are MWh after 
    % considering the losses
    Energy(p) =...
        charge.*x1+(discharge.*x2); 

    % Hourly energy flows before affecting efficiency
    Power(p) = x1-x2; 

    NetLoad_after_storage(p) =...
        electricity_demand(p)' ...
        +Power(p);

    % the SOC of the current hour is the sum of Energy 
    % up to the current hour plus the initial SOC
    optimal_state_of_charge(p) =...
        cumsum(Energy(p)) ...
        +SOC_previous(1); 
    
    % prepare for next iteration:
    % going from the current to the next interval, 
    % subtracting the overlap
    starting_hour = starting_hour + interval_size; 

    SOC_previous =...
        optimal_state_of_charge(starting_hour-1);
end

%% Export variables
NetLoad_after_storage = NetLoad_after_storage(1:8760);

% for SOC accounting purposes
Energy = Energy(1:8760); 
% energy flows from the perspective of the battery; the
% battery loses part of the energy flowing in when it
% charges (hence Energy < Power for charging), and
% loses another part of the energy flowing out when it
% discharges (hence Energy > Power for discharging)

% for net load and profit accounting purposes
Power = Power(1:8760);
% energy flows from the perspective of the grid; energy
% purchases are accounted fully (no effiency loss) but
% energy sales are affected by efficiency loss

optimal_state_of_charge =...
    optimal_state_of_charge(1:8760);
profit = -(Power*Clearing_price');
% runtime = toc;
end