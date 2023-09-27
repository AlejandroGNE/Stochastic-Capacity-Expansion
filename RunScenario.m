%% Simple power sector model to compare inaction vs mitigation scenarios

% locate yourself in the right directory
clear
rootdir= pwd;
%% Define project, set assumptions
s.region =                3;
s.Renewable_standard =    0;
s.Storage_standard =      0;
s.Renewable_subsidies =   0;
s.Storage_subsidies =     0;
s.grid_data =             2022;
s.new_plants_data =       0;
s.ramping =               0;
s.storage =               1;
s.price_maker =           0;
s.monte =                 1;
s.CVaR =                  1; 
s.mutation_shrink =       0.95;
s.fixed_costs_options =   4;
s.mutation_scale =        0.05;
s.Mesh_Tolerance =        1e-3;
s.PopulationSize =        600;
s.FunctionTolerance_PS =  1e-3;
s.MaxGenerations =        1000;
s.learning =              3;
s.fuel_data =             2020;
s.cpuspertask =           str2double(getenv('NUMBER_OF_PROCESSORS'))/2;
s.timelimit =             1000;
s.penalty_slack =         1.2;
s.alpha =                 0; % 106
s.slurmid =               string(randi(1e6));
s.plot_solvers=           1;
s.scenario_ID =           strcat(date,'_t');
s.devDebug=0;
%% Set up cluster
if s.devDebug==0
cl = parcluster();
slurmid = s.slurmid;
current_folder = cl.JobStorageLocation; 
storage_folder = strcat(current_folder,'/',slurmid);
mkdir(storage_folder);
cl.JobStorageLocation = storage_folder;
cores = feature('NumCores');
cl.NumWorkers = cores; 
parpool(cl,cores) 
end
%% Run model

% name your scenario
    new_ID= 'ref2';
%     new_ID = readmatrix(strcat(rootdir,'/','ID_count_C_S_V.csv'));
%     new_ID = new_ID+1;
%     scenario_ID = strcat(string(new_ID),"_",string(datetime('yesterday')));
    scenario_ID = strcat(string(new_ID),"_",date);   
    outputdir = strcat(rootdir,'\outputs\outputs_',scenario_ID);
    mkdir(outputdir);
    [outputs] = solver(s);
    save(strcat(outputdir,"/",scenario_ID,".mat"),"outputs","s")
%     clear outputs
%     close all

% close parallel pool
poolobj = gcp('nocreate');
delete(poolobj);