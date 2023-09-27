function [outputs] = solver(s)

% Takes data inputs describing an electricity system, 
% along with modeling choices and assumptions, and 
% tries to find the set of additions and retirements 
% that minimize the cost of the electricity system over 
% some period.

% Track runtime
tic

% Initialize grid data for selected assumptions 
[s] = initializer(s);

% Define objective function as function handle to
% parameterized_objective, which calls the 
% accounting model, cashflows
ObjectiveFunction = @(x) parameterized_objective(s,x);

if s.use_two_stage == 1
    s.opt_stage = 1;
    [~,~,~,~,pop_stage1] = ...
        ga(ObjectiveFunction,...
        s.stage1.nvars ...
        ,[],[],[],[], ...
        s.stage1.lb, ...
        s.stage1.ub,[], ...
        s.stage1.ga_options);
    s.opt_stage = 2;
    init_population = zeros(size(pop_stage1,1),s.num_tech*s.Periods);
    nt = s.num_tech;
    p = s.Periods;
    for t= 1:nt
        start = t*p-(p-1);
        finish = t*p;
        init_population(:,start:finish) = ...
            repmat(pop_stage1(:,t)/7,1,p);
    end
    s.ga_options = optimoptions( ...
        s.ga_options, ...
        'InitialPopulationMatrix',init_population);
    ObjectiveFunction = @(x) parameterized_objective(s,x);
end

% Run the solver
[OPT_x,system_cost,~,output] = ga(ObjectiveFunction,...
    s.nvars,[],[],[],[],s.lb,s.ub,[],s.ga_options);

% Store solver outputs
Optimal_DeltaMW = OPT_x;
outputs.Optimal_DeltaMW = Optimal_DeltaMW;
outputs.system_cost=system_cost;
outputs.generations = output.generations;

%% Export subsidy cost (government spending)

% The cost of the subsidy is estimated ex-post, and it
% is not considered by the solver, only its
% cost-reduction benefits are considered, in the
% accounting module. One can model many subsidies, then 
% sum their resulting optimal system cost and their 
% subsidy cost to compare between these subsidy 
% scenarios, then make an optimal choice by analysis 
% ("by hand").

Popbuild = zeros(s.num_tech,s.BPeriods*5);
OPT_x = reshape(OPT_x,s.num_tech,s.Periods);
peryear = OPT_x/5;
Popbuild(:,1:s.Periods*5) = repelem(peryear,1,5);
subsidy_cost=repmat(s.subsidies,1,s.BPeriods*5);
subsidy_cost(Popbuild<0)=0;
subsidy_cost=subsidy_cost.*Popbuild*1000;
subsidy_cost= sum(subsidy_cost);
real=(1/(1+s.interest)).^(0:s.BPeriods*5-1);
subsidy_cost=sum(subsidy_cost...
    .*repmat(real,1,s.monte))/10^9;
outputs.subsidy_cost = subsidy_cost;

%% Export runtimes
runtime = toc;
runtime = seconds(runtime);
runtime.Format = 'hh:mm';
outputs.runtime = runtime;