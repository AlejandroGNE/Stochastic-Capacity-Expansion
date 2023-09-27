function [st,ops,changed] = checkpoint(ops,st,flag)

% File name where the checkpoint is stored
global scenario_ID outputdir
% filename = strcat(scenario_ID,'_checkpoint.mat');

% We defined scenario_ID, which gives this particular scenario a name or
% ID, and outputdir, where the output will be stored, as global variables,
% so that they can be defined from the beginning of the job and used
% throughout its execution. We do this because checkpoint.m is a function
% called from the matlab function ga, and I don't know how to pass these
% variables as inputs through ga. An alternative is to store scenario_ID
% and outputdir in csv files, and read them every time checkpoint.m is
% called, but I think keeping them as global variables rather than having
% to re-read the csv files is more efficient. Probably not that significant
% of a difference given the checkpoint.m function is only called every few
% generations, not as part of the fitness function.

%% determine if the function was called by ga or ps
try
    gaa= st.Generation; % if st.Generation exists, then it's ga
    gaa=1; ps=0;
catch
    % if the above fails, then it's pattern search
    ps=1; gaa=0;
end

if gaa
    filename = strcat(scenario_ID,'_checkpoint_ga.mat');
end

if ps
    filename = strcat(scenario_ID,'_checkpoint_ps.mat');
end

%% determine case
% 1) first-ever job for this scenario
% 2) resuming interrupted job for this scenario
% 3) running job trying to save progress -checkpoint
% 4) job completed

switch flag
    case 'init'
        if exist(filename,'file')
            Case=2; % resuming interrupted job for this scenario
        else
            Case=1; % first-ever job for this scenario
        end
    case 'iter'
        Case=3; % running job saving progress -checkpoint creation
    case 'done'
        Case=4; % job completed
end

%% execute cases for ga
if gaa
    switch Case
        case 1
            % do nothing; let first-ever job progress
        case 2
            % if resuming from a checkpoint
            % load checkpointed solver state
            % replace current with checkpointed state
            load(filename,'checkpoint');
            st.Generation = checkpoint.st.Generation;
            st.Best = checkpoint.st.Best;
            st.FunEval = checkpoint.st.FunEval;
            st.Expectation = checkpoint.st.Expectation;
            st.Selection = checkpoint.st.Selection;
            st.Score = checkpoint.st.Score;
            st.Population = checkpoint.st.Population;
            st.LastImprovement = checkpoint.st.LastImprovement;
        case 3
            % checkpoint every nth generation
            n=5;
            if ismember(st.Generation, n:n:1000)
                checkpoint.st = st;
%                 checkpoint.ops = ops;
%                 checkpoint.flag = flag;
                save(filename, 'checkpoint')
                copyfile(filename,outputdir)
                % the checkpoint file is saved to the local disk
                % a copy is also saved to the outputdir, outside of local disk,
                % so we can keep track of the job by looking at its checkpoints
            end
        case 4
            % if this job is done... what? delete the checkpoint? export other data? what?
    end
    changed= false;
end

%% execute cases for ps
if ps
    optimvalues=ops; options=st; % use ps lingo to avoid confusion
    switch Case
        case 1
            % do nothing; let first-ever job progress
        case 2
            % if resuming from a checkpoint
            % load checkpointed funccount
            % subtract funccount from MaxFunctionEvaluations
            load(filename,'checkpoint');
            new_options= options;
            past_funccount= checkpoint.optimvalues.funccount;
            new_MaxFunctionEvaluations= options.MaxFunctionEvaluations - past_funccount;
            new_options.MaxFunctionEvaluations= new_MaxFunctionEvaluations;
            ops= new_options;
            % we cannot alter the pattern search "iteration" like we can
            % alter "Generation" for ga, instead, to keep track of our
            % progress across jobs for the same scenario, every time a new
            % job resumes from a checkpoint, we check how many function
            % evaluations we have accumulated, and subtract that number
            % from the maximum; this yields a new maximum, which is the
            % stopping criteria often reached when using ps
        case 3
            % checkpoint every nth iteration
            n=3;
            if ismember(optimvalues.iteration, n:n:1000)
                checkpoint.optimvalues = optimvalues;
                save(filename, 'checkpoint')
                copyfile(filename,outputdir)
            end
        case 4
            % do nothing (for now)
    end
    st= false;
    switch Case
        case {1, 3, 4}
            changed= false;
        case 2
            changed= true;
    end
end
end