function [plant,Var_Energy,Load2] = Region_data(~,~)
     
    % load regional, grid-specific data

    % existing generation fleet
    load('sequoiaExistingFleet.mat','sequoiaExistingFleet')
    plant= sequoiaExistingFleet;
    
    % renewable capacity factors
    load('NYISO_2018.mat','Var_Energy')

    % electricity demand
    load('Load_NY_deter.mat','Load_NY_deter');
    Load2=Load_NY_deter;