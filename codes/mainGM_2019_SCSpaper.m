%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB source code for creating scenarios for the paper:
%% * Lezama et. al: Flexibility management model of home appliances to support DSO requests in smart grids. submitted to Sustainable Cities and Society. 2019
%% Contact ing.flezama@gmail.com for details in the implementation


clear;clc;close all; 
%addpath('DEalg') %Participants should add to the path the folder with the code of their algorithms 
noRuns = 30; % Number of trials here

for alg_test=1:6
tTotalTime=tic; % lets track total computational time

addpath('alg_HyDEDF')
% Default parameters
DEparameters

Select_Algorithm=alg_test
%1: DE; 2: DEbest; 3: VS; 4: HyDE; 5: HyDE-DF; 6: PSO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data base 
caseStudyData=callDatabase(4);
%4: 20 devices each (20*7=140) With assumptions done by GECAD in 16/12/2019

%caseStudyData.parameterData.DSO_R=caseStudyData.parameterData.DSO_R*0;
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load MH parameters (e.g., get MH parameters from DEparameters.m file)
% algorithm='DE_rand'; %'The participants should include their algorithm here'
% DEparameters %Function defined by the participant
% No_solutions=deParameters.I_NP;

%% Load MH parameters (e.g., get MH parameters from DEparameters.m file)
switch Select_Algorithm
    case 1
        filename='SCSPaperR/Results_DE_Sol';
        algorithm='DE'; %'The participants should include their algorithm here'
        deParameters.I_strategy=1;
        deParameters.I_itermax= 2e5;

        deParameters.I_NP=20;
        deParameters.F_weight=0.9;
        deParameters.F_CR=0.5; 
        No_solutions=deParameters.I_NP;
    case 2
        filename='SCSPaperR/Results_DEbest_Sol';
        deParameters.I_strategy=2;
        deParameters.I_itermax= 2e5;

        deParameters.I_NP=20;
        deParameters.F_weight=0.9;
        deParameters.F_CR=0.5;
        No_solutions=deParameters.I_NP;
      case 3
        filename='SCSPaperR/Results_Vortex_Sol';
        deParameters.I_strategy=3;
        deParameters.I_strategyVersion=1;
        deParameters.I_itermax= 2e5;

        deParameters.I_NP=20;
        deParameters.F_weight=0.9;
        deParameters.F_CR=0.5;
        No_solutions=deParameters.I_NP;
      case 4
        filename='SCSPaperR/Results_HyDE_Sol';
        deParameters.I_strategy=3;
        deParameters.I_strategyVersion=3;
        deParameters.I_itermax= 2e5;

        deParameters.I_NP=20;
        deParameters.F_weight=0.9;
        deParameters.F_CR=0.5;
        No_solutions=deParameters.I_NP;
      case 5
        filename='SCSPaperR/Results_HyDEDF_Sol';
        deParameters.I_strategy=3;
        deParameters.I_strategyVersion=2;
        deParameters.I_itermax= 2e5;

        deParameters.I_NP=20;
        deParameters.F_weight=0.9;
        deParameters.F_CR=0.5;
        No_solutions=deParameters.I_NP;
        
      case 6
        filename='SCSPaperR/Results_PSO_Sol';
        addpath('PSO_Calg')
        algorithm='IPSO5NP_Dlb'; %'The participants should include their algorithm here'
        PSO_Cparameters %Function defined by the participant
        No_solutions=pso_CParameters.nParticles; %Notice that some algorithms are limited to one individual
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set other parameters
 %No_solutions=deParameters.I_NP;
 otherParameters =setOtherParameters(caseStudyData,No_solutions);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set lower/upper bounds of variables 
[lowerB,upperB] = setVariablesBounds(caseStudyData,otherParameters);
otherParameters.lowerlimit=lowerB;
otherParameters.upperlimit=upperB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some parameters that can be modified by the user
otherParameters.DirectMEthod=2; %1:without direct repair 2:With direct repairs (No violations guarantee)
otherParameters.ensPenalty=100; % Penalty factor:insufficient generation / energy not supplied


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call the MH for optimizationclear 
ResDB=struc([]);
    parfor iRuns=1:noRuns %Number of trails
        tOpt=tic;
        rand('state',iRuns)% ensure stochastic indpt trials

        
        switch Select_Algorithm
             case 6
                 [ResDB(iRuns).Fit_and_p, ...
                ResDB(iRuns).sol, ...
                ResDB(iRuns).fitVector, ...
                ResDB(iRuns).Best_otherInfo] =...
                 PSO_C(pso_CParameters,caseStudyData,otherParameters,lowerB,upperB);
            otherwise
                [ResDB(iRuns).Fit_and_p, ...
                ResDB(iRuns).sol, ...
                ResDB(iRuns).fitVector, ...
                ResDB(iRuns).Best_otherInfo] =...
                HyDE(deParameters,caseStudyData,otherParameters,lowerB,upperB);
           
        end
        
        ResDB(iRuns).tOpt=toc(tOpt); % time of each trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save the results and stats
       %x=ResDB.sol;
       %[S_val, other_info]=feval('fitnessFun_GM2019',x,caseStudyData, otherParameters)
         iRuns
    end
tTotalTime=toc(tTotalTime); %Total time
%% End of MH Optimization
save(filename)




end
