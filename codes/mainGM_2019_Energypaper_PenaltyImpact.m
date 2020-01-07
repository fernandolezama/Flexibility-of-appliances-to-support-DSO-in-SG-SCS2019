%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB source code for creating scenarios for the paper:
%% * Lezama et. al: Flexibility management model of home appliances to support DSO requests in smart grids. submitted to Sustainable Cities and Society. 2019
%% Contact ing.flezama@gmail.com for details in the implementation

clear;clc;close all; 
%addpath('DEalg') %Participants should add to the path the folder with the code of their algorithms 
noRuns = 20; % Number of trials here

for alg_test=3:3
tTotalTime=tic; % lets track total computational time

addpath('alg_HyDEDF')
% Default parameters
DEparameters

Select_Algorithm=alg_test
%1: DE algorithm


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
 
      case 3
        filename='SCSPaperR/Results_VS_PenaltyImpact';
        deParameters.I_strategy=3;
        deParameters.I_strategyVersion=1;
        deParameters.I_itermax= 2e5;

        deParameters.I_NP=20;
        deParameters.F_weight=0.9;
        deParameters.F_CR=0.5;   
   
   
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set other parameters
 No_solutions=deParameters.I_NP;
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
    for iRuns=1:noRuns %Number of trails
        
       
        switch iRuns
             case 1
                otherParameters.DSO=0.05;
             case 2
                otherParameters.DSO=0.1;
            case 3
                otherParameters.DSO=0.2;
            case 4
                otherParameters.DSO=0.4;
            case 5
                otherParameters.DSO=0.6;
            case 6
                otherParameters.DSO=0.8;
            case 7
                otherParameters.DSO=1;
        end
        
        
        
        tOpt=tic;
        rand('state',123)% ensure stochastic indpt trials

                [ResDB(iRuns).Fit_and_p, ...
                ResDB(iRuns).sol, ...
                ResDB(iRuns).fitVector, ...
                ResDB(iRuns).Best_otherInfo] =...
                HyDE(deParameters,caseStudyData,otherParameters,lowerB,upperB);
         
         
        ResDB(iRuns).tOpt=toc(tOpt); % time of each trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save the results and stats
       x=ResDB(iRuns).sol;
       [Fit,Other_info,NewA,NewB,NewTotal,Flex_agg,FlexA, FlexB,Rem_A,Rem_B,Penalty,Total] = feval('fitnessFun_GM2019_retrieveProfiles',x,caseStudyData,otherParameters);
       FlexA_M(iRuns,:)=FlexA;
       FlexB_M(iRuns,:)=FlexB;
       FlexTotal_M(iRuns,:)=abs(FlexA)+abs(FlexB);
       Rem_A_M(iRuns,:)=Rem_A;
       Rem_B_M(iRuns,:)=Rem_B;
       Penalty_M(iRuns,:)=Penalty;
       Total_M(iRuns,:)=Total;
       
       
         iRuns
    
       % Save_results
    end
tTotalTime=toc(tTotalTime); %Total time
%% End of MH Optimization

DSO_request=repmat(sum(abs(caseStudyData.parameterData.DSO_R)),iRuns,1);
FlexA_M=sum(abs(FlexA_M),2);
FlexB_M=sum(abs(FlexB_M),2);
FlexTotal_M=sum(FlexTotal_M,2);

Table=[DSO_request./1000 FlexA_M./1000 FlexB_M./1000 FlexTotal_M./1000  Rem_A_M Rem_B_M Penalty_M Total_M]


save(filename)

end
