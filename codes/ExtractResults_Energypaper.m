%%Extracting the results and graphics

clear;clc;close all; tTotalTime=tic; % lets track total computational time
addpath('CallDataBases') 
%addpath('DEalg') %Participants should add to the path the folder with the code of their algorithms 
noRuns = 30; % Number of trials here

for alg_test=1:6

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
         filename='SCSPaperR/Results_PSO_Sol';    
        
    case 2
        filename='SCSPaperR/Results_DE_Sol';
     case 3
        filename='SCSPaperR/Results_DEbest_Sol';
      
      case 4
        filename='SCSPaperR/Results_Vortex_Sol';
     
      case 5
        filename='SCSPaperR/Results_HyDE_Sol';
    
      case 6
        filename='SCSPaperR/Results_HyDEDF_Sol';
    
 end

load(filename,'caseStudyData','otherParameters','ResDB')

%Calculate min, max, average, std, time solutions
Values=zeros(noRuns,5);
   

for iRuns=1:noRuns 
    otherParameters.DSO=0.2;
   best=ResDB(iRuns).sol
   [S_val, other_info]=feval('fitnessFun_GM2019',best,caseStudyData, otherParameters)
    
   Values(iRuns,1)=S_val; %Fitness
   Values(iRuns,2)= other_info.Fit; %Remuneration
   Values(iRuns,3)= other_info.Penalty; %Penalties
   Values(iRuns,4)= ResDB(iRuns).tOpt; %Time
   Convergence(iRuns,:)=ResDB(iRuns).fitVector;  
   
end

% min, max, mean, std, time, mean remuneration, mean penalties
Table(Select_Algorithm,:)=[min(Values(:,1)) max(Values(:,1)) mean(Values(:,1)) std(Values(:,1)) mean(Values(:,2)) mean(Values(:,3)) mean(Values(:,4))./60]
Conver(Select_Algorithm,:)=mean(Convergence);
end

plot(Conver')
xlabel('Iterations')
ylabel('Fitness')
legend('PSO_{imp}','DE_{rand}', 'DE_{current-to-best}', 'VS', 'HyDE', 'HyDE-DF')
savefig('Fig_6_Conver.fig')


[value,ind_alg]=min(Table(:,1)); %Finding the best solution

switch ind_alg
   case 1
         filename='SCSPaperR/Results_PSO_Sol';    
        
    case 2
        filename='SCSPaperR/Results_DE_Sol';
     case 3
        filename='SCSPaperR/Results_DEbest_Sol';
      
      case 4
        filename='SCSPaperR/Results_Vortex_Sol';
     
      case 5
        filename='SCSPaperR/Results_HyDE_Sol';
    
      case 6
        filename='SCSPaperR/Results_HyDEDF_Sol';
end
load(filename)

for iRuns=1:noRuns 
    fitness(iRuns)=ResDB(iRuns).Fit_and_p(1);
end

[~,ind]=min(fitness);
x=ResDB(iRuns).sol;
%[Fit,Other_info] = fitnessFun_GM2019(x,caseStudyData,otherParameters)
 otherParameters.DSO=0.2;
[Fit,Other_info,NewA,NewB,NewTotal,Flex_agg,FlexA, FlexB] = fitnessFun_GM2019_retrieveProfiles(x,caseStudyData,otherParameters);



DSO_request=caseStudyData.parameterData.DSO_R
Baseline=otherParameters.Baseline_total
BLA=sum(otherParameters.Baseline_A)
BLB=sum(otherParameters.Baseline_B)+sum(otherParameters.Baseline_C)




