% function [Fit_and_p,bestGlobalPos, fitMaxVector, objMaxVector, BotherParameters] = ...
%     PSO_C(psoParameters,caseStudyData,otherParameters,minPos,maxPos,initialSolution)

function [Fit_and_p,bestGlobalPos, fitMaxVector,other_info]=...
    PSO_C(psoParameters,caseStudyData,otherParameters,minPos,maxPos)

%Fixed memory parameters determined by the user

% Particle Swarm Optimization in MATLAB
% GECAD 2012
% Author: Joao Soares - MATLAB SOFTWARE ENGINEER & DEVELOPER
% JOAOBA2@GMAIL.COM
% new programming design of PSO - 2012
% PSO is a ...
% description of inputs and outputs
% to be added latter **

%IF 1 activated mutation of parameters
improved_PSO=1;


%signaling_flag=otherParameters.signActivated;
%otherParameters.fnc= otherParameters.WCCI_2020_funct;

%% init PSO basic variables
nParticles=psoParameters.nParticles;
nVariables=numel(maxPos);
psoParameters.nVariables=nVariables;
maxIterations=psoParameters.max_iterations;
%threshold=psoParameters.threshold;
psoParameters.maxPositions = maxPos;
psoParameters.minPositions = minPos;
% PSO variables fitMax
psoParameters.fitMax = inf; % inf for minimization problems. -inf for maximization problems
psoParameters.objMax = inf;
psoParameters.swarmFitnessIndBest = inf(nParticles,1);

BRM=psoParameters.I_bnd_constr;

% calculate velocity limits according to James Blondin PSO Tutorial 2009 - NOT USED IN QUANTUM-BEHAVED PSO
% vmax = k * (xmax - xmin) / 2 ; [-vmax; vmax]
psoParameters.maxVelocities = psoParameters.velocityClampingFactor * (maxPos-minPos) / 2 ;
psoParameters.minVelocities = - psoParameters.maxVelocities ;
% 
% for i = 1:OtherParameters.periods
%     idV2G = OtherParameters.ids.idsV2G ... 
%         + (PSOparams.nVariables/OtherParameters.periods)*(i-1);   
%     PSOparams.maxVelocities(idV2G) = PSOparams.maxVelocities(idV2G) * 0.05  ;
%     PSOparams.minVelocities(idV2G) = PSOparams.minVelocities(idV2G) * 0.05;
% end

% movement equation weight coefficients
% PSOparams.inertiaWeight = 1; % Kennedy et al 1995
% PSOparams.memoryWeightC1 = 2;  % Kennedy et al 1995
% PSOparams.groupWeightC2 = 2;  % Kennedy et al 1995

% vectorization of maxPositions and minPositions
maxPositionsMatrix = psoParameters.maxPositions(ones(1,nParticles),:);
minPositionsMatrix = psoParameters.minPositions(ones(1,nParticles),:);

% vectorization of max and min velocities
maxVelocitiesMatrix = psoParameters.maxVelocities(ones(1,nParticles),:);
minVelocitiesMatrix = psoParameters.minVelocities(ones(1,nParticles),:);

%% init PSO initial population and velocities
%rand('state',otherParameters.iRuns) %Guarantee same initial population
swarmPositions = unifrnd(minPositionsMatrix,maxPositionsMatrix,nParticles,nVariables);
if nargin>5
    noInitialSolutions = size(initialSolution,1);
    swarmPositions(1:noInitialSolutions,:)=initialSolution;
end

% %% Correct the solution to avoid selling to a cheaper price than the marginal cost
% T=caseStudyData.General.numPeriods; %Periods
% numType1=caseStudyData.General.numN1; %Consumers only
% numType2=caseStudyData.General.numN2; %Prosumers
% 
% ind_price=numel(minPos)/2;
% Q_lead=-1*swarmPositions(:,T*(numType1+numType2)+1:ind_price);
% b_chp=caseStudyData.Type4.MC; % cost factor for the CHP generation
% c_CHP=@(b_chp,g_chp) ((b_chp.*sqrt(g_chp))./g_chp+1e-10);
% Marginal_Cost=c_CHP(b_chp(1,1),Q_lead);
% 
% Corrected_bid= unifrnd(Marginal_Cost,maxPositionsMatrix(:,ind_price+1+T*(numType1+numType2):end))   ;
% ind_to_zero=isnan(Corrected_bid);
% Q_lead(ind_to_zero)=0;
% Corrected_bid(ind_to_zero)=caseStudyData.General.cf;
% 
% swarmPositions(:,T*(numType1+numType2)+1:ind_price)=-1*Q_lead;
% swarmPositions(:,ind_price+1+T*(numType1+numType2):end)=Corrected_bid;

swarmVelocities = unifrnd(minVelocitiesMatrix,maxVelocitiesMatrix,nParticles,nVariables);

if improved_PSO==1
    %%Improved PSO
    swarmGen = rand(nParticles, 4); % set initial swarmGen - strategic weights of movement equation
else
    swarmGen = [repmat([psoParameters.inertiaWeight, psoParameters.memoryWeightC1,psoParameters.groupWeightC2],nParticles,1), rand(nParticles, 1)]; % set initial swarmGen - strategic weights of movement equation
end

%% set initial global best and individual best
swarmBestPositions = swarmPositions;
bestGlobalPos = swarmPositions(1,:);

% store matrix for further use in other functions
psoParameters.maxPositionsMatrix=maxPositionsMatrix;
psoParameters.minPositionsMatrix=minPositionsMatrix;
psoParameters.maxVelocitiesMatrix=maxVelocitiesMatrix;
psoParameters.minVelocitiesMatrix=minVelocitiesMatrix;


%% PSO loop
% pre-allocation of loop variables
fitMaxVector = nan(1,maxIterations+1);
%objMaxVector = nan(3,maxIterations);
%alpha_initial = 0.3; % Quantum PSO alpha parameter
fitIterationGap = inf;
noIterationsToGap = psoParameters.noIterationsToGap;
iIteration = 1;
% count time of the loop only
%tPSOloop = tic;
while iIteration <= maxIterations+1 %%&&  fitIterationGap >= threshold
    
    % evaluate fitness for each particle
    [swarmPositions, swarmBestPositions, bestGlobalPos, psoParameters,other_info] = ...
        Evaluation(swarmPositions, swarmBestPositions, bestGlobalPos, caseStudyData, psoParameters, otherParameters);
    
    % mutation of strategic parameters of movement equation
    if improved_PSO==1
        %%Improved PSO
        [swarmGen] = Mutation(swarmGen, psoParameters);
    end
    
    % particles' movement
    [swarmPositions,swarmVelocities] = ...
        Movement(swarmPositions, swarmVelocities, swarmGen, swarmBestPositions, bestGlobalPos, psoParameters);
   
        
    % correct particles' positions if out of bounds
    swarmPositions = correctLimitSwarm(swarmPositions, psoParameters,BRM);
    

   
    
    % store fitness evolution and obj fun evolution as well
    fitMaxVector (iIteration) = psoParameters.fitMaxVector(1);

  %fprintf('Fitness value: %f\n',fitMaxVector(1,iIteration) )
  %  fprintf('Generation: %d\n',iIteration)
    
%     plotConvergence(iIteration,fitMaxVector)
%     pause(0.01) 
     iIteration = iIteration + 1;
    
end
%toc(tPSOloop)
p1=0;
Fit_and_p=[fitMaxVector(1,iIteration-1); 0]; %;p2;p3;p4]


%% M function
function[swarmGen] = Mutation(swarmGen, psoParameters)
swarmGen = swarmGen + randn(psoParameters.nParticles,4) * psoParameters.tau; % gaussian mutation
%swarm_gen=swarm_gen.*(lognrnd(0,1,4,num_particles).^tau); % lognormal mutation
swarmGen(swarmGen<0)=0;
swarmGen(swarmGen>1)=1;


%% correct limits functions
function[newVelocity] = correctVelocity(newVelocity, psoParameters)
% correct swarm vel to uppper and lower limits
[idx] = find(newVelocity>psoParameters.maxVelocitiesMatrix);
newVelocity(idx)=psoParameters.maxVelocitiesMatrix(idx);
[idx] = find(newVelocity<psoParameters.minVelocitiesMatrix);
newVelocity(idx)=psoParameters.minVelocitiesMatrix(idx);

function[newSwarmPos] = correctLimitSwarm(newSwarmPos, psoParameters,BRM)
switch BRM

    case 1 %Our method
        %[popsize,dim]=size(p);
       [idx] = find(newSwarmPos>psoParameters.maxPositionsMatrix);
        newSwarmPos(idx)=psoParameters.maxPositionsMatrix(idx);
        [idx] = find(newSwarmPos<psoParameters.minPositionsMatrix);
        newSwarmPos(idx)=psoParameters.minPositionsMatrix(idx);
    case 2 %Random reinitialization
        [idx] = [find(newSwarmPos<psoParameters.minPositionsMatrix);find(newSwarmPos>psoParameters.maxPositionsMatrix)];
        replace=unifrnd(psoParameters.minPositionsMatrix(idx),psoParameters.maxPositionsMatrix(idx),length(idx),1);
        newSwarmPos(idx)=replace;
end




%% Movement function
function[swarmPositions,swarmVelocities] = ...
    Movement(swarmPositions, swarmVelocities, swarmGen, swarmBestPositions,bestGlobalPos, psoParameters)
% movement function of PSO core
nParticles=  psoParameters.nParticles;
nVariables = psoParameters.nVariables;
% inertiaWeight=PSOparams.inertiaWeight;
% memoryWeightC1=PSOparams.memoryWeightC1;
% groupWeightC2=PSOparams.groupWeightC2;

% inertia component
inertiaTmp=swarmGen(:,1);
inertiaTmp = inertiaTmp(:,ones(1,nVariables));
inertia = swarmVelocities.*inertiaTmp;

% memory component
memoryTmp = swarmGen(:,2);
memoryTmp = memoryTmp(:,ones(1,nVariables));
memory = (swarmBestPositions - swarmPositions).*memoryTmp;

% group component
groupTmp = swarmGen(:,3);
groupTmp = groupTmp(:,ones(1,nVariables));
perturbationTmp = swarmGen(:,4); % perturbation of the global best
perturbationTmp = perturbationTmp(:,ones(1,nVariables)); % perturbation of the global best
bgMatrix = bestGlobalPos(ones(1,nParticles),:); 
bgMatrix = bgMatrix .* ( 1 + perturbationTmp .* randn(nParticles,nVariables)); % perturbation of the global best, miranda EPSO
group = (bgMatrix - swarmPositions).*groupTmp;

% velocity equation
swarmVelocities = inertia + memory + group;
% correct swarm velocities if out of bounds
swarmVelocities = correctVelocity(swarmVelocities, psoParameters);
% update particles' positions in space
swarmPositions = swarmPositions + swarmVelocities;

% % correct particles' positions if out of bounds
% swarmPositions = correctLimitSwarm(swarmPositions, PSOparams);

%% Evaluation function
function [swarmPositions, swarmBestPositions, bestGlobalPos, psoParameters,other_info] = ...
    Evaluation(swarmPositions, swarmBestPositions, bestGlobalPos, caseStudyData, psoParameters, otherParameters)


[swarmFitness, other_info]=feval(otherParameters.fnc,swarmPositions, caseStudyData, otherParameters);


[fitVal, I_best_index] = min(swarmFitness);
if fitVal < psoParameters.fitMax
    psoParameters.fitMax = fitVal;
    %psoParameters.objMax = swarmObjFun(idBestParticle,:);
    bestGlobalPos = swarmPositions(I_best_index,:);

 psoParameters.fitMaxVector=psoParameters.fitMax; %We save the mean value and mean penalty value

else
    psoParameters.fitMaxVector=psoParameters.fitMaxVector;
    %Best_otherInfo=Best_otherInfo;

end
%otherParameters = rmfield(otherParameters, 'pfDB');
idx = find(swarmFitness<psoParameters.swarmFitnessIndBest);
psoParameters.swarmFitnessIndBest(idx) = swarmFitness(idx);
swarmBestPositions(idx,:) = swarmPositions(idx,:);