pso_CParameters.nParticles=20; % population in PSO
pso_CParameters.max_iterations=2e5; % number of max iterations
%pso_CParameters.fnc='fitnessFun_DER';
pso_CParameters.noIterationsToGap=2e2; % iterations to wait for fit improve
pso_CParameters.minIterations = 2e2; % number of min iterations/epochs
pso_CParameters.velocityClampingFactor = 1; % value for max/min velocities eq
pso_CParameters.threshold=1e-9; % threshold for fitness improvement
pso_CParameters.tau=0.20;% learning parameter for mutation function

pso_CParameters.inertiaWeight=1;
pso_CParameters.memoryWeightC1=2;
pso_CParameters.groupWeightC2=2;

pso_CParameters.I_bnd_constr = 2; %Using bound constraints /is possible to change direct in DE
% 1 repair to the lower or upper violated bound (why)
% 2 rand value in the allowed range
% 3 bounce back

pso_CParameters.perScnEval=1; %0-1: 0-100
