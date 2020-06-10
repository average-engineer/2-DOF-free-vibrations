clear all
close all
clc

%parameter definition
%type of solver being used for computing bounce and pitch
solver = 'anal';

%initial damping acting as starting point for the optimization function
damping_0 = 1000;
%boundary conditions
%lower bound of damping
damping_lower = 500;
%upper bound of damping
damping_upper = 50000;

%vector containing damoing values
damping = [damping_lower:100:damping_upper];

for ii = 1:length(damping)
    cost(ii) = cost_function(damping(ii),solver);
end

%plotting the absolute pitch and bounce motion (cost) vs damping
figure(1)
plot(damping,cost,'-*')
xlabel('Damping')
ylabel('Absolute Motion')


%finding the optimum value of damping bounded by the lower and upper bound for which the absolute motion is minimum 
% damping_optimum = fmincon(@(damping,cost)cost_function(damping,solver),damping_0,[],[],[],[],damping_lower,damping_upper);

%global optimization object
gs = GlobalSearch;

%creating an optimization problem
problem = createOptimProblem('fmincon','x0',damping_0,'objective',@(damping,cost)cost_function(damping,solver),'lb',damping_lower,'ub',damping_upper);

%running the global optimization solver on the optimization problem
[damping_optimum] = run(gs,problem)