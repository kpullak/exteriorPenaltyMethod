clc;
clear all;


% declaring the variables x1 x2 and mu used in the code 
% (Note: mu is treated as a constant until the last step)
syms x1 x2 mu 

% starting point, chosen in such a way that it violates 1 of the 4
% constraints, i.e., it is outside the feasible boundary
x0 = [2;0];

% evaluating the constraint values at the the initial point x0
k = constraint_Functions(x0);
% returns the constraint functions at point (x1, x2) 
constraints = constraint_Functions([x1,x2]);


for i=1:4
    if(k(i)>0) 
        f =  Objective_Function([x1,x2]) + mu*constraints(i)^2;      
    end
end

% compute the gradient of the penalty function constructed above
h = gradient(f,[x1,x2]);

% the objective is to solve the equations given by the gradient
% in-terms of 'mu'; as we have two variables x1,x2 and two equations
% from gradient, it can be done like below:
eqs = [h(1),h(2)];
vars = [x1,x2];

% this will solve the given equations in terms of the variables as
% input by treating the rest of the terms e.g: 'mu' as constants.
[solx1,solx2] = solve(eqs,vars) 

% The final points (x1 & x2) are where the limit of the solution
% expressions found above with mu tends to infinity are evaluated:
x_1 = limit(solx1, mu, inf);
x_2 = limit(solx2, mu, inf);

formatSpec = 'The solution points are: x1* -> %4.5f  & x2* -> %4.5f\n';
fprintf(formatSpec,x_1,x_2)

% the minimum value of the objective function is:
min_value = Objective_Function([x_1,x_2]);

formatSpec = 'The minimum value of the function is -> %4.5f\n';
fprintf(formatSpec,min_value)

