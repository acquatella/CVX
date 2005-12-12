% SEPARATE_PT_POLY  Separating a point from a polyhedron
% Sec. 8.1.1, Boyd & Vandenberghe "Convex Optimization" 
% Joelle Skaf - 10/09/05
%
% The goal is to produce a hyperplane separating x0 and the polyhedron 
% defined as {x | Ax <= b}
%           minimize    mu'*x0 - b'*lambda
%                       A'*lambda = mu
%                       norm(mu)* <= 1
%                       lambda >= 0

cvx_quiet(true)
% Input data
randn('seed',0);
n  = 10;
m  = 2*n;
x0 = randn(n,1);
A  = randn(m,n);
b  = rand(m,1);

% CVX solution
fprintf(1,'Finding a separating hyperplane between the 2 polyhedra...');

cvx_begin
    variables muu(n) lambda(m)
    maximize ( muu'*x0 - b'*lambda )
    A'*lambda == muu
    norm(muu) <= 1
    lambda >= 0
cvx_end

fprintf(1,'Done! \n');

% Verification
disp('------------------------------------------------------------------');
disp('Note that 0 is in {x | Ax <= b} by construction...' );
disp('Verifying that x0 is separated from {x | Ax <= b} i.e. mu^T*x0 > 0');
disp([' mu^T*x0 = ' num2str(muu'*x0) ]);
