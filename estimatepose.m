function [muout error] = estimatepose(Ein, K, XX, detection, niter)
%This file does optimisation using the exponential map for rotation
%parameterisation, now with local frame.
mu = [0 0 0 0 0 0]';

R = [];
R(:,:,1) = [0 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 0];
R(:,:,2) = [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0];
R(:,:,3) = [0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0];
R(:,:,4) = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
R(:,:,5) = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
R(:,:,6) = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];



npoints = size(XX,2);
nparams = 6;

J = zeros(2*npoints,nparams);
r = zeros(2*npoints,1);



error = 100000000000;
dp = 10000;
iter = 0;
lambda = 10;
projection = detection;
while iter < niter
    iter = iter + 1;

    Rn = R(:,:,1)*mu(1) + R(:,:,2)*mu(2) + R(:,:,3)*mu(3) + R(:,:,4)*mu(4) + R(:,:,5)*mu(5) + R(:,:,6)*mu(6);
    E = expm(Rn);
    XX2 = E*Ein*XX;
    E = eye(4,4);
     
        
    for i = 1:npoints
        X = XX2(1,i);
        Y = XX2(2,i);
        Z = XX2(3,i);
        Xn = E(1,1)*X + E(1,2)*Y + E(1,3)*Z + E(1,4);
        Yn = E(2,1)*X + E(2,2)*Y + E(2,3)*Z + E(2,4);
        Zn = E(3,1)*X + E(3,2)*Y + E(3,3)*Z + E(3,4);
        Xnn = K(1,1)*Xn + K(1,3)*Zn;
        Ynn = K(2,2)*Yn + K(2,3)*Zn;
        Znn = Zn;
        
        x = Xnn/Znn;
        y = Ynn/Znn;
        u = detection(1,i);
        v = detection(2,i);
        
        for j = 1:nparams
            [dXnn_dp dYnn_dp dZnn_dp] = expdiffXn(X,Y,Z,E,R(:,:,j),K);
            J(1+2*(i-1),j) = (dXnn_dp*Znn - dZnn_dp*Xnn)*(Znn^2);
            J(2*i,j) = (dYnn_dp*Znn - dZnn_dp*Ynn)*(Znn^2);
        end
        
        r(1+2*(i-1)) = (x-u);
        r(2*i) = (y-v);
    end
    
    error = norm(r)^2;
    
 
    left = J'*J + lambda*diag(diag(J'*J));
    right = J'*r;
    pn = left\right;
    
    nmu = mu - pn*dp;
   
    
    
    Rn = R(:,:,1)*nmu(1) + R(:,:,2)*nmu(2) + R(:,:,3)*nmu(3) + R(:,:,4)*nmu(4) + R(:,:,5)*nmu(5) + R(:,:,6)*nmu(6);
    E = expm(Rn);
    XX2 = E*Ein*XX;
    
    E = eye(4,4);
    
    for i = 1:npoints

        X = XX2(1,i);
        Y = XX2(2,i);
        Z = XX2(3,i);
        Xn = E(1,1)*X + E(1,2)*Y + E(1,3)*Z + E(1,4);
        Yn = E(2,1)*X + E(2,2)*Y + E(2,3)*Z + E(2,4);
        Zn = E(3,1)*X + E(3,2)*Y + E(3,3)*Z + E(3,4);
        Xnn = K(1,1)*Xn + K(1,3)*Zn;
        Ynn = K(2,2)*Yn + K(2,3)*Zn;
        Znn = Zn;
        
        x = Xnn/Znn;
        y = Ynn/Znn;
        projection(1,i) = x;
        projection(2,i) = y;
        u = detection(1,i);
        v = detection(2,i);
        r(1+2*(i-1)) = (x-u);
        r(2*i) = (y-v);
    end
    nerror = norm(r)^2;
    
    
    
      
    if nerror <= error
        mu = nmu;
        lambda = lambda * (1-0.1);
    else
        lambda = lambda * (1+0.1);
    end
    
    
    
    
    

  

    
    
    
    
    
    
   
    
%     clc;
%     display(error);
%     display(nerror);
%     display(lambda);
%     display(dp);
%     display(npoints);
%     display(iter);

%     f = figure(1);
%     clf(f);
%     hold on;
%     plot(detection(1,:),detection(2,:),'bx');
%     plot(projection(1,:),projection(2,:),'rx');

    
end
    
muout = mu;



function x = ForwardSub(a,b)
% The function solves a system of linear equations ax=b 
% where a is lower triangular by using forward substitution.
% Input variables:
% a The matrix of coefficients.
% b A column vector of constants.
% Output variable:
% x A colum vector with the solution.

n = length(b);
x(1,1) = b(1)/a(1,1);
for i = 2:n 
    x(i,1)=(b(i)-a(i,1:i-1)*x(1:i-1,1))./a(i,i);
end
    
    
    
    
    
    
