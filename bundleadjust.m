function [ outPTAM ] = bundleadjust(PTAM)
%BUNDLEADJUST Does bundle adustment on the PTAM model.

ncameras = size(PTAM.KeyFrames,2) - 1;
npoints = size(PTAM.Map.points,2);
% ncameras = 3;
% npoints = 6;


camparams = 6;
pointparams = 3;


R = [];
R(:,:,1) = [0 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 0];
R(:,:,2) = [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0];
R(:,:,3) = [0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0];
R(:,:,4) = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
R(:,:,5) = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
R(:,:,6) = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];





%Calculate the current error

r = zeros(2*ncameras*npoints,1);
param = zeros(6*ncameras + 3*npoints,1);

J = zeros(2*ncameras*npoints,6*ncameras + 3*npoints);
dp = 8000;
niter = 20;
iter = 0;
lambda = 10;
while iter < niter
    iter = iter + 1;
    
    row = -1;
    for i = 1:ncameras
        K = PTAM.KeyFrames(i+1).Camera.K;
        Ein = PTAM.KeyFrames(i+1).Camera.E;
        mu = param(3*npoints + 6*(i-1) + 1:3*npoints + 6*(i-1) + 6);
        Rn = R(:,:,1)*mu(1) + R(:,:,2)*mu(2) + R(:,:,3)*mu(3) + R(:,:,4)*mu(4) + R(:,:,5)*mu(5) + R(:,:,6)*mu(6);
        Ediff = expm(Rn);
        E = Ediff*Ein;
        
        for j = 1:npoints
            row = row + 2;
            id = PTAM.Map.points(j).id;
            
            x1 = findimagepoint(id,PTAM.KeyFrames(i+1));
            
            if ~isempty(x1)
                
                
                              
                XX = PTAM.Map.points(j).location;
                XX1 = XX;
                XX1 = XX1 + [param(3*(j-1) + 1: 3*(j-1) + 3); 0];                        
                XX2 = E*XX1;
                E2 = eye(4,4);
                Xn = XX2(1);
                Yn = XX2(2);
                Zn = XX2(3);
                
                Xnn = K(1,1)*Xn + K(1,3)*Zn;
                Ynn = K(2,2)*Yn + K(2,3)*Zn;
                Znn = Zn;
                
                x = Xnn/Znn;
                y = Ynn/Znn;
                
                u = x1(1);
                v = x1(2);
                
                r(row) = (x-u);
                r(row + 1) = (y-v);
                
                for p = 1:pointparams
                    [dXnn_dp dYnn_dp dZnn_dp] = diffXn3D(E,K,p);
                    J(row,p + 3*(j-1)) = (dXnn_dp*Znn - dZnn_dp*Xnn)*(Znn^2);
                    J(row + 1,p + 3*(j-1)) = (dYnn_dp*Znn - dZnn_dp*Ynn)*(Znn^2);
                end
                for c = 1:camparams
                    [dXnn_dp dYnn_dp dZnn_dp] = expdiffXn(Xn,Yn,Zn,E2,R(:,:,c),K);
                    J(row,3*npoints + c + 6*(i-1)) = (dXnn_dp*Znn - dZnn_dp*Xnn)*(Znn^2);
                    J(row + 1,3*npoints + c + 6*(i-1)) = (dYnn_dp*Znn - dZnn_dp*Ynn)*(Znn^2);
                end
                
                
            end
            
            
            
            
            
        end
    end
    
    error = norm(r)^2;
    
    left = J'*J + lambda*diag(diag(J'*J));
    right = J'*r;
    tic 
    pn = left\right;
    toc
    nparam = param - dp*pn;
    
    
    row = -1;
    for i = 1:ncameras
        for j = 1:npoints
            K = PTAM.KeyFrames(i+1).Camera.K;
            Ein = PTAM.KeyFrames(i+1).Camera.E;
            mu = nparam(3*npoints + 6*(i-1) + 1:3*npoints + 6*(i-1) + 6);
            Rn = R(:,:,1)*mu(1) + R(:,:,2)*mu(2) + R(:,:,3)*mu(3) + R(:,:,4)*mu(4) + R(:,:,5)*mu(5) + R(:,:,6)*mu(6);
            Ediff = expm(Rn);
            E = Ediff*Ein;
            row = row + 2;
            id = PTAM.Map.points(j).id;
            
            x1 = findimagepoint(id,PTAM.KeyFrames(i+1));
            
            if ~isempty(x1)
                XX = PTAM.Map.points(j).location;
                XX1 = XX;
                XX1 = XX1 + [nparam(3*(j-1) + 1: 3*(j-1) + 3); 0];                        
                XX2 = E*XX1;
                E2 = eye(4,4);
                Xn = XX2(1);
                Yn = XX2(2);
                Zn = XX2(3);
                
                Xnn = K(1,1)*Xn + K(1,3)*Zn;
                Ynn = K(2,2)*Yn + K(2,3)*Zn;
                Znn = Zn;
                
                x = Xnn/Znn;
                y = Ynn/Znn;
                
                u = x1(1);
                v = x1(2);
                
                r(row) = (x-u);
                r(row + 1) = (y-v);
            end
        end
    end
    nerror = norm(r)^2;
    
                
            
        
    if nerror <= error
        error = nerror;
        param = nparam;
        lambda = lambda * (1-0.1);
    else
        lambda = lambda * (1+0.1);
    end
    
   

    clc;
    display(error);
    display(iter);
    display(norm(param));
    
end


for i = 1:ncameras
    K = PTAM.KeyFrames(i+1).Camera.K;
    Ein = PTAM.KeyFrames(i+1).Camera.E;
    mu = param(3*npoints + 6*(i-1) + 1:3*npoints + 6*(i-1) + 6);
    Rn = R(:,:,1)*mu(1) + R(:,:,2)*mu(2) + R(:,:,3)*mu(3) + R(:,:,4)*mu(4) + R(:,:,5)*mu(5) + R(:,:,6)*mu(6);
    Ediff = expm(Rn);
    PTAM.KeyFrames(i+1).Camera.E = Ediff*Ein;
end
for i = 1:npoints
   PTAM.Map.points(i).location = PTAM.Map.points(i).location + [nparam(3*(i-1) + 1: 3*(i-1) + 3); 0];
end




outPTAM = PTAM;

end

function x = findimagepoint(id, KeyFrame)
x = [];
for i = 1:size(KeyFrame.ImagePoints,2)
    if (KeyFrame.ImagePoints(i).id == id)
        x = KeyFrame.ImagePoints(i).location;
    end
    
end

end
