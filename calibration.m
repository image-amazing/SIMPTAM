function [ F ] = calibration(matches1, matches2,K)
%CALIBRATION Calculates the projection matrices of cameras from two sets of
% correspondences


T1 = normalise2d(matches1);
T2 = normalise2d(matches2);

[nmatches1, T11] = normalise2dpts(matches1);
[nmatches2, T21] = normalise2dpts(matches1);

nmatches1 = T1*matches1;
nmatches2 = T2*matches2;



A = [];
for i = 1:size(matches1,2)
    x = matches1(1,i);
    y = matches1(2,i);
    x_dash = nmatches2(1,i);
    y_dash = nmatches2(2,i);
    
    newrow = [x_dash*x x_dash*y x_dash y_dash*x y_dash*y y_dash x y 1];
    A = [A; newrow];
end

[U S V] = svd(A);
f = V(:,size(V,2));
F = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];



%Constraint Enforcement
[U S V] = svd(F);
S(3,3) = 0;
F = U*S*V';

%Denormalisation
F = T2'*F*T1;



% 
% 
% E = K'*F*K;
% t = null(E');
% t = t/t(1)*2;
% 
% [U, S, V] = svd(E);
% W = [0 -1 0; 1 0 0; 0 0 1];
% R1 = U*W*V';
% R2 = U*W'*V';
% 
% P = K*[1 0 0 0; 0 1 0 0; 0 0 1 0];
% P1 = K*[R1 t];
% P2 = K*[R2 t];
% 
% error1 = [];
% error2 = [];
% for i = 1:size(matches1,2)
%     X1 = linearreproject(matches1(:,i),matches2(:,i),P,P1);
%     X2 = linearreproject(matches1(:,i),matches2(:,i),P,P2);
%     error1 = [error1 X1(3)<0];
%     error2 = [error2 X2(3)<0];
% end
% 
% error1 = sum(error1);
% error2 = sum(error2);
% 
% 
% display(error1);
% display(error2);
% 
% 
% if (error1 < error2)
%     P2 = P1;
%     display(R1);
% else
%     display(R2);
% end
% end

