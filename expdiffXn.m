function [dXnn dYnn dZnn] = expdiffXn(X,Y,Z,E,R,K)
dE_dp = R*E;


dXn = dE_dp(1,1)*X + dE_dp(1,2)*Y + dE_dp(1,3)*Z + dE_dp(1,4);
dYn = dE_dp(2,1)*X + dE_dp(2,2)*Y + dE_dp(2,3)*Z + dE_dp(2,4);
dZn = dE_dp(3,1)*X + dE_dp(3,2)*Y + dE_dp(3,3)*Z + dE_dp(3,4);

dXnn = K(1,1)*dXn + K(1,3)*dZn;
dYnn = K(2,2)*dYn + K(2,3)*dZn;
dZnn = dZn;


end

