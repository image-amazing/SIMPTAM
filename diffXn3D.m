function [dXnn dYnn dZnn] = diffXn3D(E,K,param)
dXn = E(1,1)*(param==1) + E(1,2)*(param==2) + E(1,3)*(param==3);
dYn = E(2,1)*(param==1) + E(2,2)*(param==2) + E(2,3)*(param==3);
dZn = E(3,1)*(param==1) + E(3,2)*(param==2) + E(3,3)*(param==3);



dXnn = K(1,1)*dXn + K(1,3)*dZn;
dYnn = K(2,2)*dYn + K(2,3)*dZn;
dZnn = dZn;


end

