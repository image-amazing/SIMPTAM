
detection = [2 2 1]';
X = [10 10 1 1]';
 
tx = 0;
ty = 0;
tz = 0;

for i = 1:100000

P_vanilla = [1 0 0 0; 0 1 0 0; 0 0 1 0];
Ext = [1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];
projection = P_vanilla*Ext*X;
projection = projection ./ projection(3);
x = detection(1);
y = detection(2);
px = projection(1);
py = projection(2);


error = (px-x)^2 + (py-y)^2;
display(error);




dExt_dtx = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
dExt_dty = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
dExt_dtz = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];

dp_tx = P_vanilla*dExt_dtx*X;
dp_ty = P_vanilla*dExt_dty*X;
dp_tz = P_vanilla*dExt_dtz*X;


dpx_dtx = dp_tx(1);
dpx_dty = dp_ty(1);
dpx_dtz = dp_tz(1);

dpy_dtx = dp_tx(2);
dpy_dty = dp_ty(2);
dpy_dtz = dp_tz(2);



de_tx = 2*(px-x)*dpx_dtx + 2*(py-y)*dpy_dtx;
de_ty = 2*(px-x)*dpx_dty + 2*(py-y)*dpy_dty;
de_tz = 2*(px-x)*dpx_dtz + 2*(py-y)*dpy_dtz;

tx = tx - de_tx*0.0001;
ty = ty - de_ty*0.0001;
tz = tz - de_tz*0.0001;

end






