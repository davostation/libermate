clc;clear; close all
$% input ode
g = 9.8; % m/s^2
R = 6.37*10^6; % m
dy = @(t,v, x) -g* (R^2/ (R+x) ^2);
dx = @ (t,v,x) v;
%8 RK-4 method
t(1) = 0;
v(1) = 1500;
x(1) = 0;
h = 0. 001;
interval = (0, 180];
n = (interval (2) - interval (1) ) /h t 1; % no. of point
for i = l:n
% delta t
kl = dv(t (i), v(i), x(i)):
k2 = dv (t (i) + (h/2 ), v(i) + h * (kl/2) , x(i) t h
k3 = dv(t(i) + (h/2 ), v(i) + h * (k2/2 ), x(i) + h
k4 = dv(t(i) +h, v(i) + h*k3, x(i) + h*k3);
v(i+l) = v(i) + (h/6) * (kl + 2*k2 + 2*k3 + k4);
kl = dx(t(i), v(i), x(i));
k2 = dx (t(i) + (h/2), v(i) t h * (kl/2), x(i) + h
k3 = dx (t (i) + (h/2 ), v(i) + h * (k2/2) , x(i) t h
k4 = dx (t(i) +h, v(i) + h*k3, x(i) + h*k3) ;
x(i+1 ) = x(i) + (h/6) * (kl + 2*k2 + 2*k3 + k4 );
t(i+l ) = t(i) + h;
end
find the index when v = 0. At v=0, x wil1 be max.
v zero = find (v>0); v zero = v zero(end);
X max= x(v_zero)
$% compare with the ode 45 (MATLAB inbuilt func)
[t, y] = ode45 ( @myfunc, [0 160], [1500, 0]):
V = y(:,1) ;
x = y(:,2);
% find the index when v = 0. At v = 0, x will be max
v zero = find (v>0 ); v zero = v zero(end) ;
X max = x(V_zero)
function eq = myfunc(t,y)
V = y (1);
x = y (2);
g= 9. 8; % m/s^2
R= 6.37*10^6; % m
eq(1,1) = -g* (R^2/ (R+x) ^2);
eq(2,1) = v;
end

