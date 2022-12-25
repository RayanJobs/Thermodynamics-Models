Pem =5.0
Peh = 5.0
Bc = 10^-5
Bt = 1.0
BBT = 15.62
gamma= 22.14
theta = 1.0

m = 0;
%%


% syms x(t) y(t)

% ode1 = diff(x,2)== Pem * (diff(x)+Bc*exp(gamma*y/(1+y))*x)
% ode2 = diff(y,2)== Peh*(diff(y) +Bt*Bc*exp(gamma*y/(1+y))*x + BBT*(theta-y))
% odes = [ode1;ode2];
% S = dsolve(odes)

%%
xmesh = linspace(0,1, 100);
tspan = linspace(0,100, 100);
sol = pdepe(m,@pdefun,@pdeic,@pdebc, xmesh, tspan);
plot(tspan,sol(:,100, 2))
xlabel('Time')
ylabel("Dimensionless Temperature")
%surf(xmesh,tspan,sol(:,:,2));
%surf(xmesh,tspan,sol(:,:,2));
xlabel('x'), ylabel('Time'), zlabel('Temperature')

u1_init = 0;
u2_init = 1;

function [c,f,s] = pdefun(x,t,u,dudx)
Pem =5.0;
Peh = 5.0;
Bc = 10^-5;
Bt = 1.0;
BBT = 15.62;
gamma= 22.14;
theta = 1.0;
c = [1; 1];
%c = [0; 0];
f =  [1/Pem; 1/Peh].* dudx + [-1; -1];
y = u(1) - u(2);
s = [-Bc*exp(gamma*u(2)/(1+u(2))) *u(1); -Bt*Bc*exp(gamma*u(2)/(1+u(2))) *u(1) + BBT*(theta-u(2))];
end

function u0 = pdeic(x)
Pem =5.0;
Peh = 5.0;
Bc = 10^-5;
Bt = 1.0;
BBT = 15.62;
gamma= 22.14;
theta = 1.0;
u0 = [0; 1];
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
Pem =5.0;
Peh = 5.0;
Bc = 10^-5;
Bt = 1.0;
BBT = 15.62;
gamma= 22.14;
theta = 1.0;
pl = [ul(1)-1; ul(2)];
ql = [-1/Pem; -1/Peh];
pr = [0; 0];
qr = [1; 1];
end

