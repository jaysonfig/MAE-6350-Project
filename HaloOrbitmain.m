%% Halo Orbit Calculation
% Code takes 3 initial conditions, xi, zi, and ydoti, and integrates the
% trajectory forward in time until it finds where the orbit has crossed the
% x-z plane. Once the 1/2 the period has been found, integrates the
% trajectory for one period and outputs stability parameters and a plot of
% the trajectory. If decent initial conditions are found, can use ziFixed.m
% or xiFixed.m to calculate corrections for the initial conditions. 

close all
clear
m1 = 1988500e24; %Mass of Sun kg
m2 = 0.64171e24; %Mass of Mars kg
M = m1+m2;
mu = m2/M; %CR3BP parameter
gamma = SML1(mu);
L1 = 1-mu-gamma; %Location of L1

xi    = 0.994755418540000;   %Initial conditions, starts off as guesses
zi    = 0.00152105500000000; %near L1. Only 3 parameters are needed since
ydoti = 0.00486076932000000; %y, xdot, and zdot = 0

C = conserved([xi,0,zi,0,ydoti,0],mu); %Jacobian constant
%Initial state, first 6 variables are position and velocity
%Next 36 variables are Psi matrix, which initially is the identity matrix
U0 = [xi;0;zi;0;ydoti;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;0;0;1;];
Ui = U0; %For the loop to change
inc = .1; %Time increment
tfin = 10; %Final time
tspan = 0:inc:tfin; %Time span
options = odeset('RelTol',1e-13,'AbsTol',1e-20);
%Loops until the trajectory crosses the x-z plane (y=0), iterates until it
%finds the time where y < 1e-11
ye = 1;
while ye > 1e-11
    [t,U] = ode45(@(t,U) DiffEoM(U,mu),tspan,Ui,options);
    for i = 2:length(U(:,2))
        if sign(U(i,2)*U(i-1,2)) < 0 %Takes the two points where y > 0 and 
            ti = t(i-1);             %y < 0, and sets these times as limits
            ii = i-1;                %of the time span for the next iteration
            tf = t(i);
            fi = i;
            break
        end
    end
    inc = (tf-ti)/10; %Shrinks time increment
    tspan = ti:inc:tf; %Sets time span
    Ui = U(ii,:); %Sets the state before crossing y to the initial state
    ye = U(ii,2);
    T12 = ti; %Saves half period
end
tspan = linspace(0,T12*2,1000); %Sets time to full orbit period
[~,U] = ode45(@(t,U) DiffEoM(U,mu),tspan,U0,options); %Integrates full orbit
UT = U(end,:); %Tracks end state of the orbit
M = Mcalc(UT,U0,mu); %Calculates M matrix to describe stability of orbit
[V,D] = eig(M);
lambda1 = D(1,1);
lambda2 = D(2,2);
lambda3 = D(3,3);
lambda4 = D(4,4);
nu1 = 1/2*(lambda1 + lambda2); %Two stability parameters, want magnitude
nu2 = 1/2*(lambda3 + lambda4); %to be less than 1 for stability

hold on %Plots the orbit in 3D
plot3(U(:,1),U(:,2),U(:,3)) %Orbit
plot3(1-mu,0,0,'ro') %Mars
%plot3(-mu,0,0,'go') %Sun
plot3(L1,0,0,'bo') %SML1
axis equal
hold off

% Differential Correction algorithms
% [deltazi,deltaydoti] = xiFixed(UT,mu);
% [deltaxi,deltaydoti] = ziFixed(UT,mu);
% Take outputs of xiFixed.m or ziFixed.m and add to appropriate initial
% conditions, ex. zi + deltazi = newzi

function r = SML1(mu)
% Finds location of L1 based on mu value
    r =(mu-3)^(1/3);%set starting value
    nmax=100;%set max number of iterations
	eps=1;%initialize error bound eps
    rvals=r;%initialize array of iterates
    n=0; %initialize n (counts iterations)
    while eps>=1e-14 && n<=nmax %set while-conditions
        f=r-(r^5-(3-mu)*r^4+(3-2*mu)*r^3-mu*r^2+2*mu*r-mu)/ ... 
            (2*mu-2*mu*r-3*r^2*(2*mu-3)+5*r^4+4*r^3*(mu-3));%compute next iterate
        rvals=[rvals;f]; %write next iterate in array
        eps=abs(f-r); %compute error
        r=f;n=n+1; %update x and n
    end  %end of while-loop
    r = real(r);
end

function C = conserved(X,mu)
%Calculates Jacobian constant using state and mu parameter
    x = X(:,1);
    y = X(:,2);
    z = X(:,3);
    xdot = X(:,4);
    ydot = X(:,5);
    zdot = X(:,6);
    d1 = sqrt((x+mu).^2+y.^2+z.^2);
    d2 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);
    U = 1/2*(x.^2+y.^2) + (1-mu)./d1 + mu./d2;
    C = 2*U-(xdot.^2+ydot.^2+zdot.^2);
end

function dUdt = DiffEoM(U,mu)
%Differential Equations of Motion, 42 first order diff eqs
    x = U(1);
    y = U(2);
    z = U(3);
    xdot = U(4);
    ydot = U(5);
    zdot = U(6);
    dUdx = x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) ...
        - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2));
    dUdy = y - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + ... 
        (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
    dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
        - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
    dUdt(1) = xdot;
    dUdt(2) = ydot;
    dUdt(3) = zdot;
    dUdt(4) = 2*ydot + dUdx;
    dUdt(5) = dUdy - 2*xdot;
    dUdt(6) = dUdz;
    X = [x;y;z;xdot;ydot;zdot];
    Psi(1,1:6) = U(7:12);
    Psi(2,1:6) = U(13:18);
    Psi(3,1:6) = U(19:24);
    Psi(4,1:6) = U(25:30);
    Psi(5,1:6) = U(31:36);
    Psi(6,1:6) = U(37:42);
    F = MatrixPartials(X(:),mu); %Calls different function for F matrix
    dPsidt = F*Psi;
    k = 7;
    for i = 1:6
        for j = 1:6
            dUdt(k) = dPsidt(i,j);
            k = k+1;
        end
    end
    dUdt = dUdt(:);
end

function F = MatrixPartials(X,mu)
%F matrix for Psi matrix differential equations
    x = X(1);
    y = X(2);
    z = X(3);
    Uxx = [(mu-1)/((mu+x)^2+y^2+z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1,(3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
           (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1, (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
           (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)];
    Lambda = [0 1 0;
             -1 0 0;
              0 0 0];
    F = [ zeros(3) , eye(3);
        Uxx , 2*Lambda];
end

function M = Mcalc(U,U0,mu)
% M matrix calculation to calculate orbit stability, based on Breakwell and
% Brown 1979
    x0 = U0(1);
    y0 = U0(2);
    z0 = U0(3);
    x0dot = U0(4);
    y0dot = U0(5);
    z0dot = U0(6);
    k = 7;
    for i = 1:6
        for j = 1:6
            Psi(i,j) = U(k);
            k = k+1;
        end
    end
    x = U(1);
    y = U(2);
    z = U(3);
    xdot = U(4);
    ydot = U(5);
    zdot = U(6);
    dUdx = x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) ...
        - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2));
    dUdy = y - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + ... 
        (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
    dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
        - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
    xddot = 2*ydot + dUdx;
    yddot = dUdy - 2*xdot;
    zddot = dUdz;
    dCdy0dot = -2*y0dot;
    dCdx0 = 2*x0 + ((2*mu + 2*x0)*(mu - 1))/((mu + x0)^2 + y0^2 + z0^2)^(3/2) - ...
        (mu*(2*mu + 2*x0 - 2))/((mu + x0 - 1)^2 + y0^2 + z0^2)^(3/2);
    dCdz0 = (2*z0*(mu - 1))/((mu + x0)^2 + y0^2 + z0^2)^(3/2) - ...
        (2*mu*z0)/((mu + x0 - 1)^2 + y0^2 + z0^2)^(3/2);
    M1 = [Psi(1,1) Psi(1,3) Psi(1,4) Psi(1,6);
        Psi(3,1) Psi(3,3) Psi(3,4) Psi(3,6);
        Psi(4,1) Psi(4,3) Psi(4,4) Psi(4,6);
        Psi(6,1) Psi(6,3) Psi(6,4) Psi(6,6)];
    M2 = 1/dCdy0dot*[Psi(1,5);Psi(3,5);Psi(4,5);Psi(6,5)]*[dCdx0 dCdz0 0 0];
    M3 = 1/y0dot*[0;0;xddot;zddot]*([Psi(2,1) Psi(2,3) Psi(2,4) Psi(2,6)] - ...
        Psi(2,5)/dCdy0dot*[dCdx0 dCdz0 0 0]);
    M = M1 - M2 - M3;
end

function [deltazi,deltaydoti] = xiFixed(U,mu)
% Algorithm that takes final state of a complete orbit, and outputs the
% changes in initial conditions to correct the orbit. Assumes xi will not
% change. Outputs changes needed for zi and ydoti
    x = U(1);
    y = U(2);
    z = U(3);
    xdot = U(4);
    ydot = U(5);
    zdot = U(6);
    dUdx = x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) ...
        - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2));
    dUdy = y - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + ... 
        (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
    dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
        - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
    xddot = 2*ydot + dUdx;
    yddot = dUdy - 2*xdot;
    zddot = dUdz;
    Psi(1,1:6) = U(7:12);
    Psi(2,1:6) = U(13:18);
    Psi(3,1:6) = U(19:24);
    Psi(4,1:6) = U(25:30);
    Psi(5,1:6) = U(31:36);
    Psi(6,1:6) = U(37:42);
    Mat1 = [Psi(4,3),Psi(4,5);Psi(6,3),Psi(6,5)];
    Mat2 = 1/ydot*[xddot;zddot]*[Psi(2,3),Psi(2,5)];
    delta = (Mat1 - Mat2)\[-xdot;-zdot];
    deltazi = delta(1);
    deltaydoti = delta(2);
end

function [deltaxi,deltaydoti] = ziFixed(U,mu)
% Algorithm that takes final state of a complete orbit, and outputs the
% changes in initial conditions to correct the orbit. Assumes zi will not
% change. Outputs changes needed for xi and ydoti
    x = U(1);
    y = U(2);
    z = U(3);
    xdot = U(4);
    ydot = U(5);
    zdot = U(6);
    dUdx = x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) ...
        - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2));
    dUdy = y - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + ... 
        (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
    dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
        - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
    xddot = 2*ydot + dUdx;
    yddot = dUdy - 2*xdot;
    zddot = dUdz;
    Psi(1,1:6) = U(7:12);
    Psi(2,1:6) = U(13:18);
    Psi(3,1:6) = U(19:24);
    Psi(4,1:6) = U(25:30);
    Psi(5,1:6) = U(31:36);
    Psi(6,1:6) = U(37:42);
    Mat1 = [Psi(4,1),Psi(4,5);Psi(6,1),Psi(6,5)];
    Mat2 = 1/ydot*[xddot;zddot]*[Psi(2,1),Psi(2,5)];
    delta = (Mat1 - Mat2)\[-xdot;-zdot];
    deltaxi = delta(1);
    deltaydoti = delta(2);
end