%% Patched Conics
% Code takes an initial and final date to search for a trajectory between
% Earth and Mars, then calculates a trajectory between Earth and SML1.
% Assumes that orbit about Earth and Mars is a circular orbit at an
% altitude of 400 km. For orbit about SML1, code assumes the orbit is
% merely at SML1, meaning it orbits about the Sun on the same orbital plane
% as Mars' orbit.
% Important note: Code needs lambert.m file to run correctly.
% Open Source Lambert problem solver is available at 
% https://github.com/rodyo/FEX-Lambert

close all
clear
t0 = juliandate(datetime('2022-01-01 00:00:00')); %Search start date
tf = juliandate(datetime('2028-01-01 00:00:00')); %Search end date
inc = 30; %Date search increment; ex: if inc = 7, searches every week
travmin = 150; %minimum allowed travel time
travmax = 250; %maximum allowed travel time
%Flight times were chosen based on plotting the following function after 
%the code finished running for a larger range of flight times and finding 
%a minimum spread of dV near this range
%plot3(eSML1_td,eSML1_ta-eSML1_td,eSML1_dv,'.')

k = 1; %counter
%For loop calculates starting from the start date up to the end date,
%inputing travel times starting from minumum to maximum travel times
for i = 0:inc:tf-t0
    for j = travmin:inc:travmax 
        em_dv(k) = PatchedConicsEtoM(t0+i,t0+i+j); %Records deltaV
        em_td(k) = t0+i; %Records launch date
        em_ta(k) = t0+i+j; %Records arrival date
        k = k+1; %updates counter
    end
end
[em_dV,em_i] = min(em_dv); %finds minimum delta v trajectory
em_TD = em_td(em_i); %min dV launch date
em_TA = em_ta(em_i); %min dV arrival date

datetime(em_TD,'convertfrom','juliandate')
datetime(em_TA,'convertfrom','juliandate')

k = 1; %counter
%For loop calculates starting from the start date up to the end date,
%inputing travel times starting from minumum to maximum travel times
for i = 0:inc:tf-t0
    for j = travmin:inc:travmax 
        eSML1_dv(k) = PatchedConicsEtoSML1(t0+i,t0+i+j); %Records deltaV
        eSML1_td(k) = t0+i; %Records launch date
        eSML1_ta(k) = t0+i+j; %Records arrival date
        k = k+1; %updates counter
    end
end
[eSML1_dV,eSML1_i] = min(eSML1_dv); %finds minimum delta v trajectory
eSML1_TD = eSML1_td(eSML1_i); %min dV launch date
eSML1_TA = eSML1_ta(eSML1_i); %min dV arrival date

datetime(eSML1_TD,'convertfrom','juliandate')
datetime(eSML1_TA,'convertfrom','juliandate')

function dv = PatchedConicsEtoM(te,tm)
%% Patched Conics Delta V Calculation
% Code takes a given Earth launch date and a given Mars arrival date and
% steps through calculations to find the total dV needed for the
% trajectory.
%% 0 Important parameters
GM_S = 1.32712440018E20/(1000^3); %All constants obtained from JPL Horizons
GM_E = 398600.435436;
GM_M = 42828.375214;
Re = 6371.01;
Rm = 3389.92;

he = 400; %Earth and Mars orbit altitudes, set arbitarily to 400 km
hm = 400;
%% 1 Lambert Problem
dt = tm-te; %Travel time
%Gathers Earth and Mars heliocentric positions at the given times
[re,ve,rm,vm,~,~] = EarthMarsSML1Positions(te,dt);
%Open Source Lambert problem solver, available at https://github.com/rodyo/FEX-Lambert
%Finds needed velocities for desired trajectories
[VE, VM, ~, exitflag] = lambert(re', rm', dt, 0, GM_S);
if exitflag ~= 1
    dv = NaN;
    return
end

%% 2
%Defines the required Earth velocity of the spacecraft relative to the Sun
V1PoS = norm(VE); 

%% 3
%Defines the velocity of Earth relative to the Sun
V1EoS = norm(ve);
%Finds the velocity of the spacecraft relative to Earth
V1PoE = V1PoS - V1EoS;

%% 4
%Sets velocity of the spacecraft relative to Earth equal to v-inf to find
%escape orbit
C1 = V1PoE^2/2;
V1esc = sqrt(2*(C1 + GM_E/(Re+he)));

%% 5
%Finds velocity of the spacecraft in Earth orbit
V1orb = sqrt(GM_E*(2/(Re+he) - 1/(Re+he)));

%% 6
%Finds dV of the first burn
dv1 = V1esc-V1orb;

%% 7
%Defines the required Mars velocity of the spacecraft relative to the Sun
V2PoS = norm(VM);

%% 8
%Defines the velocity of Mars relative to the Sun
V2MoS = norm(vm);
%Finds the velocity of the spacecraft relative to Mars
V2PoM = V2PoS - V2MoS;

%% 9
%Sets velocity of the spacecraft relative to Earth equal to v-inf to find
%escape orbit
C2 = V2PoM^2/2;
V2esc = sqrt(2*(C2 + GM_M/(Rm+hm)));

%% 10
%Finds velocity of the spacecraft in Mars orbit
V2orb = sqrt(GM_M*(2/(Rm+hm) - 1/(Rm+hm)));

%% 11
%Finds dV of the second burn
dv2 = V2esc-V2orb;

%% 12
%Sums both dVs of the burns
dv = dv2+dv1;

end

function dv = PatchedConicsEtoSML1(te,tSML1)
%% Patched Conics Delta V Calculation
% Code takes a given Earth launch date and a given SML1 arrival date and
% steps through calculations to find the total dV needed for the
% trajectory.
%% 0 Important parameters
GM_S = 1.32712440018E20/(1000^3); %All constants obtained from JPL Horizons
GM_E = 398600.435436;
Re = 6371.01;

he = 400; %Earth orbit altitudes, set arbitarily to 400 km
%% 1 Lambert Problem
dt = tSML1-te; %Travel time
%Gathers Earth and Mars heliocentric positions at the given times
[re,ve,~,~,rSML1,vSML1] = EarthMarsSML1Positions(te,dt);
%Open Source Lambert problem solver, available at https://github.com/rodyo/FEX-Lambert
%Finds needed velocities for desired trajectories
[VE, VSML1, ~, exitflag] = lambert(re', rSML1', dt, 0, GM_S);
if exitflag ~= 1
    dv = NaN;
    return
end

%% 2
%Defines the required Earth velocity of the spacecraft relative to the Sun
V1PoS = norm(VE); 

%% 3
%Defines the velocity of Earth relative to the Sun
V1EoS = norm(ve);
%Finds the velocity of the spacecraft relative to Earth
V1PoE = V1PoS - V1EoS;

%% 4
%Sets velocity of the spacecraft relative to Earth equal to v-inf to find
%escape orbit
C1 = V1PoE^2/2;
V1esc = sqrt(2*(C1 + GM_E/(Re+he)));

%% 5
%Finds velocity of the spacecraft in Earth orbit
V1orb = sqrt(GM_E*(2/(Re+he) - 1/(Re+he)));

%% 6
%Finds dV of the first burn
dv1 = V1esc-V1orb;

%% 7
%Defines the required SML1 velocity of spacecraft relative to the Sun
V2PoS = norm(VSML1);

%% 8
%Defines the velocity of SML1 relative to the Sun
V2SML1oS = norm(vSML1);

%% 9
%Finds dV of the second burn
dv2 = abs(V2SML1oS-V2PoS);

%% 10
%Sums both dVs of the burns
dv = dv1+dv2;

end

function [re,ve,rm,vm,rSML1,vSML1] = EarthMarsSML1Positions(t0,dt)
%Based on existing code from MAE 4060
%Originally the code was used for tracking positions of Earth and Venus, so
%certain parameter changes had to be made, namely changing Venus to Mars.
%Code now also finds the position of SML1.

%INPUTS
%   t0 - Initial time (JD)
%   dt - Cruise duration (days)
%
%OUTPUTS
%   re - 3x1 col vector - Heliocentric position of the Earth (in km) at t0
%   ve - Heliocentric velocity magnitude of Earth (in km/s) at t0
%   rm - 3x1 col vector - Heliocentric position of Mars (in km) at t0+dt
%   vm - Heliocentric velocity magnitude of Mars (in km/s) at t0+dt
%   rSML1 - 3x1 col vector - Heliocentric position of SML1 (in km) at
%   t0+dt
%   vSML1  - Heliocentric velocity magnitude of SML1(in km/s) at t0+dt

%Earth and Mars orbital elements
%All parameters obtained from JPL Horizons
a_e =  9.990472190345268E-01;       %sma (AU)
a_m = 1.523788234861207E+00;
e_e = 1.578453256417450E-02;        %eccentricity
e_m = 9.350867317100472E-02;
I_e = 2.365870868513173E-03*pi/180; %Inclination (rad)
I_m = 1.847930297853360E+00*pi/180; 
w_e = 2.525270024725202E+02*pi/180; %arg. of periapsis (rad)
w_m = 2.867104377270767E+02*pi/180;
O_e = 2.109010444915891E+02*pi/180; %long. of ascending node (rad)
O_m = 4.948938336201466E+01*pi/180;
t_p_e = 2459582.827623156831;       % time of periapsis passage (JD)
t_p_m = 2459752.043297053780; 
GM_e = 398600.435436;               %G*Mass of Planet (km^3/s^2)
GM_m = 42828.375214;
GM_s = 1.32712440018E20/(1000^3);   %G*Mass of sun (km^3/s^2)
kmAU = 149597870.700;               %1 AU in km 

tJ = 86400.0; %Seconds/Julian day
mue = (GM_s+GM_e)/(kmAU^3)*(tJ^2);
mum = (GM_s+GM_m)/(kmAU^3)*(tJ^2);

[xe,ye,ve] = kepler_2body(a_e,e_e,mue,t0,t_p_e); %Ecliptic with E at t0
[xm,ym,vm] = kepler_2body(a_m,e_m,mum,t0+dt,t_p_m); %Ecliptic with M at t0+dt

re = [xe,ye,0]';
rm = [xm,ym,0]';

re = heliocentricConv(re,I_e,w_e,O_e); %Converts to heliocentric using
rm = heliocentricConv(rm,I_m,w_m,O_m); %orbital parameters

%Need to calculate location of SML1
M_mars = .64171e24; %Mass of Mars [kg]
M_sun = 1988500e24; %Mass of Sun [kg]
mu = M_mars/(M_mars+M_sun); %CR3BP parameter
options = optimoptions('fsolve','Display','none');

Rsm = norm(rm); %Distance between Sun and Mars
r0 = Rsm*(M_mars/(3*M_sun))^(1/3); %Rough estimate using Hill Radius
f = @(r) -M_sun/((Rsm-r)^2)+(M_mars/(r^2)) + (M_sun/(M_sun+M_mars)*Rsm - r) ... 
    * (M_sun+M_mars)/(Rsm^3); %A bit more exact of a calculation
gamma = (1 - mu - fsolve(f,r0,options)); %Ratio to locate SML1
rSML1 = rm*gamma;
vSML1 = vm*gamma;

re = re(:)*kmAU;
ve = ve(:)*kmAU./tJ;
rm = rm(:)*kmAU;
vm = vm(:)*kmAU./tJ;
rSML1 = rSML1(:)*kmAU;
vSML1 = vSML1(:)*kmAU./tJ;
end
    
function r_h = heliocentricConv(r,I,w,O)
%Heliocentric conversion using orbit parameters
    Rw = [cos(w)  sin(w) 0;
          -sin(w) cos(w) 0;
          0       0      1];
    RI = [1 0 0;
        0 cos(I) sin(I);
        0 -sin(I) cos(I)];
    RO = [cos(O) sin(O) 0;
        -sin(O) cos(O) 0
        0 0 1];
    r_h = (Rw*RI*RO)'*r;
end

function [x,y,v] = kepler_2body(a,e,mu,t,Tp) 
% Two-body problem using Newton-Raphson inversion of Kepler's
% time equation
% Inputs:
%   a - semi-major axis
%   e - eccentricity
%   mu - gravitational parameter (with same distance units as a)
%   t - day
%   Tp - passage of periapsis
%
% Output:
%   x,y - Orbital radius components in the perifocal
%         frame (e,q) directions.

n = sqrt(mu/(a^3));       %mean motion
h = sqrt((1-e^2)*mu*a); %angular momentum
M = mod(n*(t-Tp),2*pi); %mean anomaly

%Newton-Raphson to find eccentric anomaly. Iterate to 
%machine-percision tolerance tol=eps(2*pi)
tol = eps(2*pi);
    if M/(1-e) < sqrt(6*(1-e)/e)
        EccAnom0 = M/(1-e);
    else
        EccAnom0 = (6*M/e)^(1/3);
    end
    EccAnomT(1) = EccAnom0;
    EccAnomT(end+1) = EccAnomT(end) - (M - EccAnomT(end) + e*sin(EccAnomT(end)))/(e*cos(EccAnomT(end)) - 1);
    while abs(EccAnomT(end)-EccAnomT(end-1)) > tol
        EccAnomT(end+1) = EccAnomT(end) - (M - EccAnomT(end) + e*sin(EccAnomT(end)))/(e*cos(EccAnomT(end)) - 1);
    end
    EccAnom = EccAnomT(end);
%calculate x and y positions in the perifocal frame using Kepler's
%equations
x = a.*(cos(EccAnom)-e);
y = a*sqrt(1-e^2).*sin(EccAnom);
x = x';
y = y';
%calcualte the x and y velocities in the perifocal frame
vx = mu/h.*sin(EccAnom);
vy = mu/h.*(e+cos(EccAnom));
v = (vx.^2+vy.^2)^(1/2);
end