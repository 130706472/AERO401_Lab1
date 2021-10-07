clc;clear;close all;

% DATA = readtable('SSWT_DATA_2021.xlsx');
DATA = load('DATA.mat').DATA;

t = table2array(DATA(:,1));

% Pitot tube pressure [psig]
P_pitot = table2array(DATA(:,2));

% Plenum pressure [psig] -> stagnation pressure of the test section
P_plen = table2array(DATA(:,3));

% Tank pressure [psig]
P_tank = table2array(DATA(:,4));

% Static pressure [psig] -> at the pitot tube
P_stat = table2array(DATA(:,5));

%% Problem 1
figure(1)
plot(t,P_pitot,'DisplayName', 'Pitot Tube Pressure')
hold on
plot(t,P_plen,'DisplayName', 'Plenum Pressure')
plot(t,P_tank,'DisplayName', 'Tank Pressure')
plot(t,P_stat,'DisplayName', 'Static Pressure')
grid on
legend()
xlabel('Time [s]')
ylabel('Pressure [psig]')

%% Problem 2

% Throat width [inch]
width_th = 1.02;
% Throat area [in^2] -> critical area
A_throat = width_th.*4.8;

% Test area [in^2] -> Target Area
A_test = 4.8^2; 

% Test Section Mach number
M_ts = critAinv(A_test/A_throat, 1.4);
M_ts = M_ts(2);

% To do: find how accurate is the result.

%% Problem 3





%% Functions (Gas Dynamic)
function [Mach] = critAinv(AoAst,gamma)
%   AoA: A/A_star, ratio of area and critical area
    syms M 
    eqn = AoAst == (1/M)*(((2/(gamma+1))*(1+((gamma-1)/2)*(M^2)))^((gamma+1)/(2*(gamma-1))));
    Mach = real(double(solve(eqn,M)));
end


