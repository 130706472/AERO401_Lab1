clc;clear;close all;

DATA = load('DATA.mat').DATA;
global t P_pitot P_plen P_tank P_stat gamma

t = table2array(DATA(:,1));

% Pitot tube pressure [psig]
P_pitot = table2array(DATA(:,2));

% Plenum pressure [psig] -> stagnation pressure of the test section
P_plen = table2array(DATA(:,3));

% Tank pressure [psig]
P_tank = table2array(DATA(:,4));

% Static pressure [psig] -> at the pitot tube
P_stat = table2array(DATA(:,5));

gamma = 1.4;

% problem1()
% problem4()


%% Problem 1
function [] = problem1()
global t P_pitot P_plen P_tank P_stat
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
end

%% Problem 2
function [] = problem2()
global gamma
% Throat width [inch]
width_th = 1.02;
% Throat area [in^2] -> critical area
A_throat = width_th.*4.8;

% Test area [in^2] -> Target Area
A_test = 4.8^2; 

% Test Section Mach number
M_ts = critAinv(A_test/A_throat, gamma);
M_ts = M_ts(2);

% To do: find how accurate is the result.
end

%% Problem 3
function [] = problem3()
global t P_pitot P_plen gamma

    P02oP01 = (P_pitot + 14.7)./(P_plen + 14.7);
    M1_shock = [];
    % M2 ---> Mach number after shock, i.e. on the pitot tube
    M2_shock = [];
    % Calculate M2 in shock relations
    for i = 1:length(P02oP01)
        M1_shock = [M1_shock; stagPres(P02oP01(i), gamma)];
        M2_shock = [M2_shock; mach_change(M1_shock(i), gamma)];
        fprintf(1,['t = %2.3f \t M1 = %1.4f \t M2 = %1.4f\n'],t(i),M1_shock(i),M2_shock(i))
    end
end

%% Problem 4
function [] = problem4()
global t P_pitot P_stat gamma
    P02oP2 = (P_pitot + 14.7)./(P_stat + 14.7);
    M2_isen = [];
    M1_shock = [];
    for i = 1:length(P02oP2)
        M2_isen = [M2_isen; isentroPress(P02oP2(i), gamma)];
        M1_shock = [M1_shock; mach_change_inverse(M2_isen(i),gamma)];
        fprintf(1,['t = %2.3f \t M1 = %1.4f \t M2 = %1.4f\n'],t(i),M1_shock(i),M2_isen(i))
    end
end

%% Functions (Gas Dynamic)
function [Mach] = critAinv(AoAst,gamma)
%   Given ration of A/A_critical, ratio of area and critical area
    syms M 
    eqn = AoAst == (1/M)*(((2/(gamma+1))*(1+((gamma-1)/2)*(M^2)))^((gamma+1)/(2*(gamma-1))));
    Mach = real(double(solve(eqn,M)));
end

function [M1] = stagPres(P02oP01, gamma)
%   Given P_stag2/P_stag1, solve for upstream mach number
    syms M
    eqn = P02oP01 == (((gamma+1).*(M.^2))./((gamma-1).*(M.^2) + 2))...
        .^(gamma/(gamma-1)).*((gamma+1)./((2.*gamma.*M^2) - (gamma-1)))...
        .^(1/(gamma-1));
    M1 = abs(max(real(double(vpasolve(eqn, M)))));
end

function [M2] = mach_change(M1,gamma)
%   Normal shock realtion, solve downstream Mach by given upstream
    M2 = sqrt(((gamma-1).*(M1.^2) + 2)./(2.*gamma.*(M1^2) - (gamma - 1)));
end

function [M1] = mach_change_inverse(M2,gamma)
%   Normal Shock Relation, Solve upstream Mach by given downstream Mach
    syms M
    eqn = M2.^2 == ((gamma-1).*(M.^2) + 2)./(2.*gamma.*(M^2) - (gamma - 1));
    M1 = abs(min(real(double(vpasolve(eqn, M)))));
end

function [M2] = isentroPress(P02oP2, gamma)
%   Given P_stag/P_static, solve for flow mach number
    syms M
    eqn = P02oP2 == ((1 + ((gamma - 1)/2).*(M^2)).^(gamma./(gamma - 1)));
    M2 = abs(max(real(double(solve(eqn, M)))));
end


