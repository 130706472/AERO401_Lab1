clc;clear;close all;

% DATA = readtable('SSWT_DATA_2021.xlsx');
DATA = load('DATA.mat').DATA;

t = table2array(DATA(:,1));

% Pitot tube pressure [psig]
P_pitot = table2array(DATA(:,2));

% Plenum pressure [psig] -> stagnation pressure
P_plen = table2array(DATA(:,2));