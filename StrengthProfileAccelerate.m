%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rheol.m
% 
% Laurent Montesi, 11/17/2014
% Tools for manipulating strength profiles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load planet data
load planet
% load rock data
load rock
% Some constant
Celsius=273.15;


% Define model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choose_planet
choose_th
choose_def
choose_model

calc_strength
ifig=1; plot_strength

Stotal=IntegrateStrength(model,e,thid,Temperature,rock,did);
% Second model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_initial=model;
clear model;
choose_model
efinal=10^(bisection(@(le)IntegrateStrength(model,10^le,thid,Temperature,rock,did)-Stotal,e));
e=efinal;
calc_strength
ifig=2; plot_strength

