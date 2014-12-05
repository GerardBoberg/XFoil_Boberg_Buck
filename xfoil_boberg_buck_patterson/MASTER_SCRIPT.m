
%% xfoil project Aero 306 -- MASTER_SCRIPT
% Gerard Boberg, Trevor Buck, Zane Patterson
%
% 4 Dec 2014
%
% This script is equivelent to running both the "RunSingleAnalysis" and the
% "SweepAngleAttack" scripts in short succession.
%    

clc
clear all
close all

% Make the text output look nice

disp( ' ' )
disp( ' ' )
disp('   Running Single Analysis')
disp('------------------------------')
disp( ' ' )
run RunSingleAnalysis  % Actually runs the script. Go into the script to
disp( ' ' )            %    modify parameters
disp( ' ' )
disp('------------------------------')
disp( ' ' )
disp( ' ' )
disp('      Running AoA Sweep')
disp('------------------------------')
disp( ' ' )
run SweepAngleAttack % Actually runs the script. Go into the script to
disp( ' ' )          %      modify parameters
disp( ' ' )
disp('------------------------------')
disp( ' ' )
disp( ' ' )

% End of File