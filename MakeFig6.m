% MakeFig6
%
% Generates Figure 6 of the paper.
%
% See also:
% Exercise5, Exercise5SearchSequence
%
% ..............................................................................
%
% Created: April 21, 2014 by Vasco Curdia
%
% Copyright 2014 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Preamble
clear all
tic
ttic = toc();

%% Generate Simulations
Exercise5('dSP',4,'PersSP',90)
Exercise5SearchSequence(...
    'FileNameSuffix','_dSP_4_Pers_90',...
    'Shocks2Plot',{'hXitil'})

%% Generate Plot
IRFPlotCompareExercise5(...
    'FileNameSuffix','_dSP_4_Pers_90',...
    'Shocks2Plot',{'hXitil'},...
    'FigShape',{3,3},...
    'FigPrint',1,...
    'FigPrefix','Fig_6_')

%% Elapsed time
disp(' '), vctoc(ttic), disp(' ')

%% -----------------------------------------------------------------------------
