% MakeFig7
%
% Generates Figure 7 of the paper.
%
% See also:
% Exercise4, Exercise4SearchSequence
%
% ..............................................................................
%
% Created: April 22, 2014 by Vasco Curdia
%
% Copyright 2014 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Preamble
clear all
tic
ttic = toc();

%% Generate Simulations
Exercise4('dSP',4,'PersSP',90)
Exercise4SearchSequence(...
    'FileNameSuffix','_dSP_4_Pers_90',...
    'Shocks2Plot',{'hXitil'})

%% Generate Plot
IRFPlotCompareExercise4(...
    'FileNameSuffix','_dSP_4_Pers_90',...
    'Shocks2Plot',{'hXitil'},...
    'FigShape',{3,3},...
    'FigPrint',1,...
    'FigPrefix','Fig_7_')

%% Elapsed time
disp(' '), vctoc(ttic), disp(' ')

%% -----------------------------------------------------------------------------
