% MakeFig5
%
% Generates Figure 5 of the paper.
%
% See also:
% Exercise2, Exercise3
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
Exercise2('dSP',12,'PersSP',90,'isIgnoreZLB',0)
Exercise2('dSP',12,'PersSP',90,'isIgnoreZLB',1)
Exercise3('dSP',12,'PersSP',90,'isIgnoreZLB',0)
Exercise3('dSP',12,'PersSP',90,'isIgnoreZLB',1)

%% Generate Plot
IRFPlotCompareExercise3(...
    'FileNameSuffix','_dSP_12_Pers_90',...
    'FigPrint',1,...
    'FigPrefix','Fig_5_')

%% Elapsed time
disp(' '), vctoc(ttic), disp(' ')

%% -----------------------------------------------------------------------------
