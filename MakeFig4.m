% MakeFig4
%
% Generates Figure 4 of the paper.
%
% See also:
% Exercise2, IRFPlotCompareExercise2
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

%% Generate simulations
Exercise2('dSP',4,'PersSP',90)

%% Generate plot
IRFPlotCompareExercise2(...
    'FileNameSuffix','_dSP_4_Pers_90',...
    'FigPrint',1,...
    'FigPrefix','Fig_4_')

%% Elapsed time
disp(' '), vctoc(ttic), disp(' ')

%% -----------------------------------------------------------------------------
