% MakeFig9
%
% Generates Figure 9 of the paper.
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
Exercise4('dSP',4,'PersSP',90,'BindBp',10)
Exercise4SearchSequence(...
    'FileNameSuffix','_dSP_4_Pers_90_Bind_10bp',...
    'Shocks2Plot',{'hchitil'})

%% Generate Plot
IRFPlotCompareExercise4(...
    'FileNameSuffix','_dSP_4_Pers_90_Bind_10bp',...
    'Shocks2Plot',{'hchitil'},...
    'FigShape',{3,3},...
    'FigPrint',1,...
    'FigPrefix','Fig_9_')

%% Elapsed time
disp(' '), vctoc(ttic), disp(' ')

%% -----------------------------------------------------------------------------
