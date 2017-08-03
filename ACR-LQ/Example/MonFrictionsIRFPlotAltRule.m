% MonFrictionsIRFPlotAltRule
%
% Plots the IRF in the Monetary Frictions model for both the optimal policy
% rule and the alternative rule
%
% Required mat file: name of mat file is defined by user below
% Assumes that the LQ solutions were previously obtained and IRF generated
% and they are contained in a variables called 'IRF' with dimensions 
% ( (ny+ncsi) x nHorizon x ncsi )
%
% See also:
% LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar, LQAltRule, 
% MonFrictions, MonFrictionsIRFPlot, MonFrictionsIRFPlotComp, SmetsWouters
%
% .........................................................................
%
% Copyright 2004-2008 by Filipo Altissimo, Vasco Curdia and Diego Rodriguez 
% Palenzuela.
% Created: January 21, 2008
% Updated: January 29, 2008

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help MonFrictionsIRFPlotAltRule
% or
%    doc MonFrictionsIRFPlotAltRule

%% ------------------------------------------------------------------------

%% preamble
clear all
tic

%% designate mat file
MatFileName = 'MonFrictions_LQ_DistortionsCashless';
load(MatFileName,'y','IRF','IRF_Rule')
y = {y(:).name};

%% designate and label the variables to plot and scale
var_plot = {'Y','n','Pi','R'};
var_label = {'Y','n','\Pi','R'};
scale = [1,1,4,4]; % annualize inflation and interest rates

%% set scale and label of shocks
shock_label = {'tau','mu','G','BARC','BARH','A'};
ncsi = length(shock_label);
% change the following if you want to use different shock sizes other than 1%
shock_size = ones(ncsi,1);

%% Plot IRFs
tid = 0:1:size(IRF,2)-1;
for j=1:ncsi
    figure('Name',['Responses to a shock of ',num2str(shock_size(j)),'% in ',shock_label{j}])
    for jj=1:4
        [tf,var_pos] = ismember(var_plot{jj},y);
        subplot(2,2,jj)
        plot(tid,shock_size(j)*scale(jj)*IRF(var_pos,:,j),'-b')
        hold on
        plot(tid,shock_size(j)*scale(jj)*IRF_Rule(var_pos,:,j),'--r')
        title(var_label{jj})
        plot(tid,zeros(size(tid)),'k:')
        hold off
    end
    legend('Optimal','Rule')
    % can change, delete or comment the folloing command. Is simply
    % converts plot to an eps file.
    eval(['print -depsc2 IRF_',MatFileName,'_AltRule_',shock_label{j},'.eps'])
end

%% ------------------------------------------------------------------------

