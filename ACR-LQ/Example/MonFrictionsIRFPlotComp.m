% MonFrictionsIRFPlotComp
%
% Plots comparisons of the IRF in the Monetary Frictions model 
%
% Assumes that the LQ solutions were previously obtained and IRF generated
% and they are contained in a variables called 'IRF' with dimensions 
% ( (ny+ncsi) x nHorizon x ncsi )
%
% Required mat files: 
%   - LQDistortionsCashless
%   - LQNoDistortionsCashless       
%   - LQDistortionsFrictions
%   - LQNoDistortionsFrictions
%
% See also:
% LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar, LQAltRule, 
% MonFrictions, MonFrictionsIRFPlot, MonFrictionsIRFPlotAltRule,
% SmetsWouters
%
% .........................................................................
%
% Copyright 2004-2008 by Filipo Altissimo, Vasco Curdia and Diego Rodriguez 
% Palenzuela.
% Created: August 16, 2007
% Updated: January 29, 2008

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help MonFrictionsIRFPlotComp
% or
%    doc MonFrictionsIRFPlotComp

%% ------------------------------------------------------------------------

%% preamble
clear all
tic

%% designate and label the variables to plot and scale
var_plot = {'Y','n','Pi','R'};
var_label = {'Y','n','\Pi','R'};
scale = [1,1,4,4]; % annualize inflation and interest rates

%% set scale and label of shocks
shock_label = {'tau','mu','G','BARC','BARH','A'};
% change this if you want to use different shock sizes other than 1%
shock_size = ones(6,1);

%% set comparison to make
Comparison = {'DistortionsFrictions','DistortionsCashless'};
% Comparison = {'DistortionsFrictions','NoDistortionsFrictions'};
% Comparison = {'DistortionsCashless','NoDistortionsCashless'};
CompLegend = {'frictions','cashless'};
% CompLegend = {'with distortions','no distortions'};

%% Plot IRFs
load(['MonFrictions_LQ_',Comparison{1}],'y','IRF')
y = {y(:).name};
lineStyle = {'-b','--r'};
tid = 0:1:size(IRF,2)-1;
for j=1:6
    figure('Name',['Responses to a shock of ',num2str(shock_size(j)),'% in ',shock_label{j}])
    for jj=1:4
        [tf,var_pos] = ismember(var_plot{jj},y);
        subplot(2,2,jj)
        for jC=1:2
            load(['MonFrictions_LQ_',Comparison{jC}],'IRF')
            plot(tid,shock_size(j)*scale(jj)*IRF(var_pos,:,j),lineStyle{jC})
            hold on
        end
        xlim([0,tid(end)])
        set(gca,'XTick',0:4:tid(end))
        title(var_label{jj})
        plot(tid,zeros(size(tid)),'k:')
        hold off
    end
    legend(CompLegend{:})
    eval(['print -depsc2 MonFrictions_LQ_',Comparison{1},'_',Comparison{2},'_',shock_label{j},'.eps'])
end

%% ------------------------------------------------------------------------
