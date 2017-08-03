% LQGenSymVar
%
% Creates symbolic variables for y, csi, epsilon, GLM and FLM, which are needed to
% run the LQ codes
%
% this script also counts variables and generates
%   ny:      number of y variables
%   ncsi:    number of csi variables
%   nG:      number of G constraints (and therefore the number of lagrange
%            multipliers associated to G cocnstraints)
%   nF:      number of F constraints (and therefore the number of lagrange
%            multipliers associated to F cocnstraints)
%
% Requires cell arrays:
%   y:      with 2 columns, one variable per row, first column contains 
%           strings with name of variables and second column contains the 
%           guess values for the steady state
%   csi:    simple cell array with strings corresponding to the various
%           exogenous variables
%   GLM:    with 2 columns, one lagrange multiplier per row, first column 
%           contains strings with name of variables and second column 
%           contains the guess values for the steady state
%   FLM:    with 2 columns, one lagrange multiplier per row, first column 
%           contains strings with name of variables and second column 
%           contains the guess values for the steady state
%
% See also:
% LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar, LQAltRule, 
% LQWEval, MonFrictions, MonFrictionsIRFPlot, MonFrictionsIRFPlotComp, 
% MonFrictionsIRFPlotAltRule, SmetsWouters
%   
% .........................................................................
%
% Copyright 2004-2009 by Filipo Altissimo, Vasco Curdia and Diego Rodriguez 
% Palenzuela.
% Created: August 16, 2007
% Updated: September 18, 2009

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help LQGenSymVar
% or
%    doc LQGenSymVar

%% ------------------------------------------------------------------------

%% endogenous variables
y = cell2struct(y,{'name','guess'},2);
ny = length(y);
y_ss = []; y_t = []; y_tF = []; y_tL = [];
for j=1:ny
    eval(['syms ',y(j).name,'_ss ',y(j).name,'_t ',y(j).name,'_tF ',y(j).name,'_tL'])
    y_ss  = [y_ss,  eval([y(j).name,'_ss'])];
    y_t   = [y_t,   eval([y(j).name,'_t'])];
    y_tF  = [y_tF,  eval([y(j).name,'_tF'])];
    y_tL = [y_tL, eval([y(j).name,'_tL'])];
end

%% exogenous variables
ncsi = length(csi);
csi_t = []; csi_tF = []; csi_tL = []; eps_t = [];
for j=1:ncsi
    eval(['syms ',csi{j},'_t ',csi{j},'_tF ',csi{j},'_tL'])
    csi_t  = [csi_t,  eval([csi{j},'_t'])];
    csi_tF = [csi_tF, eval([csi{j},'_tF'])];
    csi_tL = [csi_tL, eval([csi{j},'_tL'])];
    eval(['syms eps_',csi{j},'_t'])
    eps_t = [eps_t, eval(['eps_',csi{j},'_t'])];
end
csi_ss = sym(zeros(1,ncsi));

%% lagrange multipliers for F constraints
FLM = cell2struct(FLM,{'name','guess'},2);
nF = length(FLM);
FLM_ss = []; FLM_t = []; FLM_tF = [];
for j=1:nF
    eval(['syms ',FLM(j).name,'_ss ',FLM(j).name,'_t ',FLM(j).name,'_tF'])
    FLM_ss = [FLM_ss, eval([FLM(j).name,'_ss'])];
    FLM_t  = [FLM_t,  eval([FLM(j).name,'_t'])];
    FLM_tF = [FLM_tF, eval([FLM(j).name,'_tF'])];
end

%% lagrange multipliers for G constraints
GLM = cell2struct(GLM,{'name','guess'},2);
nG = length(GLM);
GLM_ss = []; GLM_t = []; GLM_tF = []; GLM_tL = [];
for j=1:nG
    eval(['syms ',GLM(j).name,'_ss ',GLM(j).name,'_t ',...
        GLM(j).name,'_tF ',GLM(j).name,'_tL'])
    GLM_ss  = [GLM_ss,  eval([GLM(j).name,'_ss'])];
    GLM_t   = [GLM_t,   eval([GLM(j).name,'_t'])];
    GLM_tF  = [GLM_tF,  eval([GLM(j).name,'_tF'])];
    GLM_tL = [GLM_tL, eval([GLM(j).name,'_tL'])];
end

%% ------------------------------------------------------------------------

