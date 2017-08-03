function W=LQWEval(LQmat,REELQ,REERule,VarShocks)

% LQWEval
%
% This function uses the LQ approximation method proposed by Benigno and
% Woodford (2008) to evaluate the welfare of a given policy.
%
% NOTE: solution assumes stationary solution under given rule.
% 
% Usage:
%
%   W=LQWEval(LQmat,REELQ,REERule,VarShocks)
%
% Inputs:
%   - LQmat: structure with numerical matrices for the LQ approx
%   - REELQ: structure with REE solution under LQ
%   - REERule: structure with REE solution under rule (if welfare of LQ
%       desired, simply use REERule=REELQ
%   - VarEps: covariance matrix of shocks
%
% Outputs:
%   - LQW: the value of welfare of a given rule
%
% Required Matlab routines:
%   - Symbolic Toolbox
%   - LQ toolbox
%   - lyapcsd.m and lyapcs.m, by Chris Sims, available in his website
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
% Created: March 24, 2009
% Updated: March 27, 2009

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help LQWEval
% or
%    doc LQWEval

%% ------------------------------------------------------------------------

%% check gensys eu
if ~all(REERule.eu==1), W = -inf; return, end

%% normalize system to extract Lagrange multipliers
ny = size(LQmat.A0,1);
ncsi = size(LQmat.B0,2);
nF = size(LQmat.C0,1);
nG = size(LQmat.D0,1);
nz = length(REERule.z_t);
isvphi = isfield(REERule,'Svphi');
if isvphi
    Phi1 = REERule.Phi1;
    Phi2 = REERule.Phi2;
    Sy = REERule.Sy;
    Scsi = REERule.Scsi;
    Svphi = REERule.Svphi;
else
    Phi1 = [REERule.Phi1,zeros(nz,nG);REELQ.Phi_vphi_z,zeros(nG,nz-ny-ncsi),REELQ.Phi_vphi_vphi];
    Phi2 = [REERule.Phi2;REELQ.Phi_vphi_eps];
    Sy = [REERule.Sy,zeros(ny,nG)];
    Scsi = [REERule.Scsi,zeros(ncsi,nG)];
    Svphi = [zeros(nG,nz),eye(nG)];
end

%% basic terms
W1 = Sy'*LQmat.A0*Sy+2*Sy'*(LQmat.B0*LQmat.S+LQmat.B1)*Scsi;
W2 = Sy'*LQmat.A1*Sy+2*Sy'*LQmat.B2*Scsi;
W3 = Svphi'*LQmat.D0*Sy;

%% Solve for covariance
V = real(lyapcsd(Phi1,Phi2*VarShocks*Phi2'));

%% Solve for expectational terms
E_1_P_y_plus = real(lyapcsd(LQmat.beta^(1/2)*Phi1,Phi2*VarShocks*Phi2'));
E_mu_P_y_cyc = real(lyapcsd(LQmat.beta^(1/2)*Phi1,(1-LQmat.beta)*V));

%% combine welfare terms
W = 1/2/(1-LQmat.beta)*trace(W1*(E_1_P_y_plus+Phi1*E_mu_P_y_cyc))+...
    1/2/(1-LQmat.beta)*trace(W2*(LQmat.beta*E_1_P_y_plus*Phi1'+E_mu_P_y_cyc*Phi1'))+...
    1/LQmat.beta*trace(W3*Phi1*V);

%% ------------------------------------------------------------------------

