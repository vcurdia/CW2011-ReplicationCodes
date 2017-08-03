function REE=LQAltRuleLinear(Rule,LQmat,LQS,LQz_t,NewVars)

% LQAltRule
%
% Generates state space and REE solution for alternative rule that is linear.
%
%   Out=LQAltRuleLinear(Rule,LQmat,LQS,LQz_t,NewVars)
%
% Inputs:
%   - Rule: symbolic vector with equations defining the rule(s)
%     [See below remarks about rule]
%   - LQmat: structure with numerical matrices for the LQ approx
%   - LQS: numeric autocorrelation matrix for exogenous variables
%   - LQz_t: symbolic vector with symbolic variables of state space
%   - NewVars: symbolic vector with additional symbolic variables (dated t)
%
% Outputs:
%   - REE: structure with the following fields:
%     * REE.Phi1 and Out.Phi2: matrices for the REE solution such that:
%       z_t = Phi1*z_tL + Phi2*eps_t
%     * REE.z_t: symbolic vector with symbolic variables used in the state
%       space characterizing the rule.
%
% The code will issue a warning message if the solution to the REE is not
% normal in any way.
%
% Remarks about rule(s):
% ----------------------
%
% This code is designed to receive as an input linear equations describing
% a rule, based on the initial variables introduced in LQ and possibly new 
% ones.
% 
% General form of rule:
%   Rule(z_t,z_tF,z_tL)=0
%
% Each equation can contain forward looking and backward looking terms, 
% but not both at the same time. If both are desired then artificial 
% variables should be created, and additional equations introduced (for 
% this example: MonFrictions). 
%
% NOTE: This rule should be formulated so that it is consistent with the 
%       optimal steady state 
%
% Convention for refering to variables in equations:
%   x_t   refers to x{t}
%   x_tF  refers to x{t+1}
%   x_tL refers to x{t-1}
%   x_ss  refers to the steady state level of 'x'
%   (where x refers to some variable with name 'x')
%
% Required m-files:
%   - symbolic toolbox
%   - gensys, available in Chris Sims's website
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
% Created: May 1, 2008
% Updated: September 18, 2009

%% ------------------------------------------------------------------------

%% Define some variables
nR = length(Rule);
ny = size(LQmat.A0,1);
nF = size(LQmat.C0,1);
nG = size(LQmat.D0,1);
ncsi = size(LQmat.B0,2);

%% Analyze variables in the Rule and augment LQz_t (if necessary)
RLz_t = LQz_t(1:ny+ncsi);
if exist('NewVars','var')
    RLz_t = [RLz_t,NewVars];
end
nz = length(RLz_t);

%% check number of rules
if nR~=(nz-nF-nG-ncsi)
    error('Number of rules is %.0f, but number of rules needed is %.0f!',nR,nz-nF-nG)
end

%% ------------------------------------------------------------------------

%% Define z_tF and z_tL in case they are needed
RLz_tF = RLz_t; RLz_tL = RLz_t;
for j=1:nz
    RLz_tF(j) = subs(RLz_t(j),RLz_t(j),[char(RLz_t(j)),'F']);
    RLz_tL(j) = subs(RLz_t(j),RLz_t(j),[char(RLz_t(j)),'L']);
end

%% ------------------------------------------------------------------------

%% generate state space with rule

% introduce the constraints
%   0 = LQmat.D0*E_t[hy_tF]+LQmat.D1*hy_t-LQmat.D2*csi_t
%   0 = LQmat.C0*hy_t+LQmat.C1*hy_tL-LQmat.C2*csi_t
%   0 = LQS*csi_tL+eps_t-csi_t
G0 = [-LQmat.D0 zeros(nG,nz-ny);
      -LQmat.C0 LQmat.C2 zeros(nF,nz-ny-ncsi);
      zeros(ncsi,ny) eye(ncsi) zeros(ncsi,nz-ny-ncsi)];
G1 = [LQmat.D1 -LQmat.D2 zeros(nG,nz-ny-ncsi);
      LQmat.C1 zeros(nF,nz-ny);
      zeros(ncsi,ny) LQS zeros(ncsi,nz-ny-ncsi)];
G2 = [zeros(nG+nF,ncsi);
      eye(ncsi);
      zeros(nR,ncsi)];
G3 = [eye(nG);
      zeros(nF+ncsi,nG)];

% in this framework the constant is zero
const = zeros(nz,1);

for jr=1:nR
    % check rule type
    gtt = double(jacobian(Rule(jr),RLz_tF));
    gt = double(jacobian(Rule(jr),RLz_t));
    gt1 = double(jacobian(Rule(jr),RLz_tL));
    if ~all(gtt==0)
        if ~all(gt1==0)
            error('Cannot have both forward and backward looking variables in same rule!')
        end
        G0 = [G0; -gtt];
        G1 = [G1; gt];
        G3 = [G3, zeros(size(G3,1),1); zeros(1,size(G3,2)) 1];
    else
        G0 = [G0; -gt];
        G1 = [G1; gt1];
        G3 = [G3; zeros(1,size(G3,2))];        
    end
end

%% ------------------------------------------------------------------------

%% Solve for REE solution 
[Phi1,Const,Phi2,fmat,fwt,ywt,gev,eu] = gensys(G0,G1,const,G2,G3);
if any(eu~=1)
    fprintf('Warning: eu = (%.0f,%.0f)',eu)
    if all(eu)==-2
        fprintf(' - Coincident zeros!!!')
    elseif eu(1)~=1
        fprintf(' - Solution does not exist!!!')
    elseif eu(2)~=1
        fprintf(' - Solution is not unique!!!')
    end
    fprintf('\n')
end

%% ------------------------------------------------------------------------

%% Prepare output
REE.Phi1 = Phi1;
REE.Phi2 = Phi2;
REE.z_t = RLz_t;
REE.eu = eu;
REE.Sy = [eye(ny),zeros(ny,nz-ny)];
REE.Scsi = [zeros(ncsi,ny),eye(ncsi),zeros(ncsi,nz-ny-ncsi)];

%% ------------------------------------------------------------------------
