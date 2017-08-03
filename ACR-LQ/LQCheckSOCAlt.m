function SOC=LQCheckSOCAlt(LQBETTA,LQS,LQmat,NumPrecision)

% LQCheckSOCAlt
%
% This function checks for optimal policy second order conditions using an
% alternative method
% 
%   SOC=LQCheckSOCAlt(LQBETTA,LQS,LQmat,NumPrecision)
%
% Inputs:
%   - LQBETTA: time discount factor
%   - LQS: autocorrelation matrix for exogenous variables
%   - LQmat: structure with numerical matrices containing the LQ approx
%   - NumPrecision: (optional) numerical precision to evaluate the SOC. If
%       ommitted it is set to 1e-6.
%
% Outputs:
%   - SOC: flag for SOC conditions (optional)
%       * first element refers to the first part of the conditions
%       * second element refers to the second part of the conditions
%
% The code will issue a warning message if the solution to the REE is not
% normal in any way.
%
% The code will check for second order conditions and issue a status
% message if any part of it fails.
%
% Required Matlab routines:
%   - gensys.m, by Chris Sims, available in his website
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
% Created: July 26, 2004
% Updated: March 27, 2009

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help LQCheckSOCold
% or
%    doc LQCheckSOCold

%% ------------------------------------------------------------------------

%% Setup some background information

if nargin<3, NumPrecision = 1e-6; end

ny = size(LQmat.A0,1);
nF = size(LQmat.C0,1);
nG = size(LQmat.D0,1);
ncsi = size(LQmat.B0,2);
nk = ny+ncsi+nF+nG+ny+ncsi+nG;

G0 = LQmat.G0;
G1 = LQmat.G1;
G2 = LQmat.G2;
G3 = LQmat.G3;

SOC = 0;

%% ------------------------------------------------------------------------

%% Check: S matrix
% the S matrix needs to have eigenvalues with modulus lower than
% beta^(-1/2)
if ~all(abs(eig(LQS))<LQBETTA^(-1/2))
    disp('WARNING: Eigenvalues of S need to have modulus lower than beta^(-1/2)!')
    disp('         SOC verification might be spurious.')
end

%% Check: rank of constraints
% rank([C0;D0]) = nF+ng
C0D0 = [LQmat.C0;LQmat.D0];
if rank(C0D0)<size(C0D0,1)
    disp('WARNING: Rank of [C0;D0] lower than number of constraints!')
    disp('         SOC verification might be spurious.')
end

%% ------------------------------------------------------------------------

%% Solve the REE using gensys
% % Check to see if there are any A0 rows corresponding to the FOC with
% % no expectations, due to the parameter values used (even if symbolically
% % they could be non-zero)
% cv = find(all(G0(1:ny,:)==0,2)~=0);
% G0(cv,:) = -G1(cv,:);
% G1(cv,:) = 0;
% G3(:,cv) = [];

% In this framework const is set to zero
const = zeros(nk,1);

[B1,Const,B2,fmat,fwt,ywt,gev,eu] = gensys(G0,G1,const,G2,G3);

if any(eu~=1)
    eu
    if all(eu)==-2
        warning('Coincident zeros. REE solution cannot be computed!!!')
    elseif eu(1)~=1
        warning('REE solution does not exist!!!')
    elseif eu(2)~=1
        warning('REE solution is not unique!!!')
    end
end

%% ------------------------------------------------------------------------

%% Build polynomial M(L)
LQM0 = zeros(nF+nG+ny); LQM1 = LQM0; LQM2 = LQM0; 
LQM0(nF+1:nF+nG,nF+nG+1:end) = LQBETTA^(-1/2)*LQmat.D0;
LQM0(nF+nG+1:end,1:nF) = LQBETTA^(1/2)*LQmat.C1';
LQM0(nF+nG+1:end,nF+nG+1:end) = 1/2*LQBETTA^(1/2)*LQmat.A1';
LQM1(1:nF,nF+nG+1:end) = LQmat.C0;
LQM1(nF+1:nF+nG,nF+nG+1:end) = LQmat.D1;
LQM1(nF+nG+1:end,1:nF) = LQmat.C0';
LQM1(nF+nG+1:end,nF+1:nF+nG) = LQmat.D1';
LQM1(nF+nG+1:end,nF+nG+1:end) = 1/2*(LQmat.A0+LQmat.A0');
LQM2(1:nF,nF+nG+1:end) = LQBETTA^(1/2)*LQmat.C1;
LQM2(nF+nG+1:end,nF+1:nF+nG) = LQBETTA^(-1/2)*LQmat.D0';
LQM2(nF+nG+1:end,nF+nG+1:end) = 1/2*LQBETTA^(1/2)*LQmat.A1;

%% check first part of SOC
LQSOC1 = 1;
for LQtheta=-pi:2*pi/100:pi
    Mbar = LQM0*exp(LQtheta*i)+LQM1+LQM2*exp(-LQtheta*i);
    for p=2*(nF+nG)+1:nF+nG+ny
        if sign(real(det(Mbar(1:p,1:p))))~=(-1)^(p-nF-nG)
            fprintf('principal minor %.0f: %.0g (should have sign %.0f)\n',...
                p,det(Mbar(1:p,1:p)),(-1)^(p-nF-nG));
            LQSOC1 = 0;
            break
        end
    end
    if ~LQSOC1, break, end
end

%% Check second part of SOC
fprintf('\nSecond Order Conditions Checks:')
fprintf('\n===============================\n')
if LQSOC1 && nG==0
    fprintf('First set of conditions passed.\n')
    fprintf('No need to evaluate second set of conditions (no forward looking constraints).\n')
elseif ~LQSOC1
    fprintf('Warning: First set of conditions failed.\n')
else
    % First part was successfull
    fprintf('First set of conditions passed.\n')
    
    % Check second part
    k_y   = 1:ny;
    k_csi = ny+1:ny+ncsi;
    k_GLM = ny+ncsi+nF+1:ny+ncsi+nF+nG;

    B1_yy   = B1(k_y,k_y);
    B1_yGLM = B1(k_y,k_GLM);
    B1_ycsi = B1(k_y,k_csi);
    B2_y    = B2(k_y,:);

    B1_GLMy   = B1(k_GLM,k_y);
    B1_GLMGLM = B1(k_GLM,k_GLM);
    B1_GLMcsi = B1(k_GLM,k_csi);
    B2_GLM    = B2(k_GLM,:);

    T0 = [B1_GLMGLM,B1_GLMy;B1_yGLM,B1_yy];
    Psi0 = [B2_GLM;B2_y];
    Psi1 = [B1_GLMcsi;B1_ycsi]-Psi0*LQS;
    T = LQBETTA^(1/2)*T0;

    % form J matrix:
    %   J = T'*X*T + T'*J*T
    % with
    %   X = S1'*(AL0+AL0')*S1 + beta^(1/2)*T'*S1'*AL1*S1 + beta^(1/2)*S1'*AL1'*S1*T
    %   S1 = [zeros(Ny,NG) eye(Ny)]
    S1 = [zeros(ny,nG) eye(ny)];
    X = S1'*(LQmat.A0+LQmat.A0')*S1 + LQBETTA^(1/2)*T'*S1'*LQmat.A1*S1 +...
        LQBETTA^(1/2)*S1'*LQmat.A1'*S1*T;
    TXT = T'*X*T;
%     J = inv(eye((ny+nG)^2)-kron(T',T'))*TXT(:);
%     J = reshape(J,ny+nG,ny+nG);
    J = lyapcsd(T',TXT);
    J11 = real(J(1:nG,1:nG));

    LQSOC2 = 1;
    dd2 = [];
    for p = 1:nG
        Jpp = J11(1:p,1:p);
        dd2(p,1) = det(Jpp);
        if round(det(Jpp)/NumPrecision)*NumPrecision==0
            LQSOC2 = 0;
        elseif sign(det(Jpp))~=(-1)^p
            LQSOC2 = 0;
        end
    end

    if LQSOC2
        fprintf('Second set of conditions passed.\n')
%             disp('list of 2nd part of SOC conditions (dd2):'), disp(dd2)
    else
        fprintf('Warning: Second set of conditions failed.\n')
        [V,D] = eig(J11);
        LQmu = LQmat.D0*S1*T*[eye(nG);zeros(ny,nG)]*V;
        idx = find(round(diag(D)/NumPrecision)*NumPrecision==0);
        DevComm = ~all(round(LQmu(:,idx)/NumPrecision)*NumPrecision==0);
        if ~DevComm
            fprintf('There are no optimal deviations from commitment.\n') 
        else
            fprintf('Optimal deviations from commitment might exist. Check results below.\n')
            fprintf('\nList of 2nd part of SOC conditions (dd2):\n'), disp(dd2)
            disp('If failure of SOC is due to semi-definite J11 then it is possible')
            disp('to check whether we have a problem or not. In particular check')
            disp('that for the eigenvalues that are zero we get actual deviations')
            disp('from the commitment (LQmu different from zero)')
            disp(' ')
            disp('eigenvalues of J11 (D):'), disp(D)
            disp('eigenvectors of J11 (V):'), disp(V)
            disp('deviations from commitment, mu vector (LQmu):')
            disp('(only the ones corresponding to zero eigenvalues are relevant)')
            disp(LQmu)
        end
    end
    disp(' ')
end

%% ------------------------------------------------------------------------
