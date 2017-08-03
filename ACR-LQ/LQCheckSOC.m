function SOC=LQCheckSOC(LQBETTA,LQS,LQmat,NumPrecision,MaxIterations)

% LQCheckSOC
%
% This function checks for optimal policy second order conditions.
% 
%   SOC=LQCheckSOC(LQBETTA,LQS,LQmat,NumPrecision,MaxIterations)
%
% Inputs:
%   - LQBETTA: time discount factor
%   - LQS: autocorrelation matrix for exogenous variables
%   - LQmat: structure with numerical matrices containing the LQ approx
%   - NumPrecision: (optional) numerical precision to evaluate the SOC. If
%       ommitted it is set to 1e-6.
%   - MaxIterations: (optional) maximum number of iterations for P11
%       convergence. (if MaxIterations is set, then NumPrecision also
%       needs to be set). Default value: 1000.
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
% Remark:
% In order to proceed the code needs to solve for matrix P11 (analyzed in
% Benigno and Woodford (2007). If it fails in solving for that matrix, the
% code will issue a warning message and switch to an alternative method
% (suggested in older versions of the method) and call LQCheckSOCAlt.
%
% Required Matlab routines:
%   - LQCheckSOCold (possibly)
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
%    help LQCheckSOC
% or
%    doc LQCheckSOC

%% ------------------------------------------------------------------------

%% Setup some background information

if nargin<4, NumPrecision = 1e-6; end
if nargin<5, MaxIterations = 1000; end

DispIterations = 0; % Show the evolution of iteration convergence in P11

ny = size(LQmat.A0,1);
nF = size(LQmat.C0,1);
nG = size(LQmat.D0,1);
ncsi = size(LQS,1);

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

%% Build support matrices
H = [zeros(ny,ncsi); LQmat.D2; LQS; eye(ncsi)];
% H_hat = [zeros(ny+nG+2*ncsi,ny+nG), H, zeros(ny+nG+2*ncsi,ncsi)];
G1 = [1/2*LQmat.A1;LQmat.C1;LQmat.D1];
G2 = [zeros(ny+nF,nG); -eye(nG)];
G4 = [LQmat.B2; zeros(nF+nG,ncsi)];

%% Solve for P11
M0 = [LQmat.A0, LQmat.C0', LQmat.D0';
      LQmat.C0, zeros(nF,nF+nG);
      LQmat.D0, zeros(nG,nF+nG)];
P11Conv = 0;
MP11 = zeros(ny+nF+nG);
P11old = eye(ny); %eye(ny); %zeros(ny);
itCount = 1;
while itCount<=MaxIterations
    MP11(1:ny,1:ny)=P11old;
    Mj = M0 + LQBETTA*MP11;
    P11 = -G1'*inv(Mj)*G1;
    P11Conv = abs(P11-P11old);
    P11ConvMax = max(P11Conv(:));
    if DispIterations
        disp(sprintf('Iteration %.0f - P11 convergence: %.10f',itCount,P11ConvMax))
    end
    P11Conv = all(P11Conv(:)<NumPrecision);
    if P11Conv
        break
    else
        P11old = P11;
        itCount = itCount+1;
    end
end
% Check solution
if ~P11Conv
    disp(' ')
    disp('WARNING: Solution for P11 did not converge!')
    disp(sprintf('         Iteration %.0f - P11 convergence: %.10f',itCount-1,P11ConvMax))
    disp(' ')
    disp('Cannot check for SOC with this method!')
    disp('Switching to alternative method...')
    SOC = LQCheckSOCAlt(LQBETTA,LQS,LQmat,NumPrecision);
    return
end

%% build matrix M
MP11 = zeros(ny+nF+nG);
MP11(1:ny,1:ny)=P11;
M = M0 + LQBETTA*MP11;
Minv = inv(M);

%% Solve for P12 and P14
P12 = -G1'*Minv*G2;
P14 = -G1'*Minv*G4;

%% Solve for P13
P13bar = -G1'*Minv*[LQmat.B0*LQS+LQmat.B1+LQBETTA*(P12*LQmat.D2+P14);-LQmat.C2;zeros(nG,ncsi)];
P13hat = -LQBETTA*G1'*Minv*[eye(ny);zeros(nF+nG,ny)];
P13 = reshape(inv(eye(ncsi*ny)-kron(LQS',P13hat))*P13bar(:),ny,ncsi);
G3 = [LQmat.B0*LQS+LQmat.B1+LQBETTA*(P12*LQmat.D2+P13*LQS+P14);-LQmat.C2;zeros(nG,ncsi)];

%% solve  for P22, P24, P44, P23, P34
P22 = -G2'*Minv*G2;
P24 = -G2'*Minv*G4;
P44 = -G4'*Minv*G4;
P23 = -G2'*Minv*G3;
P34 = -G3'*Minv*G4;

%% solve for P33
P33bar = -G3'*Minv*G3+LQBETTA*LQmat.D2'*(P22*LQmat.D2+P23*LQS+P24)...
         +LQBETTA*(P24'*LQmat.D2+P34'*LQS+P44)+LQBETTA*LQS'*(P23'*LQmat.D2+P34);
P33hat = LQBETTA^(1/2)*LQS';
P33 = reshape(inv(eye(ncsi^2)-kron(P33hat,P33hat))*P33bar(:),ncsi,ncsi);

%% Generate P
P = [ P11,  P12,  P13, P14;
     P12',  P22,  P23, P24;
     P13', P23',  P33, P34;
     P14', P24', P34', P44];

%% Create other support matrices
G = [G1, G2, G3, G4];
Phi1 = [-eye(ny), zeros(ny,nF+nG)]*Minv*G;
Phi11 = Phi1(:,1:ny);
Phi12 = Phi1(:,ny+1:ny+nG);
Phi13 = Phi1(:,ny+nG+1:ny+nG+ncsi);
Phibar = [Phi11, Phi1*H; zeros(ncsi,ny), LQS];
Psibar = [Phi13-Phi12*inv(P22)*P23;eye(ncsi)];

%% ------------------------------------------------------------------------

%% Check FOC
eigPhi11 = eig(Phi11);
FOCFail = sum((eig(Phi11)>=LQBETTA^(-1/2)));
if FOCFail>0
    error('Optimal solution is not determinate! Cannot check for SOC!')
end

fprintf('\nSecond Order Conditions Checks:')
fprintf('\n===============================\n')

%% Check SOC1
for j=(nF+nG+1):ny
    SOC1 = (sign(round(det(M((ny-j+1):end,(ny-j+1):end))/NumPrecision)*NumPrecision)==(-1)^j);
    if ~SOC1
        fprintf('Warning: First set of conditions failed.\n\n')
        break
    end
end

%% Check SOC2 (only if SOC1 verified)
SOC2 = 0;
if SOC1
    fprintf('First set of conditions passed.\n')
    for j=1:nG
        SOC2 = (sign(round(det(P22(1:j,1:j))/NumPrecision)*NumPrecision)==(-1)^j);
        if ~SOC2
            fprintf('Warning: Second set of conditions failed.\n\n')
        end
    end
    if SOC2
        fprintf('Second set of conditions passed.\n\n')
    end
end

%% ------------------------------------------------------------------------

%% return outputs
SOC = [SOC1; SOC2];

%% ------------------------------------------------------------------------
