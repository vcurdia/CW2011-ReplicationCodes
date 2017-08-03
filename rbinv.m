function Ainv=rbinv(A,crit)

% function that computes robust inverses, using svd
%
% ..............................................................................
%
% Created: October 29, 2002 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% 
% Copyright 2002-2011 by Vasco Curdia

if nargin==1,crit=1e-12;end
[u,d,v]=svd(A);
d=diag(d);
first0=find(d<crit,1);
if isempty(first0),first0=size(A,1)+1;end
u=u(:,1:first0-1);
v=v(:,1:first0-1);
d=d(1:first0-1);
Ainv=v*diag(1./d)*u';