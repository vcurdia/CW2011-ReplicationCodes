function out=vctoc(ttic,ttoc)

% vctoc
%
% This file analyzes and reports elapsed time in more detailed fashion, as
% compared to simple toc command.
% NOTE: it assumes that toc command was used before.
%
% Usage:
%
%   vctoc
%       uses global time variable toc as reference
%
%   vctoc(ttic)
%       uses ttic as the start point, where ttic is toc value saved earlier
%
%   vctoc(ttic,ttoc)
%       uses ttic as the start point, and ttoc for end point. If ttic=[]
%       then it will be the elapsed time as measured in ttoc
%
%   out=vctoc
%   out=vctoc(ttic)
%   out=vctoc(ttic,ttoc)
%       same as above but saves txt in out, instead of displaying it
%
% .........................................................................
%
% Created: October 19, 2002
% Updated: November 10, 2008
% by Vasco Curdia

% -------------------------------------------------------------------------

if nargin==1
    sElapsed = toc-ttic;
elseif nargin==2
    if isempty(ttic)
        sElapsed = ttoc;
    else
        sElapsed = ttoc-ttic;
    end
else
    sElapsed = toc;
end

hElapsed = fix(sElapsed/(60^2));
sElapsed = sElapsed - hElapsed*(60^2);
mElapsed = fix(sElapsed/60);
sElapsed = sElapsed - mElapsed*60;
rElapsed = sElapsed - fix(sElapsed);
sElapsed = sElapsed - rElapsed;
rElapsed = fix(100*rElapsed);

txt = 'Elapsed time:';
if hElapsed>0, txt = sprintf('%s %.0fh',txt,hElapsed); end
if (mElapsed>0)||(hElapsed>0), txt = sprintf('%s %.0fm',txt,mElapsed); end
txt = sprintf('%s %.0fs %02.0f',txt,sElapsed,rElapsed);

if nargout==1
    out = txt;
else
    disp(txt)
end
