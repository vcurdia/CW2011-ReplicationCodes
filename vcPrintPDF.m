function []=vcPrintPDF(FileName,KeepEPS,OpenPDF)

% vcPrintPDF
%
% Prints current figure into eps and then convert it to pdf and deletes the
% original eps file.
%
% Arguments:
%   - FileName: (string) the name of the file (without extension)
%   - KeepEPS: (optional)(logical) Keep both the pdf and eps files? (if
%              this one is specified then ColorFigure needs to be specified
%              too)
%
% Requirements:
%   - epstopdf
%
% ..............................................................................
%
% Created: August 9, 2007 by Vasco Curdia
% Updated: July 26, 2011 by Vasco Curdia
% Updated: October 23, 2013 by Vasco Curdia
% 
% Copyright 2007-2013 by Vasco Curdia

% -------------------------------------------------------------------------

%% Check Arguments
if nargin==1 || isempty(KeepEPS)
  KeepEPS = 0;
end
if nargin<3
  OpenPDF = 0;
end

%% print eps file
print('-depsc2',[FileName,'.eps'])

%% convert to pdf
eval(['!epstopdf ',FileName,'.eps'])

%% delete eps file
if ~KeepEPS
  eval(sprintf('delete %s.eps', FileName))
end

%% open pdf
if OpenPDF
  open([FileName,'.pdf'])
end
    