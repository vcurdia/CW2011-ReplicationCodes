function setpath

% setpath
%
% Set path to codes needed.
%
% Created: August 18, 2016 by Vasco Curdia
% Copyright 2016 by Vasco Curdia

pathBase = '../../';
pathList = {...
    'ACR-LQ',...
    'Sims-Gensys',...
    'Sims-Optimize',...
    };
for j=1:length(pathList)
    pathAdd{j} = [pathBase,pathList{j}];
end
addpath(pathAdd{:})
