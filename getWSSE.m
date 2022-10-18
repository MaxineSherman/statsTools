% function WSSE = getWSSE(data)
% 
% This computes the within-subject standard error. Use this for repeated-measures designs
% 
% INPUTS:
% data                   An nsubjects x nconditions matrix. All conditions should be repeated-measures
%
% OUTPUTS:
% WSSE                   A 1 x nconditions vector containing the within-subject standard error for each condition
%
%
% =========================================================================
% v1.0 
% For questions/comments please email Maxine Sherman
% m.sherman@sussex.ac.uk / maxinesherman@gmail.com
% =========================================================================
                

function WSSE = getWSSE(data)


[nsubj,ncond] = size(data);

% Get grand mean
grand_mean = nanmean(nanmean(data)); 

% Get participant means
p_mean = nanmean(data,2);

% Get adjustment factors
adj_factor =  - p_mean + grand_mean;

% Adjust data
data_adj = [];
for p = 1:nsubj
    data_adj(p,:) = data(p,:) + adj_factor(p);
end

% get SE
WSSE=[];
for c = 1:ncond
    WSSE(1,c)    = nanstd(data_adj(:,c))/sqrt(nsubj);
end

