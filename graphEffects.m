% function graphEffects(data,Levels,varargin)
%
% INPUTS:
%
% (1) data: an nsubj x ncond matrix of RM-data to plot
%
% (2) Levels: for a 3 x 2 x 2 design you'd input [3,2,2]
%
% (3) Varargin takes the form {FactorsToPlot,FactorName,DV} ...
%     (i)   FactorsToPlot  = e.g. [2,3] or [1,2,3]
%     (ii)  FactorName     = e.g. {'Expectation','Attention}
%     (iii) DV             = 'criterion'
%     (iv)  levelNames     = { {'low','med','high'} , {'full','diverted'} }
%
% by maxine

function graphEffects(data,Levels,varargin)

try;
    DV='DV';
    % formatting
    fs         = 14;
    lw         = 3;
    ticklength = [0 0];
    
    % Count the number of plots
    nPlots = numel(varargin);
    
    % Count number of factors & rm-conditions
    nFactors = length(Levels);
    nCells   = prod(Levels);
    
    % Compute indices
    F = computeIndices(Levels,nFactors,nCells);
    
    % Run through plots
    for iplot = 1:nPlots
        
        plotInfo      = varargin{iplot};
        factorsToPlot = plotInfo{1};
        factorName    = plotInfo{2};
        DV            = 'DV';%plotInfo{4};
        levelNames    = plotInfo{3};
        nWay          = length(factorsToPlot);
        
        % Prepare colours
        nfac = factorsToPlot(end);
        for i = 1:nfac; 
            x = [0,0,0];
            x(i) = 0.5;
            cols{i}=x;
        end
        
        figure;
        
        % Is this a 1-way 2-way or 3-way plot?
        switch nWay
            
            case 1
                X = [];
                for iLevel = 1:Levels(factorsToPlot)
                    idx = F{factorsToPlot}{iLevel};
                    X(:,iLevel) = nanmean(data(:,idx),2);
                end
                M = nanmean(X,1); SE = getWSSE(X);
                
                %barplot_grouped(M,SE,{},factorName,'inputmeans','colors',cols);
                errorbarplot(M,SE,levelNames,factorName,DV,cols);
                
            case 2
                for i1 = 1:Levels(factorsToPlot(1))
                    for i2 = 1:Levels(factorsToPlot(2))
                        idx = intersect(F{factorsToPlot(1)}{i1},F{factorsToPlot(2)}{i2});
                        X(:, i1,i2) = nanmean(data(:,idx),2);
                    end
                end
                M  = squeeze(nanmean(X,1));
                SE = getWSSE(X);
                errorbarplot(M,SE,levelNames,factorName,DV,cols);
                
            case 3
                for i1 = 1:Levels(factorsToPlot(1))
                    subplot(Levels(factorsToPlot(1)),1,i1);
                    X = [];
                    for i2 = 1:Levels(factorsToPlot(2))
                        for i3 = 1:Levels(factorsToPlot(3))
                            idx = intersect(F{factorsToPlot(2)}{i2},intersect(F{factorsToPlot(1)}{i1},F{factorsToPlot(3)}{i3}));
                            X(:, i2, i3) = nanmean(data(:,idx),2);
                        end
                    end
                    M  = squeeze(nanmean(X,1));
                    SE  = getWSSE(X);
                    errorbarplot(M,SE,levelNames,factorName,DV,cols);
                    title(factorName{1}{i1},'FontSize',fs);
                end
                
        end
        xlabel('condition','FontSize',fs);

        
    end
    
    
catch err;save graphEffectserror.mat; end
end

%% ---------------------------------------------------------------------%
function F = computeIndices(Levels,nFactors,nCells)

F = cell(nFactors,1);

% Do last factor first
F{nFactors} = cell(Levels(end),1);

for iFactor = nFactors:-1:1
    
    % make it the right size
    F{iFactor} = cell(Levels(1),1);
    cellLength = nCells/Levels(iFactor);
    numlev     = Levels(iFactor);
    
    % Is it the first level?
    if iFactor == 1
        for iLevel = 1:numlev
            F{iFactor}{iLevel} = [1 + cellLength*(iLevel-1):cellLength*iLevel];
        end
        
        % Is it the last level?
    elseif iFactor == nFactors
        for iLevel = 1:numlev
            F{iFactor}{iLevel} = [iLevel:numlev:nCells];
        end
        
        % If neither, we need to take the first 1/Levels(iFactor) from each of the factors in the level above and merge them.
    else
        lengthOfEachCell    = nCells/numlev; % each of the numlev cells will have this many elements
        nreps               = lengthOfEachCell/Levels(iFactor+1);
        
        for iLevel = 1:numlev
            X = [];
            for irep = 1:nreps
                basemult = (irep-1)*(Levels(iFactor)*Levels(iFactor+1));
                x = [1:Levels(iFactor+1)] + Levels(iFactor+1)*(iLevel-1) + basemult;
                X = [X,x];
            end
            F{iFactor}{iLevel} = X;
        end
    end
end
end

%% ----------------------------------------------------

function WSSE = getWSSE(x)
[nsubj,l1,l2] = size(x);
data          = reshape(x,nsubj,l1*l2);

[nsubj,ncond] = size(data);

% Get grand mean
grand_mean = nanmean(mean(data));

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
    WSSE(c)    = nanstd(data_adj(:,c))/sqrt(nsubj);
end

WSSE = reshape(WSSE,l1,l2);
end

%% -------------------------------------------------------------------%%
function errorbarplot(M,SE,levelNames,factorName,dv,cols);

try
    
    % check everything is the right way up
    if ismember(0,sum(size(M)==size(SE)));SE=SE';end
    
    % get plotting params
    lw  = 2.5;
    fs  = 12;
    s   = 100;
    
    % get sizes
    [nL1 nL2] = size(M);
    
    % Get x tick labels for...
    
    % a one factor plot
    if nL1 == 1 || nL2==1;
        xnames    = levelNames{1};
        xcondname = factorName{1};
        
    % a 2 factor plot ( 1 > 2)    
    elseif numel(levelNames{1}) >= numel(levelNames{2});
        xnames      = levelNames{1};
        legendNames = levelNames{2};
        xcondname   = factorName{1};
    
    % a 2 factor plot (2 > 1)    
    else xnames     = levelNames{2};
        legendNames = levelNames{1};
        xcondname   = factorName{2};
    end
    
    % ---------------------------------------------------------------------
    % If we're plotting F1 on the x-axis...
    % ---------------------------------------------------------------------
    if nL1 >= nL2

        % plot a line bar for each level of F2
        for i = 1:nL2
            l=errorbar(1:nL1,M(:,i),SE(:,i),SE(:,i));set(l,'color',cols{i},'LineWidth',lw);
            hold on;
        end
        
        % plot the legend
        if nL1>1&nL2>1
            hleg = legend(legendNames);
        end

        % add dots onto the lines
        for i = 1:nL2
            scatter(1:nL1,M(:,i),s,'MarkerFaceColor',cols{i},'MarkerEdgeColor','k');
            hold on;
        end

        % format
        set(gca,'XTickLabel',xnames,'XTick',1:nL1,'LineWidth',lw,'FontSize',fs,'XLim',[0,nL1+1],'TickLength',[0 0])
        set(gca,'XLim',[0 numel(xnames)+1]);
        hleg.String(end) = [];
        
    % ---------------------------------------------------------------------
    % If we're plotting F2 on the x-axis...
    % ---------------------------------------------------------------------    
    else
        
        % plot a line bar for each level of F1
        for i = 1:nL1
            l=errorbar(1:nL2,M(i,:),SE(i,:),SE(i,:));set(l,'color',cols{i},'LineWidth',lw);
            hold on;
        end

        % plot the legend
        if nL1>1&nL2>1
            hleg = legend(legendNames);
        end

         % add dots onto the lines
        for i = 1:nL1
            scatter(1:nL2,M(i,:),s,'MarkerFaceColor',cols{i},'MarkerEdgeColor','k');
            hold on;
        end

        % format
        set(gca,'XTickLabel',xnames,'XTick',1:nL2,'LineWidth',lw,'FontSize',fs,'XLim',[0,nL2+1],'TickLength',[0 0])
        set(gca,'XLim',[0 numel(xnames)+1]);
        hleg.String(end) = [];
    end
    
    hleg.String(end) = [];
    %hleg.String(end) = [];
    ylabel(dv,'FontSize',fs);
    xlabel(xcondname,'FontSize',fs);
    
    if nL1>1 & nL2>1
        title([factorName{1} ' x ' factorName{2}],'FontSize',fs);
    else title(['Main effect of ' factorName])
    end
    
    set(gca,'TickLength',[0,0],'XTickLabels',xnames);
    

catch err;save err_plot.mat; end;
end

