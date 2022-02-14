% function rm = rm_ANOVA( data , levels , [factorNames] , [levelNames] , [dvName] , [ploton] , [threshold] )
%
% This function:
% -computes an n-factor rmANOVA
% -displays F and p values for each main effect/interaction
% -outputs the full table
% -plots effects that are significant
% !!function can't handle more than 4-way!!
%
% Warning - if you want to use the plot function make sure you have
% graphEffects.m as well.
%
% INPUTS:
% For an m x n x p x ... design:
%
% data              = nsubj x prod(levels) matrix
% levels            = [m n p ...]                    [IN ORDER OF DATA MATRIX]
%
% OPTIONAL INPUTS:
%
% factorNames       = {'Attention','Expectation'}    [IN ORDER OF DATA MATRIX]
% levelNames        = { { 'Full' , 'Diverted' } , { 'E25' , 'E50' , 'E75' } } [IN ORDER OF DATA MATRIX]
% 
% [OPTIONAL] dvName = 'dprime' (defaults to 'DV')
% [OPTIONAL] ploton plots all sig. results (defaults to 0)
% [OPTIONAL] threshold = threshold for significance (defaults to 0.05)
%
%
%
% EXAMPLE:
%
% 5 subjects in an intervention (drug/placebo) x time (pre/post) design
%       drug_t1     drug_t2   placebo_t1    placebo_t2
%         3            1          3             2
%         5            2          6             6
%         4            3          4             5
%         7            2          5             4
%         2            3          2             3
%
% data        = [ 3 1 3 2; 5 2 6 6 ; 4 3 4 5; 7 2 5 4 ; 2 3 2 3];
% levels      = [2 2];
% factorNames = { 'intervention' , 'time' }
% levelNames  = { {'drug' , 'placebo'} , {'pre' , 'post'} }
% dvName      = 'depressionScore';
% ploton      = 1;
%


% maxine 11/8/17
% THIS IS A WRAP-AROUND

function [ranovatbl, X] = rm_ANOVA( data , levels , factorNames , levelNames , dvName , ploton , threshold );

%% how many sf to present?
nsf = 3;
nsf = ['%.' num2str(nsf) 'f'];

try;
  
    
    %% sort out arguments

    % if factor names aren't present, create them
    if nargin < 3 | isempty(factorNames)
        factorNames = cell(numel(levels),1);
        for iFactor = 1:numel(levels)
            factorNames{iFactor,1} = ['Factor' num2str(iFactor)];
        end
    end     

    % if level names aren't present, create them.
    if nargin < 4 | isempty(levelNames)
        levelNames = cell(numel(levels),1);
        for iFactor = 1:numel(levelNames)
            for jLevel = 1:levels(iFactor)
            levelNames{iFactor,1}{jLevel,1} = ['F' num2str(iFactor) 'L' num2str(jLevel)];
            end
        end
    end

    if nargin < 5 | isempty(dvName); dvName = 'DV'; end
    if nargin < 6 | isempty(ploton); ploton = 0; end
    if nargin < 7| isempty(threshold); threshold = 1; end
    if numel(levels) > 3; ploton = 0; end % can only plot up to 3.
    
    
%     %% Check levelNames are the right way around (should be row vectors)
%     for i = 1:numel(levelNames)
%         L = levelNames{i};
%         L = reshape(L,numel(L),1);
%         levelNames{i} = L;
%     end
    
    %% Check factorName are the right way around (should be column vectors)
    factorNames = reshape(factorNames,1,numel(factorNames));
    
    
    %% put our data into a table.
    for i = 1:size(data,2);
        varNames{i} = ['Y' num2str(i)];
    end
    t = array2table(data,'VariableNames',varNames);
    
        
    %% Create a table reflecting the within subject factors
    nFactors = numel(levels);
    
    switch nFactors
        
        case 1;
            within = table(levelNames,'VariableNames',factorNames);
            wModel = [factorNames{1}];
            
        case 2;
            L1 = levelNames{1}; L2 = levelNames{2};
            f1 = {};            f2 = {};
            
            for i = 1:levels(1)
                f1 = [f1 ; repmat(L1(i),levels(2),1)];
                f2 = [f2 ; L2];
            end
            within = table(f1 , f2 ,'VariableNames',factorNames);
            wModel = [factorNames{1} '*' factorNames{2}];
            
        case 3;
            L1 = levelNames{1}; L2 = levelNames{2}; L3 = levelNames{3};
            f1 = {};            f2 = {};            f3 = {};
            
            for i = 1:levels(1)
                for j = 1:levels(2)
                    f1 = [f1 ; repmat(L1(i),levels(3),1)];
                    f2 = [f2 ; repmat(L2(j),levels(3),1)];
                    f3 = [f3 ; L3];
                end
            end
            within = table(f1 , f2 , f3 , 'VariableNames',factorNames);
            wModel = [factorNames{1} '*' factorNames{2} '*' factorNames{3}];
            
        case 4;
            L1 = levelNames{1}; L2 = levelNames{2}; L3 = levelNames{3}; L4 = levelNames{4};
            f1 = {};            f2 = {};            f3 = {};            f4 = {};
            
            for i = 1:levels(1)
                for j = 1:levels(2)
                    for k = 1:levels(3)
                        f1 = [f1 ; repmat(L1(i),levels(4),1)];
                        f2 = [f2 ; repmat(L2(j),levels(4),1)];
                        f3 = [f3 ; repmat(L3(k),levels(4),1)];
                        f4 = [f4 ; L4];
                    end
                end
            end
            within = table(f1 , f2 , f3 , f4, 'VariableNames',factorNames);
            wModel = [factorNames{1} '*' factorNames{2} '*' factorNames{3} '*' factorNames{4}];
    end
    
    
    
    %% fit the repeated measures model
    modelSpec = ['Y1-Y' num2str(prod(levels)) '~1'];
    [rm] = fitrm(t,modelSpec,'WithinDesign',within);
    
    
    %% run repeated measures anova here
    ranovatbl = ranova(rm, 'WithinModel',wModel);
    
    %% display only relevant things
    X     = ranovatbl(3:2:end,4:5); % F and p value of effects
    pvals = table2array(X(:,2));
    disp(X( pvals < threshold , : ));

    %% We're only going to plot significant terms
    pvals  = table2array(X(:,2));
    plotMe = find(pvals < threshold);
    
    if ~isempty(plotMe) & ploton
        
        switch nFactors
            
            % 1 way ANOVA?
            case 1;
                graphEffects(data,levels,{1,{factorNames{1}}, {levelNames{1}}, dvName});
                title(factorNames{1});
                % 2-way ANOVA?
            case 2;
                
                % main effects?
                y = plotMe(plotMe <= nFactors);
                for i = 1:numel(y)
                    graphEffects(data,levels,{[y(i)],{factorNames{y(i)}}, {levelNames{y(i)}}, dvName});
                    title(factorNames{i});
                end
                % interaction?
                if ismember(3,plotMe)
                    graphEffects(data,levels,{[1 2],{factorNames{1} factorNames{2}},{levelNames{1} levelNames{2}}, dvName});
                    title('2-way')
                    xlabel(factorNames{2});
                end
                
                % 3 way ANOVA?
            case 3;
                
                % main effects?
                y=plotMe(plotMe < 4);
                for i = 1:numel(y)
                    graphEffects(data,levels,{[y(i)],{factorNames{y(i)}},{levelNames{y(i)}}, dvName});
                    %title(factorNames{i});
                end
                
                % 3-way interaction?
                if ismember(7,plotMe)
                    graphEffects(data,levels,{[1 2 3],{factorNames{1}, factorNames{2}, factorNames{3}},{levelNames{1}, levelNames{2}, levelNames{3}}, dvName});
                    %title('3-way');
                    %xlabel(factorNames{3});
                end
                
                % 2-way interaction?
                if ismember(4,plotMe)
                    graphEffects(data,levels,{[1 2],{factorNames{1}, factorNames{2}},{levelNames{1}, levelNames{2}}, dvName});
                    %title([factorNames{1} ' x ' factorNames{2}]);
                    %xlabel(factorNames{2});
                end
                if ismember(5,plotMe)
                    graphEffects(data,levels,{[1 3],{factorNames{1}, factorNames{3}},{levelNames{1}, levelNames{3}}, dvName});
                    %title([factorNames{1} ' x ' factorNames{3}]);
                    %xlabel(factorNames{3});                    
                end
                if ismember(6,plotMe)
                    graphEffects(data,levels,{[2 3],{factorNames{2}, factorNames{3}},{levelNames{2}, levelNames{3}}, dvName});
                    %title([factorNames{2} ' x ' factorNames{3}]);
                    %xlabel(factorNames{3});
                end
                
        end
    end


    %% report effects
    
    % ignore the 2 rows because it's the intercept
    for i_effect = 3:numel(ranovatbl.Properties.RowNames)
        
        % check it's an effect row
        if isempty(strfind(ranovatbl.Properties.RowNames{i_effect},'Error'))
            
            %% get the name
            name = ranovatbl.Properties.RowNames{i_effect}(13:end);
            
            %% get F
            F = sprintf(nsf,(ranovatbl.F(i_effect)));
            
            %% get DF
            DF = sprintf([nsf ',' nsf],[ranovatbl.DF(i_effect),ranovatbl.DF(i_effect+1)]);
            
            %% get P
            P  = sprintf('%.3f',ranovatbl.pValue(i_effect));
            
            %% get partial eta sq
            ETA = ranovatbl.SumSq(i_effect)/(ranovatbl.SumSq(i_effect)+ranovatbl.SumSq(i_effect+1));
            ETA = sprintf(nsf,ETA);
            
            %% write out effect
            disp(['<strong>' name ': </strong>']);
            disp(['F(' DF ') = ' F ', p = ' P ', etap2 = ' ETA ]);
        end
    end
            


catch err; save err_anova.mat;
end

end