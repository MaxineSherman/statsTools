% function [BF_empirical RR] = robustnessRegion( [shape_of_prior] , [data] , [prior_params], [sensitive_BF], [plotit], [timeout] )
%
% This function computes a bayes factor and robustness region for your
% data. The robustness region is estimated via brute search, and times out
% after 4 seconds. If the BF hasn't crossed the limits 1/3 or 3 after this
% time, then it's assumed that the robustness region extends to infinity.
%
% The code also plots the robustness region for you.
%
% This code is based around the BF Matlab code from Zoltan's website:
% http://www.lifesci.sussex.ac.uk/home/Zoltan_Dienes/inference/Bayes.htm
%
% WARNING: DO NOT USE A GAUSSIAN PRIOR WITH NON-ZERO MEAN. THIS PART OF THE CODE HAS NOT
% BEEN IMPLEMENTED PROPERLY.
%
% OUTPUTS:
%
% BF_empirical:   your bayes factor
% RR:             your robustness region.
%
% INPUTS:
%
% shape_of_prior (string): Inputs are 'normal', 'uniform' or 'minimal'
%                          Here, minimal means you're testing a point prior
%                          by specifying the minimal interesting difference
%
% data (matrix/vector):    You have two options here.
%                          EITHER data is a 1 x 2 vector containing the
%                          mean difference and SE of the difference...
%
%                          OR data is an nsubj x 1 vector of the difference
%                          score for each subject.
%
% prior_params             - If your prior is 'normal', this should be
%                              [prior_mean, prior_SD, number_of_tails]
%
%                          - If your prior is 'uniform', this should be
%                                   [lower_limit, upper_limit]
%
%                          - If you are testing the minimal interesting
%                          difference, this should be
%                                [minimal_interesting_difference]
%
% sensitive_BF:            What is the minimum BF you need to infer sensitive
%                          evidence for H1? The default value is 3, but
%                          you may choose e.g. 10. Note that if you select
%                          e.g. 10, then you will have a sensitive null if
%                          BF < 1/10
%
% plotit                 do you want to plot the RR? [1/0 or true/false]
%                            defaults to false
%
% timeout                how long should the search continue for.
%                            defaults to 3 seconds.
%
%
% EXAMPLE:
%
% I predict a difference of 13 in my data (in either direction), so my prior is a gaussian with mean 0 and SD 13.
% I look at my data. My mean difference is 10, and my SEdiff is 4.
% I want to use the default value for a sensitive BF (3)
% I want to plot the robustness region
% I don't want the robustness region search to go on for more than 5 seconds.
%
% So I would call the following:
% [BF, RR] = robustnessRegion( 'normal', [10, 4], [0 13 2], [], true, 5 )
%
%
% NOTE: you can run this with no inputs and the code will ask you to enter
% the relevant information at the command line
%
% =========================================================================
% v1.1 June 30th 2021 - updated so that you can edit sensitive_BF
%                     - brought back plotting feature
%                     - amended get_inputs so that these 2 parameters are
%                       collected too
%
% For questions/comments please email Maxine Sherman
% m.sherman@sussex.ac.uk / maxinesherman@gmail.com
% =========================================================================



function [BF_empirical RR] = robustnessRegion( shape_of_prior , data , prior_params, sensitive_BF, plotit, timeout  )

%% ------------------------------------------------------------------------
%   interpret data
% -------------------------------------------------------------------------

if nargin == 0
    [shape_of_prior , data , prior_params, sensitive_BF] = get_inputs;
end

if numel(data) == 2 % from summary stats
    bf_type = 'summary';
else
    bf_type = 'data';
    
    if size(data,2) > size(data,1)
        data = data';
    end
end

if nargin < 4 & nargin ~= 0; sensitive_BF = 3; end
if nargin < 5 & nargin ~= 0; plotit = false;   end
if nargin < 6;               timeout = 3;      end

%% ------------------------------------------------------------------------
%   if we're doing MID, flip to positive
% -------------------------------------------------------------------------
if strcmpi(shape_of_prior,'minimal')
    switch bf_type
        case 'summary'; if data(1) < 0;    data(1) = -data(1); end
        case 'data';    if mean(data) < 0; data = -data;       end
    end
end
%% ------------------------------------------------------------------------
%   if we're doing uniform, flip to positive
% -------------------------------------------------------------------------
if strcmpi(shape_of_prior,'uniform')
    switch bf_type
        case 'summary'; if data(1) < 0;    data(1) = -data(1);prior_params = abs(prior_params); end
        case 'data';    if mean(data) < 0; data = -data;  prior_params = abs(prior_params);     end
    end
end

%% ------------------------------------------------------------------------
%   the prior for H or N should never be negative
%   this is because the prior is in SD units.
%   if the prior + mean(data) have the same sign then make both positive.
%   if theprior + mean(data) have different signs, make the prior
%   positive and mean(data) negative.
% -------------------------------------------------------------------------
if strcmpi(shape_of_prior,'normal')
    
    switch bf_type
        case 'summary'
            if sign(data(1)) == sign(prior_params(2))
                data(1) = abs(data(1));
                prior_params(2) = abs(prior_params(2));
            else
                data(1) = -abs(data(1));
                prior_params(2) = abs(prior_params(2));
            end
            
        case 'data'
            % if both mean(data) and the prior are pos, do nothing
            if mean(data) > 0 & prior_params(2) > 0
                
                % if mean(data) < 0 and prior < 0, flip both
            elseif mean(data) < 0 & prior_params(2) < 0
                data            = -data;
                prior_params(2) = abs(prior_params(2));
                
                % if mean(data) > 0 but prior < 0, flip both (so data is in
                % opposite direction from prior)
            elseif mean(data) > 0 & prior_params(2) < 0
                data            = -data;
                prior_params(2) = abs(prior_params(2));
                
                % if mean(data) < 0 but prior > 0 then keep as is
            elseif mean(data)<0 & prior_params(2) > 0
                
            end
    end
end

%% ------------------------------------------------------------------------
%   cycle through and compute BF
% -------------------------------------------------------------------------

% compute real BF
switch bf_type
    case 'data'
        BF_empirical = getBF_2(0, shape_of_prior, data, prior_params);
    case 'summary'
        BF_empirical = getBF_from_summary_stats(0, shape_of_prior, data, prior_params);
end

switch shape_of_prior
    case 'minimal'; idx = 1; % only change MID
    case 'normal';  idx = 2; % only change SD
    case 'uniform'; idx = 2; % only change upper bound
end

% round prior to 2sf
prior_params(idx) = str2double(sprintf('%.2f',prior_params(idx)));

% -------------------------------------------------------------------------
%   Don't let the user request a non-zero mean
% -------------------------------------------------------------------------
if strcmpi(shape_of_prior,'normal')
    assert( prior_params(1)==0, 'This code only works when the mean of your Gaussian prior is 0. Sorry!');
end

%% ------------------------------------------------------------------------
%   compute RR for sensitive BFs
% -------------------------------------------------------------------------

% Check which prior makes it insensitive
if BF_empirical <= 1/sensitive_BF | BF_empirical >= sensitive_BF
    
    for i_dir = 1:2
        
        switch i_dir
            case 1; update_dir = 'decrease';
            case 2; update_dir = 'increase';
        end
        
        if BF_empirical < 1/sensitive_BF
            limit        = 1/sensitive_BF;
        elseif BF_empirical > sensitive_BF
            limit        = sensitive_BF;
        end
        
        % initialise things
        RR_prior       = prior_params;
        BF             = BF_empirical;
        current_sign   = sign(BF - limit);
        prev_sign      = sign(BF - limit);
        start_time     = tic;
        
        %% begin searching...
        while current_sign == prev_sign
            
            % update the prior
            switch update_dir
                case 'increase';  RR_prior(idx) = RR_prior(idx)*1.05;
                case 'decrease';  RR_prior(idx) = RR_prior(idx)*0.95;
            end
            
            % compute the new BF
            switch bf_type
                case 'data'
                    BF = getBF_2(0, shape_of_prior, data, RR_prior);
                case 'summary'
                    BF = getBF_from_summary_stats(0, shape_of_prior, data, RR_prior);
            end
            
            % find the error
            prev_sign             = current_sign;
            current_sign          = sign(BF - limit);
            
            % are we out of time?
            if toc(start_time) > timeout
                if i_dir == 1
                    RR(1) = -inf;
                else
                    RR(2) = inf;
                end
                break
            else
                RR(i_dir)= RR_prior(idx);
            end
            
        end
        
    end
    
    
    % Check which prior makes it sensitive
elseif BF_empirical > 1/sensitive_BF & BF_empirical < sensitive_BF
    
    % cycle through the 2 limits
    for i_dir = 1:2
        
        RR_prior = prior_params;
        
        switch i_dir
            case 1; update_dir = 'decrease';
            case 2; update_dir = 'increase';
        end
        
        % initialise things
        BF              = BF_empirical;
        current_sign   = [sign(BF - 1/sensitive_BF), sign(BF - sensitive_BF)];
        prev_sign      = [sign(BF - 1/sensitive_BF), sign(BF - sensitive_BF)];
        
        % in case it's inf, set a timer
        start_time = tic;
        
        %% searching...
        while sum(current_sign == prev_sign) == 2
            
            % update the prior
            switch update_dir
                case 'increase';  RR_prior(idx) = RR_prior(idx)*1.05;
                case 'decrease';  RR_prior(idx) = RR_prior(idx)*0.95;
            end
            
            % compute the new BF
            switch bf_type
                case 'data'
                    BF = getBF_2(0, shape_of_prior, data, RR_prior);
                case 'summary'
                    BF = getBF_from_summary_stats(0, shape_of_prior, data, RR_prior);
            end
            % find the error
            prev_sign             = current_sign;
            current_sign          = [sign(BF - 1/sensitive_BF), sign(BF - sensitive_BF)];
            
            % are we out of time?
            if toc(start_time) > timeout
                if i_dir == 1
                    RR(1) = -inf;
                else
                    RR(2) = inf;
                end
                break
            else
                RR(i_dir) = RR_prior(idx);
            end
        end
    end
    %% ------------------------------------------------------------------------
    %   other BFs???
    % -------------------------------------------------------------------------
    
else
    sprintf('BF = %.1f. cannot interpret. quitting...',BF_empirical)
end
%% ------------------------------------------------------------------------
%   Draw RR plot
% -------------------------------------------------------------------------


% get BFs
xlimit = max(10*prior_params(idx), max(~isinf(RR)));

P = linspace(0,xlimit);
for i = 1:100
    p      = prior_params;
    p(idx) = P(i);
    
    switch bf_type
        case 'data'
            bf(i)  = getBF_2(0, shape_of_prior, data, p);
        case 'summary'
            bf(i)  = getBF_from_summary_stats(0, shape_of_prior, data, p);
    end
end

% get RR for plotting
rr = RR; if isinf(abs(rr(1))); rr(1) = 0; end; if isinf(abs(rr(2))); rr(2) = xlimit; end

if nargin == 0
    plotit = input('do you want to plot your results?');
end

if plotit
    % initialise figure
    figure; hold on
    
    % move to log space?
    if max(bf) > 100
        bf     = log(bf);
        rr     = log(rr);
        bf_emp = log(BF_empirical);
        
        L1     = log(1/sensitive_BF);
        L2     = log(sensitive_BF);
        yname  = 'logBF';
    else
        yname  = 'BF';
        bf_emp = BF_empirical;
        
        L1     = 1/sensitive_BF;
        L2     = sensitive_BF;
    end

    % plot RR
    h1 = plot(linspace(rr(1),rr(2)),ones(100,1)*L1,'Color',[.8 .6 .6],'LineWidth', 8);
    plot(linspace(rr(1),rr(2)),ones(100,1)*L2,'Color',[.8 .6 .6],'LineWidth', 8);
    
    % plot reference lines
    h2 = plot(P, ones(size(P))*L1, 'k--','LineWidth',3);
         plot(P, ones(size(P))*L2, 'k--','LineWidth',3);
    
    % plot BF
    h3 = plot(P,bf,'Color',[.8 .3 .3],'LineWidth',3);
    
    % plot empirical BF
    h4 = scatter(prior_params(idx),bf_emp,200,'filled','MarkerFaceColor',[.7 .3 .3],'MarkerEdgeColor','k','LineWidth',3);
    
    % legend
    l  = legend([h4 h3 h1 h2],{'Your BF','BFs by prior','RR','Sensitive BF'},'FontSize',18,'LineWidth',3);
    
    % format
    xlabel('Prior','FontSize',20)
    ylabel(yname,'FontSize',20)
    set(gca,'FontSize',18,'TickLength',[0 0],'LineWidth', 3, 'YLim',[0 max(bf)],'XLim',[0 xlimit]);
    box on;
end


%% Report results

if nargin == 0
    disp(sprintf('Your BF is %.3f',BF_empirical));
    disp(sprintf('Your robustness region is [%.3f , %.3f]. A sensitive BF is > %d or < 1/%d',[RR(1) RR(2) sensitive_BF sensitive_BF]));
end

end

%% ========================================================================
%  [shape_of_prior , data , prior_params] = get_inputs;
%  ========================================================================

function [shape_of_prior , data , prior_params, sensitive_BF] = get_inputs

shape_of_prior = input('Is your prior normal (1), uniform (2) or minimal interesting difference (3)? : ');

switch shape_of_prior
    
    case 1
        shape_of_prior  = 'normal';
        prior_params(1,1) = input('Enter your prior mean: ');
        prior_params(1,2) = input('Enter your prior SD: ');
        prior_params(1,3) = input('Is your prior 1 or 2 tailed: ');
        
    case 2
        shape_of_prior  = 'uniform';
        prior_params(1,1) = input('What is the lower limit: ');
        prior_params(1,2) = input('What is the upper limit: ');
        
        
    case 3
        shape_of_prior  = 'minimal';
        prior_params(1) = input('Enter your minimum interesting difference: ');
        
end

% get data
data(1,1) = input('What is the mean difference in your data: ');
data(1,2) = input('What is the SE of the difference in your data: ');

sensitive_BF = input('What BF is the minimum to indicate sensitive support for H1?: ');
end




%% ========================================================================
%   [BF, L_theory, L_null] = getBF_2(print_report,test_type, data, prior)
%  ========================================================================


% [BF, L_theory, L_null] = getBF_2(print_report,test_type, data, prior)
%
% TEST TYPE:
%
% 'uniform'. varargin1 = lower bound;
%            varargin2 = upper bound
%
% 'normal'. varargin1 = prior mean
%           varargin2 = prior sd.
%           varargin3 = 1/2 tails
%
% 'minimal' varargin = minimal interesting difference.
%
%
%
% Pre-reg BFs: d' MID = 0.06
%              conf = 6%
%              mratio = U(0, 50)





function  [BF, L_theory, L_null] = getBF_2(print_report,test_type, data, prior)


%% Interpret inputs

switch test_type
    
    case 'uniform'
        
        lower = prior(1);
        upper = prior(2);
        
        isUniform = 1;
        
    case 'normal'
        
        meanoftheory    = prior(1);
        sdtheory        = prior(2);
        oneOrTwoTails   = prior(3);
        omega           = sdtheory*sdtheory;
        
        isUniform = 0;
        
    case 'minimal' % minimally interesting difference
        
        MID  = prior;
        
    otherwise
        
        error('unknown test_type. please use uniform, normal or minimal')
        
end

%% Interpret data
data_N  = size(data,1);
if size(data,2) == 1;
    data_SE = std(data)./sqrt(data_N);
    data_M  = mean(data);
else
    data_SE = std(data(:,2)-data(:,1))./sqrt(data_N);
    data_M  = mean(data(:,2)-data(:,1));
end

normaly = @(mn, variance, x) 2.718283^(- (x - mn)*(x - mn)/(2*variance))/realsqrt(2*pi*variance);
sd2     = data_SE*data_SE;
%% compute test

switch test_type
    
    case 'minimal'
        L_null   = normaly(0, sd2, data_M);
        L_theory = normaly(MID*2, sd2, data_M);
        BF       = L_theory/L_null;
        
    otherwise
        
        area = 0;
        if isUniform == 1
            theta = lower;
        else theta = meanoftheory - 5*(omega)^0.5;
        end
        if isUniform == 1
            incr = (upper- lower)/2000;
        else incr =  (omega)^0.5/200;
        end
        
        for A = -1000:1000
            theta = theta + incr;
            if isUniform == 1
                dist_theta = 0;
                if and(theta >= lower, theta <= upper)
                    dist_theta = 1/(upper-lower);
                end
            else %distribution is normal
                if oneOrTwoTails == 2
                    dist_theta = normaly(meanoftheory, omega, theta);
                else
                    dist_theta = 0;
                    if theta > 0
                        dist_theta = 2*normaly(meanoftheory, omega, theta);
                    end
                end
            end
            
            height = dist_theta * normaly(theta, sd2, data_M); %p(population value=theta|theory)*p(data|theta)
            area = area + height*incr; %integrating the above over theta
        end
        
        L_theory = area;
        L_null   = normaly(0, sd2, data_M);
        BF = L_theory/L_null;
        
        
end

%if print_report
%    sprintf('BayesFactor: %.2f | likelihood of theory: %.2f | Likelihood of null: %.2f',[BF,L_theory,L_null])
%end

end



%% ========================================================================
%  [BF, L_theory, L_null] = getBF_from_summary_stats(test_type, varargin)
%  ========================================================================

% [BF, L_theory, L_null] = getBF_from_summary_stats(test_type, varargin)
%
% TEST TYPE:
%
% 'uniform'. varargin1 = lower bound;
%            varargin2 = upper bound
%
% 'normal'. varargin1 = prior mean
%           varargin2 = prior sd.
%           varargin3 = 1/2 tails
%
% 'minimal' varargin = minimal interesting difference.
%
%
%
% Pre-reg BFs: d' MID = 0.06
%              conf = 6%
%              mratio = U(0, 50)





function  [BF, L_theory, L_null] = getBF_from_summary_stats(print_report,test_type, data, prior)


%% Interpret inputs

switch test_type
    
    case 'uniform'
        
        lower = prior(1);
        upper = prior(2);
        
        isUniform = 1;
        
    case 'normal'
        
        meanoftheory    = prior(1);
        sdtheory        = prior(2);
        oneOrTwoTails   = prior(3);
        omega           = sdtheory*sdtheory;
        
        isUniform = 0;
        
    case 'minimal' % minimally interesting difference
        
        MID  = prior;
        
    otherwise
        
        error('unknown test_type. please use uniform, normal or minimal')
        
end

%% Interpret data
data_M  = data(1);
data_SE = data(2);

normaly = @(mn, variance, x) 2.718283^(- (x - mn)*(x - mn)/(2*variance))/realsqrt(2*pi*variance);
sd2     = data_SE*data_SE;
%% compute test

switch test_type
    
    case 'minimal'
        L_null   = normaly(0, sd2, data_M);
        L_theory = normaly(MID*2, sd2, data_M);
        BF       = L_theory/L_null;
        
    otherwise
        
        area = 0;
        if isUniform == 1
            theta = lower;
        else theta = meanoftheory - 5*(omega)^0.5;
        end
        if isUniform == 1
            incr = (upper- lower)/2000;
        else incr =  (omega)^0.5/200;
        end
        
        for A = -1000:1000
            theta = theta + incr;
            if isUniform == 1
                dist_theta = 0;
                if and(theta >= lower, theta <= upper)
                    dist_theta = 1/(upper-lower);
                end
            else %distribution is normal
                if oneOrTwoTails == 2
                    dist_theta = normaly(meanoftheory, omega, theta);
                else
                    dist_theta = 0;
                    if theta > 0
                        dist_theta = 2*normaly(meanoftheory, omega, theta);
                    end
                end
            end
            
            height = dist_theta * normaly(theta, sd2, data_M); %p(population value=theta|theory)*p(data|theta)
            area = area + height*incr; %integrating the above over theta
        end
        
        L_theory = area;
        L_null   = normaly(0, sd2, data_M);
        BF = L_theory/L_null;
        
        
end

if print_report
    %    sprintf('BayesFactor: %.2f | likelihood of theory: %.2f | Likelihood of null: %.2f',[BF,L_theory,L_null])
end
end




