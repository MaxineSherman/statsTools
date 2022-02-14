function [m_posterior, s_posterior, posterior_prob] = bayesian_integration( prior , likelihood , donotplot)

if nargin == 2; donotplot = 0; end

% intepret inputs
m_prior   = prior(1);
s_prior   = prior(2);

m_data    = likelihood(1);
s_data    = likelihood(2);

% get precisions
precision_prior      = 1./s_prior;
precision_data        = 1./s_data;

% compute posterior precision
precision_posterior  = precision_prior + precision_data;
s_posterior          = 1./precision_posterior;

% compute posterior mean
m_posterior          = (m_prior*precision_prior/precision_posterior) ...
                      + (m_data*precision_data/precision_posterior);
                  
                  
posterior_prob      = makedist('normal',m_posterior,s_posterior);   
posterior_prob      = posterior_prob.pdf(m_posterior);

% plot
if ~donotplot
    
X     = linspace(-10,10);
PRIOR = normpdf(X,m_prior,s_prior);
DATA  = normpdf(X,m_data,s_data);
POST  = normpdf(X,m_posterior,s_posterior);

lw    = 2;
fs    = 14;

hold on;
plot(X,PRIOR,'Color',[0.8 0 0.8],'LineWidth',lw);
plot(X,DATA,'Color','k','LineWidth',lw);
plot(X,POST,'Color',[0 0.8 0],'LineWidth',lw*1.5);

set(gca,'LineWidth',lw,'TickLength',[0 0]);%,'YLim',[0 1.5]);
legend('Prior','Likelihood','Posterior')
end