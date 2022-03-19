function [log_post,y] = pdf(X,Obs,sd)

% Change pdf.m and log_prior_pdf.m if needed 

N = size(X,1);
log_lik = nan(N,1);
y = nan(N,length(Obs));

log_prior = log_prior_pdf(X);
parfor i = 1:N
    y(i,:) = forwardmodel(X(i,:),i);
    Err = Obs - y(i,:)';
    log_lik(i,1) = - ( length(Err) / 2) * log(2 * pi) - sum ( log( sd ) ) - 1/2 * sum ( ( Err./sd ).^2);
end

log_post = log_prior + log_lik;

end

function log_prior = log_prior_pdf(X)

% Calculate the log pdf of X in the prior distribution

N = size(X,1);
Npar = size(X,2);
log_prior = nan(N,1);
xmin = importdata('para_min.txt')';
xmax = importdata('para_max.txt')';
xmin = xmin(1:Npar);
xmax = xmax(1:Npar);

for i = 1:N
    temp = nan(Npar,1);
    for j = 1:Npar
        temp(j,1) = log(1/(xmax(j)-xmin(j)));
    end
    log_prior(i,1) = sum(temp);
end

end

