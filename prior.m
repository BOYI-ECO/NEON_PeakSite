function X = prior(N,d,xmin,xmax)

% xmin = [1.0 0.10 0.10 0.00 0.10];
% xmax = [500 2.00 0.99 0.10 0.99];

X = nan(N,d);

for i = 1:N
    X(i,:) = unifrnd(xmin,xmax);   %normrnd(xmin,xmax,1,d)
end    

end