function [Teager] = TE(x)

% Calculates mean Teager energy when given a signal x

% initialize variable
Teager = 0;
k = 0;
% Sum 
for i = 3:length(x)
    k = k+1;
    Teager = Teager + (x(i-1)^2 - x(i) * x(i-2));
end
Teager = Teager/k;

%% Take log to normalize feature
Teager = log(Teager);