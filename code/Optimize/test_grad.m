% test gradient

% initialize
xrange = linspace(1,220,101);
[X,Y] = meshgrid(xrange);
querypoints = [X(:),Y(:)]';
order = 1;
param_in = [ 0,0.0250];
param_0 = [0,0,0,0];
e_true = forward_chem(querypoints, param_in, order);

%% generate data from forward model
gradf = zeros(10, 2*order);
for i = 1:10
    gradf(i,:) = dfda(e_true,querypoints, param_in, order,10^(-2*i-4));
    display(i);
end

