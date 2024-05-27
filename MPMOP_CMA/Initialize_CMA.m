function S = Initialize_CMA(x, fit, pro)

[~, rank] = sort(fit);
x = x(rank, :);
fit = fit(rank);
seeds = nbc_seeds(x, fit);

seedX = x(seeds, :);
xk    = seedX;
sk    = 1:length(seeds);
fk    = fit(seeds);

% CMA parameter initialize
S =  struct('s',num2cell(sk)','x',num2cell(xk,2),'sigma',0.5,'C',eye(pro.D),'pc',0,'ps',0, 'bx', num2cell(xk,2), 'bf', num2cell(fk), 'valid', 1, 'cmaGen', 0, 'ter', 0);

fliterd = [S.bf] <= mean(fit);
S = S(fliterd);

end