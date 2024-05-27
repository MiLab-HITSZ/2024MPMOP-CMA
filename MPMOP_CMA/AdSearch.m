function [x, fit, newS] = AdSearch(sol, x, fit, minorGen, DM, pro, algRand, rp, NP, ur)
if ~isempty(sol)
    addpop = reinit_history(sol, pro, rp, NP, ur, false);
    addfit = -pro.GetFits(addpop);
    [addpop, addfit] = MPSearch(addpop, addfit, size(addpop, 1), minorGen, DM, pro, algRand);

    x = [x; addpop];
    fit = [fit; addfit];
    newS = Initialize_CMA(x, fit, pro);
    [~, rank] = sort(fit);
    x = x(rank(1:NP), :);
    fit = fit(rank(1:NP));
else

    [x, fit] = MPSearch(x, fit, size(x, 1), minorGen, DM, pro, algRand);
    newS = Initialize_CMA(x, fit, pro);
end

end

