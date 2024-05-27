function Sigma = CMASearch(Sigma, k, pro)

ter = [Sigma.ter];

tidx = find(ter == true);

bx = cat(1, Sigma.bx);
bf = cat(1, Sigma.bf);


bxter = bx(tidx,:);
bfter = bf(tidx);


bxdis = pdist2(Sigma(k).bx, bxter);
% bxdis(k) = inf;
[minbxdis, minidx] = min(bxdis);
m = tidx(minidx);

if isempty(Sigma(k).x)
    
    Sigma(k).ter = 1;
    return
end

if ~isempty(tidx) && (minbxdis < 0.01 && bfter(m) < Sigma(k).bf)
     Sigma(k).valid = 0;
     Sigma(k).ter = 1;
     return
end

mu = 4 + floor(3 * log(pro.D));
OffsX = mvnrnd(Sigma(k).x, Sigma(k).sigma^2 * Sigma(k).C, mu);
OffsX = boundary_check(OffsX, pro.lower, pro.upper);
OffsFit = -pro.GetFits(OffsX);

if pro.change
    return
end


[~, rank] = sort(OffsFit);
Sigma(k) = UpdateCMA(OffsX(rank, :), Sigma(k));

combineOffsFit = [OffsFit; Sigma(k).bf];
combineOffsX = [OffsX; Sigma(k).bx];
[minbf, minIdx] = min(combineOffsFit);

Sigma(k).bf = minbf;
Sigma(k).bx = combineOffsX(minIdx, :);
end