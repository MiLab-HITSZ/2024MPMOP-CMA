function Rank = MPRank(x, PopObj, DM)

N = size(PopObj,1);
numGroups = DM;

%% Detect the dominance relation between each two solutions
DominateCell = cell(1, numGroups);
for i = 1:numGroups
    DominateCell{i} = false(N);
end

for i = 1:N-1
    for j = 1:numGroups
        k_less = any(PopObj(i,(2*j-1):(2*j)) < PopObj(i+1:end,(2*j-1):(2*j)), 2);
        k_greater = any(PopObj(i,(2*j-1):(2*j)) > PopObj(i+1:end,(2*j-1):(2*j)), 2);
        DominateCell{j}(i, i+find(k_less & ~k_greater)) = true;
        DominateCell{j}(i+find(~k_less & k_greater), i) = true;
    end
end

%% judge the distance
matdis = pdist2(x, x);
matdis(logical(eye(N))) = inf;
mindis = mean(min(matdis));
d =2.0* mindis;
matdis = matdis < d;

for j = 1:numGroups
    DominateCell{j} = DominateCell{j}&matdis;
end

RankCell = cell(1, numGroups);
for j = 1:numGroups
    RankCell{j} = zeros(N, 1);
end

BeDom_np = zeros(numGroups, N);
for i = 1:N
    for j = 1:numGroups
        BeDom_np(j, i) = sum(DominateCell{j}(:, i));
        if BeDom_np(j, i) == 0
            RankCell{j}(i) = 1;
        end
    end
end

for j = 1:numGroups
    level = 1;
    while any(RankCell{j} == 0)
        get_index = (find(RankCell{j} == level))';
        for k = get_index
            sp_index = find(DominateCell{j}(k, :) == 1);
            BeDom_np(j, sp_index) = BeDom_np(j, sp_index) - 1;
            NextLevel = BeDom_np(j, sp_index) == 0;
            RankCell{j}(sp_index(NextLevel)) = level + 1;
        end
        level = level + 1;
    end
end
Rank = cell2mat(RankCell);
end

