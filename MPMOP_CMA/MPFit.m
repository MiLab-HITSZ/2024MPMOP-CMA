function mFit = MPFit(x, fit, pro, DM, minRefer, maxRefer)

N = size(x, 1);
numGroups = DM;
mFit = zeros(N, 2*numGroups);

normFit = (fit - minRefer)/(maxRefer - minRefer);

dimIdx = dimension_group(pro.D, DM);

shiftX_dimIdx = cell(1, numGroups);
p = zeros(N, numGroups);

for i = 1:numGroups
    shiftX_dimIdx{i} = (x(:, dimIdx{i}) - pro.lower(dimIdx{i})) ./ (pro.upper(dimIdx{i}) - pro.lower(dimIdx{i}));
    p(:, i) = sum(shiftX_dimIdx{i}, 2) / length(dimIdx{i});
end

normFit = normFit(:);
for i = 1:numGroups
    mFit(:, 2*i-1) = p(:, i) + normFit;
    mFit(:, 2*i) = 1 - p(:, i) + normFit;
end


end


function groups = dimension_group(D, K)
    base_num = floor(D / K);
    remaining_elements = rem(D, K);
    
    groups = cell(1, K);
    idx = 1;
    for i = 1:remaining_elements
        groups{i} = idx : idx + base_num;
        idx = idx + base_num + 1;
    end

    for i = remaining_elements + 1:K
        groups{i} = idx : idx + base_num - 1;
        idx = idx + base_num;
    end
end
