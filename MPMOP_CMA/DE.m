function offsX = DE(x, N, F, CR, pro, algRand)

D = pro.D;
 % mutation, DE/rand/1
index = zeros(N, 3);
for i = 1:N
    index(i, :) = randperm(algRand, N-1, 3);
    index(i, index(i, :) >= i) = index(i, index(i, :) >= i) + 1;
end

mutant = x(index(:, 1), :) - F .* (x(index(:, 2), :) - x(index(:, 3), :));

% crossover
cross = rand(algRand, N, D) < CR;
parent_only_idx = find(sum(cross, 2) == 0);
for i = parent_only_idx
    cross(i, randi(algRand, D)) = true;
end
offsX = cross .* mutant + (1-cross) .* x;

% bound check
offsX = boundary_check(offsX, pro.lower, pro.upper);
end
