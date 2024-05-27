function pop = reinit_history(archive, pro, rp, NP, ur, change)
    persistent data  linktable hisX
    
    if isempty(archive)
        return
    end

    if change
        hisX = [];
        [data, linktable] = history_data(archive);
        
        for i = 1: length(data)
            tempX = data{i}(:, 1:end-1)';
            tempFit = -pro.GetFits(tempX);
            if pro.change
                return
            end
           [~, minIdx] = min(tempFit);
           hisX = [hisX; tempX(minIdx, :)];
        end
    end
        

    
    x2 =[];
    if size(data{1}, 2) > 1
       diffF = 0.5;
       xt1 = cell2mat(cellfun(@(x) x(:, end)', data, 'UniformOutput', false));
       xt2 = cell2mat(cellfun(@(x) x(:, end-1)', data, 'UniformOutput', false));  
       x2 = xt1 + diffF .* (xt1 - xt2);
       x2 = boundary_check(x2, pro.lower, pro.upper);
    end
    
    x1 = cell2mat(cellfun(@(x) x(:, end)', data, 'UniformOutput', false));
    x1 = boundary_check(x1, pro.lower, pro.upper);
    
    x3 = hisX;
 
    x4 = [];
    for i = length(linktable):-1: 1
        curArchive = archive{i};
        curIdx = 1: size(curArchive, 1);
        curUnlinkIdx = setdiff(curIdx, linktable{i});
        if ~isempty(curUnlinkIdx)
            x4 = [x4; curArchive(curUnlinkIdx, :)];
        end
        
    end
    
    R = [x1; x2; x3; x4];

    variance_mat = 0.01 * eye(pro.D);    
    offsets = mvnrnd(zeros(1, pro.D), variance_mat, size(R, 1) * rp);

    pop = repmat(R, rp, 1) + offsets;
    pop = boundary_check(pop, pro.lower, pro.upper);
    pop = pop(1:min(size(pop, 1), ceil(NP * ur)), :);
end


function [current_history, linktable] = history_data(Archive)
    current_archive = Archive{end};
    num_solutions = size(current_archive, 1);
    current_history = cell(num_solutions, 1);
    linktable = cell(length(Archive), 1); 
    linktable{length(Archive)} = 1:size(current_archive, 1);
    for m = 1:num_solutions
        current_solution = current_archive(m, :);
        history = current_solution;

        for prev_t = length(Archive)-1:-1:1
            prev_archive = Archive{prev_t};        
            distances = sqrt(sum((prev_archive - current_solution).^2, 2));
            [~, nearest_index] = min(distances);
            history = [prev_archive(nearest_index, :); history];
            current_solution = prev_archive(nearest_index, :);
            linktable{prev_t} = [linktable{prev_t}, nearest_index];
        end

        current_history{m} = history';
    end
end
