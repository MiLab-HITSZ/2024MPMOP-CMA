function [peak, speak, sallpeak] = MPMOP_CMA(fn, run)

    algRand = RandStream.create('mt19937ar', 'Seed', run);
    RandStream.setGlobalStream(algRand);

    pro = DMMOP(fn);
    NP = [250, 250, 250, 250, 500, 500, 100, 500*ones(1, 9), 300, 300, 300, 300, 500, 500, 500, 500];
    NP = NP(fn);
    D = pro.D;
    x = rand(NP, D) .* (pro.upper - pro.lower) + pro.lower;
    fit = -pro.GetFits(x);
    DM = 2;

    % initilize parameters
    ar = 0.2;
    tw = 20;
    rp = 5;
    ur = 0.5;
    sol = {};
    minorGen = 3;
	

    while ~pro.Terminate()
		fprintf('MPMOP-CMA runing, funcNo:%d, run index:%d, environment:%d\n', fn, run, pro.env+1);
        rest = pro.freq - rem(pro.evaluated, pro.freq);
        maxGen = floor(rest / NP * ar);

        % multiparty multiobjective optimization
        [x, fit] = MPSearch(x, fit, NP, maxGen, DM, pro, algRand);
        

        % CMA-ES search stage
        S = Initialize_CMA(x, fit, pro);
        bx = cat(1, S.bx);
        bf = cat(1, S.bf);

        while ~pro.CheckChange(bx, -bf)
            ter = [S.ter];
            if all(ter)

                % additional search stage
                [x, fit, newS] = AdSearch(sol, x, fit, minorGen, DM, pro, algRand, rp, NP, ur);
                S = [S; newS];
            end
            
            ter = [S.ter];
            k = find(ter == false, 1);
            S = CMASearch(S, k, pro);
            bx = cat(1, S.bx);
            bf = cat(1, S.bf);
        end

        % dynamic response strategies
        valid_indices = [S.valid] == 1 & [S.cmaGen] >= 1;
        bx = cat(1, S(valid_indices).bx);
        bf = cat(1, S(valid_indices).bf);


        [~, rank] = sort(bf);
        bx = bx(rank, :);
        sol = [sol, {bx}];
        
        if length(sol) > tw
            sol(1) = [];
        end

        addpop = reinit_history(sol, pro, rp, NP, ur, true);
        x = rand(NP, pro.D) .* (pro.upper - pro.lower) + pro.lower;
        x(1:size(addpop, 1), :) = addpop;
        fit = -pro.GetFits(x);
    end

    [peak, allpeak] = pro.GetPeak();
    speak = sum(peak, 2);
    sallpeak = sum(allpeak);
end
