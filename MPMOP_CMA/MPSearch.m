function [x, fit, FrontNo] = MPSearch(x, fit, NP, maxGen, DM, pro, algRand)
    curGen = 0;
    minRefer = min(fit);
    maxRefer = max(fit);
    while curGen < maxGen
        F  = .5;
        CR = .7;
        
        % Implement DE operators to generate the offspring
        offsX = DE(x, NP, F, CR, pro, algRand);
        offsFit = -pro.GetFits(offsX);
        
        if pro.change
            return
        end
        
        if min(offsFit) < minRefer
            minRefer = min(offsFit);
        end
        
        if max(offsFit) > maxRefer
            maxRefer = max(offsFit);
        end
        
        % Parents + Offspring
        x = [x; offsX];
        fit = [fit; offsFit];
        
        mFit = MPFit(x, fit, pro, DM, minRefer, maxRefer);
        
        % Get the Rank for each DM
        rank = MPRank(x, mFit, DM);

         if DM == 1
            FrontNo = inf(size(x, 1), 1);
            idx = [];
            r = 1;
            MaxFNo = 1;
            while(length(idx) < NP)
                MaxFNo = r;
                temp = find(rank == r);
                FrontNo(temp) = r;
                idx = [idx; temp];
                r = r+1;
            end
        else
            [FrontNo,MaxFNo] = NDSort(rank,NP);
        end
        
        CrowdDis = CrowdingDistance(mFit,FrontNo);

        Next = FrontNo < MaxFNo; 
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:NP-sum(Next)))) = true;
        
        x = x(Next,:);
        fit = fit(Next, :);
        FrontNo = FrontNo(Next);
        curGen = curGen + 1;
    end
end
