function ex_mpmop()
test_func = 1:24;

if isempty(gcp('nocreate'))
    parpool(30);
end
    pr_total = [];
    fname = sprintf('./MPMOP_CMA/result/ALG');
    if ~exist(fname, 'dir')
        mkdir(fname);
    end
    for func = test_func
        spmd(30)
            [peak, speak, sallpeak] = MPMOP_CMA(func, spmdIndex);
        end
        peak = cat(1, peak{1:end});
        sallpeak = cat(1, sallpeak{1:end});
        result = cat(2, speak{1:end});
        pr = mean(result') / sallpeak(1);
        dlmwrite(sprintf('./MPMOP_CMA/result/ALG/P%d.csv', func), peak);
        dlmwrite(sprintf('./MPMOP_CMA/result/ALG/F%d.csv', func), result);
        pr_total = [pr_total; pr];
    end
    dlmwrite(sprintf(['./MPMOP_CMA' ...
        '/result/ALG/F_total.csv']), pr_total);
delete(gcp);
end
