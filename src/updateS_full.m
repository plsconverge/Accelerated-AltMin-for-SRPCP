function [r, val, status, time_sort] = updateS_full(arr, rho)
    % full sort version for sparse subproblem

    time_sort = 0;

    if (min(arr) < 0.0)
        fprintf("Input array should be nonnegative\n");
        r = -1;
        val = nan;
        status = 1;
        return;
    end
    if (rho <= 0.0)
        fprintf("Parameter rho should be positive\n");
        r = -1;
        val = nan;
        status = 1;
        return;
    end

    n = length(arr);
    a_nnz = nnz(arr);
    a_nrm2 = norm(arr);
    a_max = max(arr);

    if (a_nnz == 0)
        r = 0;
        val = a_max + 1.0;
        status = 0;
        return;
    elseif (rho >= a_max / a_nrm2)
        r = 0;
        val = a_max + 1.0;
        status = 0;
        return;
    elseif (rho <= 1.0 / sqrt(a_nnz))
        r = a_nnz;
        val = 0.0;
        status = 0;
        return;
    end


    bound = 1.0 / (rho^2);
    k_bar = floor(bound);
    if (bound - k_bar < 1e-8) 
        k_bar = k_bar - 1;
    end
    k_bar = min([k_bar, a_nnz - 1]);

    sort_st = tic;
    arr = sort(arr, "descend");
    time_sort = time_sort + toc(sort_st);

    a_sq = norm(arr(k_bar+1 : end))^2;
    suffix_sq = cumsum(arr(2:k_bar).^2, "reverse");
    bi = arr(1:k_bar-1) ./ sqrt((1:k_bar-1)' .* arr(1:k_bar-1).^2 + suffix_sq + a_sq);
    r = sum(bi > rho);
    if (r <= 0)
        r = -1;
        val = nan;
        status = 2;
        return;
    end
    val = rho * sqrt((suffix_sq(r) + a_sq) / (1 - r * rho^2));
    status = 0;
end