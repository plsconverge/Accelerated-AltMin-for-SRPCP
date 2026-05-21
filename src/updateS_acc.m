function [r, val, status, time_sort] = updateS_acc(arr, rho, k0)
% solves min \|s-a\|_2 + \rho \|s\|_1
% assumes a >= 0, rho > 0
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

    % n = length(arr);
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
    t_bar = floor(bound);
    if (bound - t_bar < 1e-8)
        t_bar = t_bar - 1;
    end
    t_bar = min([t_bar, a_nnz - 1]);

    minWidth = 100;
    maxWidth = 5e4;
    width = min([maxWidth, max([minWidth, floor(0.1 * k0)])]);
    k1 = max([1, k0 - width]);
    k2 = min([t_bar, k0 + width]);

    sort_st = tic;
    arr = biquickselect(arr, k1, k2);
    time_sort = time_sort + toc(sort_st);

    sq = sum(arr(k2+1:end).^2);
    sq1 = sum(arr(k1+1:k2).^2) + sq;
    bk1 = arr(k1) / sqrt(k1 * arr(k1)^2 + sq1);
    bk2 = arr(k2) / sqrt(k2 * arr(k2)^2 + sq);

    if (bk1 > rho && bk2 <= rho) 
        % index r involved
        sort_st = tic;
        arr(k1+1:k2-1) = sort(arr(k1+1:k2-1), "descend");
        time_sort = time_sort + toc(sort_st);

        suffix_sq = cumsum(arr(k1+1:k2).^2, "reverse");
        bi = arr(k1:k2-1) ./ sqrt((k1:k2-1)' .* arr(k1:k2-1).^2 + suffix_sq + sq);
        r = k1 - 1 + sum(bi > rho);
        val = rho * sqrt((suffix_sq(r-k1+1) + sq) / (1 - r * rho^2));
        status = 0;
        return;
    elseif (bk1 <= rho)
        % exponential expansion -- direction -> 1
        alpha = 3;
        expansion = max([200, min([1e5, floor(0.2 * k0)])]);
        resort_count = 0;
        lidx = k1;
        lidx_old = k1;
        sq = sq1;
        while lidx > 1
            resort_count = resort_count + 1;

            lidx = max([1, k1 - expansion]);
            sort_st1 = tic;
            arr(1:lidx_old-1) = quickselect(arr(1:lidx_old-1), lidx);
            time_sort = time_sort + toc(sort_st1);

            sq1 = sq1 + sum(arr(lidx+1:lidx_old).^2);
            bk = arr(lidx) / sqrt(lidx * arr(lidx)^2 + sq1);
            if (bk > rho)
                % r involved, traverse to search
                sort_st2 = tic;
                arr(lidx+1:lidx_old-1) = sort(arr(lidx+1:lidx_old-1), "descend");
                time_sort = time_sort + toc(sort_st2);

                suffix_sq = cumsum(arr(lidx+1:lidx_old).^2, "reverse");
                bi = arr(lidx:lidx_old-1) ./ sqrt((lidx:lidx_old-1)' .* arr(lidx:lidx_old-1).^2 + suffix_sq + sq);
                r = lidx - 1 + sum(bi > rho);
                val = rho * sqrt((suffix_sq(r-lidx+1) + sq) / (1 - r * rho^2));
                status = 0;
                return;
            end

            % r not found, expand
            expansion = expansion * alpha;
            lidx_old = lidx;
            sq = sq1;
        end
    else
        % exponential expansion -- direction -> n
        alpha = 3;
        expansion = max([200, min([1e5, floor(0.2 * k0)])]);
        resort_count = 0;
        ridx = k2;
        ridx_old = k2;
        sq1 = sq;
        while ridx < t_bar
            resort_count = resort_count + 1;
            
            ridx = min([t_bar, k2 + expansion]);
            sort_st1 = tic;
            arr(ridx_old+1:end) = quickselect(arr(ridx_old+1:end), ridx - ridx_old);
            time_sort = time_sort + toc(sort_st1);

            sq = sq - sum(arr(ridx_old+1:ridx).^2);
            bk = arr(ridx) / sqrt(ridx * arr(ridx)^2 + sq);
            if (bk <= rho)
                % r involved, traverse to search
                sort_st2 = tic;
                arr(ridx_old+1:ridx-1) = sort(arr(ridx_old+1:ridx-1), "descend");
                time_sort = time_sort + toc(sort_st2);

                suffix_sq = cumsum(arr(ridx_old+1:ridx).^2);
                bi = arr(ridx_old+1:ridx) ./ sqrt((ridx_old+1:ridx)' .* arr(ridx_old+1:ridx).^2 - suffix_sq + sq1);
                r = ridx_old + sum(bi > rho);
                if (r == ridx_old)
                    val = rho * sqrt(sq1 / (1 - r * rho^2));
                else
                    val = rho * sqrt((sq1 - suffix_sq(r-ridx_old)) / (1 - r * rho^2));
                end
                status = 0;
                return;
            end

            % r not found, expand
            expansion = expansion * alpha;
            ridx_old = ridx;
            sq1 = sq;
        end
    end

    % failed due to unknown error
    r = -1;
    val = nan;
    status = 2;
end
