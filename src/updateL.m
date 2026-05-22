function [L, status, k, nrm] = updateL(L, rho, mode, k0)
    % control function for low-rank subproblem solver

    % status = 0 -- converge
    %          1 -- invalid mode
    %          2 -- subproblem fail
    if (strcmp(mode, "fullsvd"))
        [L, status, k, nrm] = updateL_full(L, rho);
    elseif (strcmp(mode, "partialsvd"))
        [L, status, k, nrm] = updateL_par(L, rho, k0);
    else
        status = 1;
        k = -1;
        nrm = 0.0;
    end
end


function [L, status, k, nrm] = updateL_full(L, rho)
    % outer processor
    [U, D, V] = svd(L, "econ");
    D_diag = diag(D);

    [status, tk, k] = explicit_solver(D_diag, rho);

    D_diag = max(D_diag - tk, 0.0);
    L = U(:, 1:k) * (D_diag(1:k) .* V(:, 1:k)');
    nrm = sum(D_diag(1:k));
end

function [status, tk, k] = explicit_solver(arr, tau)
    % inner solver for UpdateL -- full version
    a_nnz = nnz(arr);
    a_nrm2 = norm(arr, 2);
    a_max = arr(1);

    if (a_nnz == 0)
        status = 0;
        tk = a_max + 1.0;
        k = 0;
        return;
    elseif (tau >= a_max / a_nrm2)
        status = 0;
        tk = a_max + 1.0;
        k = 0;
        return;
    elseif (tau <= 1.0 / sqrt(a_nnz))
        status = 0;
        tk = 0.0;
        k = a_nnz;
        return;
    end

    t2inv = 1.0 / (tau^2);
    k_bar = floor(t2inv);
    if (t2inv - k_bar < 1e-8) 
        k_bar = k_bar - 1;
    end
    a_sq = sum(arr(k_bar+1:end).^2);

    suffix_sq = cumsum(arr(2:k_bar).^2, "reverse");
    bi = arr(1:k_bar-1) ./ sqrt((1:k_bar-1)' .* arr(1:k_bar-1).^2 + suffix_sq + a_sq);
    k = sum(bi > tau);
    tk = tau * sqrt((suffix_sq(k) + a_sq) / (1 - k * tau^2));
    status = 0;

    % flag = true;
    % for k_test = k_bar : -1 : 1
    %     tk = sqrt(a_sq / (t2inv - k_test));
    %     if (arr(k_test) > tk && arr(k_test+1) <= tk)
    %         flag = false;
    %         k = k_test;
    %         status = 0;
    %         break;
    %     end
    %     a_sq = a_sq + arr(k_test)^2;
    % end
    % 
    % if (flag)
    %     fprintf("Warning: Failed to find k\n");
    %     status = 2;
    %     k = -1;
    %     tk = 0.0;
    % end
    % return;
end

function [L, status, k, nrm] = updateL_par(L, rho, k0)
    % outer processor 
    [~, n] = size(L);
    mult = 2;
    nrm_fro = norm(L, "fro");

    re_count = 0;
    while true
        [U, D, V] = svds(L, k0, "largest");
        D_diag = diag(D);

        [status, tk, k] = explicit_solver_par(D_diag, rho, k0, nrm_fro);

        if (status == 0)
            break;
        end
        re_count = re_count + 1;
        if (re_count > 5)
            fprintf("Re-solving too many times, loop terminated\n");
            status = 2;
            k = -1;
            nrm = 0.0;
            return;
        end

        k0 = k0 * mult;
        k0 = min([k0, n]);
    end

    D_diag = max(D_diag - tk, 0.0);
    nrm = sum(D_diag(1:k));
    L = U(:, 1:k) * (D_diag(1:k) .* V(:, 1:k)');
end

function [status, tk, k] = explicit_solver_par(arr, tau, k0, a_sq)
    % inner solver for UpdateL -- partial version
    a_max = arr(1);
    if (a_max < 1e-10)
        status = 0;
        tk = a_max + 1.0;
        k = 0;
        return;
    elseif (tau < 1e-10)
        status = 0;
        tk = 0.0;
        k = k0;
        return;
    elseif (tau >= a_max / a_sq)
        status = 0;
        tk = a_max + 1.0;
        k = 0;
        return;
    end

    t2inv = 1.0 / (tau^2);
    k_bar = floor(t2inv);
    if (t2inv - k_bar < 1e-8) 
        k_bar = k_bar - 1;
    end
    k_bar = min([k_bar, k0 - 1]);
    a_sq = a_sq * a_sq;

    % suffix_sq = cumsum(arr(1:k_bar-1).^2);
    % bi = arr(1:k_bar-1) ./ sqrt((1:k_bar-1)' .* arr(1:k_bar-1).^2 - suffix_sq + a_sq);
    % k = sum(bi > tau);
    % tk = tau * sqrt((a_sq - suffix_sq(k)) / (1 - k * tau^2));
    % status = 0;

    flag = true;
    for k_test = 1 : k_bar
        a_sq = a_sq - arr(k_test) * arr(k_test);
        tk = sqrt(a_sq / (t2inv - k_test));

        if (arr(k_test) > tk && arr(k_test+1) <= tk)
            flag = false;
            k = k_test;
            status = 0;
            break;
        end
    end

    if (flag)
        status = 2;
        k = -1;
        tk = 0.0;
    end
    return;
end
