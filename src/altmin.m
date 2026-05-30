function [L, S, res] = altmin(D, lambda, mu, params)
    % main function for AltMin algorithm
    % solves min ||L||_* + lambda*||S||_1 + mu*||L+S-D||_F
    % 
    % Inputs: 
    % D -- observation matrix, dim m*n
    % lambda, mu -- parameters, lambda = 1 / sqrt(m), mu = sqrt(n/2)
    % params -- controllers
    %
    % Outputs:
    % L -- recovered low-rank matrix
    % S -- recovered sparse matrix
    % res -- result messages
    %
    % implement of this algorithm partly refers to 
    % Deng, S., Li, X., & Zhang, Y. (2025). Alternating minimization for square root principal component pursuit. INFORMS Journal on Computing.


    % start timer
    time_st = tic;

    % set params
    [n1, n2] = size(D);
    tol = 1e-5;
    max_iter = 5e3;
    max_time = 18000;

    trans_flag = false;
    if (n1 < n2)
        D = D';
        fprintf("Warning: Row < Col, Use transposition.\n");
        [n1, n2] = size(D);
        trans_flag = true;
    end

    rank_guess = round(n2 / 2);
    sparse_guess = round(n1 * n2 / 2);
    buffer = 0.1;

    % default mode
    Lmod = "partialsvd";
    Smod = "partialsort";
    display = true;

    % parser
    if (isfield(params, "Lmod")) Lmod = params.Lmod; end
    if (isfield(params, "Smod")) Smod = params.Smod; end
    if (isfield(params, "tol")) tol = params.tol; end
    if (isfield(params, "disp")) display = params.disp; end
    if (isfield(params, "maxiter")) max_iter = params.maxiter; end

    % info
    if display
        fprintf("Processing Matrix of Size %dx%d...\n", n1, n2);
        fprintf(" Iter        Obj            Eta        rank    sparse\n");
    end

    % initialization & memory allocation
    L = zeros(n1, n2);
    S = zeros(n1, n2);
    res = struct("status", -1, "obj", 0.0, "time", 0.0);
    res.time_detail = struct("L", 0.0, "S", 0.0, "sort", 0.0);
    res.shist = zeros(max_iter, 1);

    for iter = 1 : max_iter
        % update S
        time_st_S = tic;
        S = D - L;
        [S, status, k, l1nrm, time_sort] = updateS(S, lambda / mu, Smod, sparse_guess);
        res.time_detail.S = res.time_detail.S + toc(time_st_S);
        res.time_detail.sort = res.time_detail.sort + time_sort;
        res.shist(iter) = nnz(S);

        if (status)
            res.status = status;
            return;
        end
        sparse_guess = k;

        % update L
        time_st_L = tic;
        L = D - S;
        [L, status, k, nuc_nrm] = updateL(L, 1.0 / mu, Lmod, rank_guess);
        res.time_detail.L = res.time_detail.L + toc(time_st_L);

        if (status)
            res.status = status;
            return;
        end
        rank_guess = round(k * (1.0 + buffer));
        rank_guess = min([rank_guess, n2]);

        % diagnostics
        [obj, eta] = diagnostics(D, L, S, lambda, mu, nuc_nrm, l1nrm);

        % info
        if display
            fprintf("%4d   %.8f    %.8f    %3d    %d\n", iter, obj, eta, k, sparse_guess);
        end

        % terminal check
        if (eta < tol)
            fprintf("Iter %d -- Converge, eta = %.8f < tol = %.8f\n", iter, eta, tol);
            res.status = 0;
            res.obj = obj;
            res.time = toc(time_st);

            if trans_flag
                L = L';
                S = S';
            end

            return;
        end

        % time check
        time_check = toc(time_st);
        if (time_check > max_time)
            fprintf("Time Limited, Solver Terminated\n");
            res.status = 2;
            res.obj = obj;
            res.time = toc(time_st);
            return;
        end
    end
    % max iterations reached
    fprintf("Max Iterations Reached, Solver Terminated\n");
    res.status = 2;
    res.obj = obj;
    res.time = toc(time_st);
    return;
end

function [obj, eta] = diagnostics(D, L, S, lambda, mu, nuc_nrm, l1nrm)
    % diagnostics, return criterion eta and objective obj
    nrm_L = norm(L, "fro");
    nrm_S = norm(S, "fro");

    t = L + S - D;
    nrm_fro = norm(t, "fro");
    t = S - mu * t / nrm_fro;
    t = sign(t) .* max(abs(t) - lambda, 0.0);

    eta = norm(S - t, 'fro') / max([1.0, nrm_S, nrm_L]);
    obj = mu * nrm_fro + nuc_nrm + lambda * l1nrm;
end
