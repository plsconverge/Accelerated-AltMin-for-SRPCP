function [S, status, k, l1nrm, time_sort] = updateS(S, rho, mode, k0)
    % control function for sparse subproblem

    % status = 0 -- converge
    %          1 -- invalid mode
    %          2 -- subproblem fail
    if (strcmp(mode, "fullsort"))
        [S, status, k, l1nrm, time_sort] = updateS_ful(S, rho);
    elseif (strcmp(mode, "partialsort"))
        [S, status, k, l1nrm, time_sort] = updateS_par(S, rho, k0);
    else
        fprintf("Warning: Invalid mode for updating S\n");
        status = 1;
        k = -1;
        l1nrm = 0.0;
    end
end


function [S, status, k, l1nrm, time_sort] = updateS_ful(S, rho)
    % outer processor
    sign_s = sign(S);
    S = abs(S);
    arr = S(:);

    [k, tk, status, time_sort] = updateS_full(arr, rho);
    if (status ~= 0)
        S = nan;
        l1nrm = nan;
        return;
    end

    % soft threshold
    S = sign_s .* max(S - tk, 0.0);
    l1nrm = norm(S(:), 1);
end

function [S, status, k, l1nrm, time_sort] = updateS_par(S, rho, k0)
    % outer processor
    sign_s = sign(S);
    S = abs(S);
    arr = S(:);

    [k, tk, status, time_sort] = updateS_acc(arr, rho, k0);
    if (status ~= 0)
        S = nan;
        l1nrm = nan;
        return;
    end

    % soft threshold
    S = sign_s .* max(S - tk, 0.0);
    l1nrm = norm(S(:), 1);
end
