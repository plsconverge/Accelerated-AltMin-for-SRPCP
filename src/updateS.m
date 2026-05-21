function [S, status, k, l1nrm, time_sort] = updateS(S, tau, mode, k0)
    % status = 0 -- converge
    %          1 -- invalid mode
    %          2 -- subproblem fail
    if (strcmp(mode, "fullsort"))
        [S, status, k, l1nrm, time_sort] = updateS_ful(S, tau);
    elseif (strcmp(mode, "partialsort"))
        [S, status, k, l1nrm, time_sort] = updateS_par(S, tau, k0);
    else
        fprintf("Warning: Invalid mode for updating S\n");
        status = 1;
        k = -1;
        l1nrm = 0.0;
    end
end


function [S, status, k, l1nrm, time_sort] = updateS_ful(S, tau)
    sign_s = sign(S);
    S = abs(S);
    arr = S(:);

    [k, tk, status, time_sort] = updateS_full(arr, tau);
    if (status ~= 0)
        S = nan;
        l1nrm = nan;
        return;
    end

    % soft threshold
    S = sign_s .* max(S - tk, 0.0);
    l1nrm = norm(S(:), 1);
end

function [S, status, k, l1nrm, time_sort] = updateS_par(S, tau, k0)
    sign_s = sign(S);
    S = abs(S);
    arr = S(:);

    [k, tk, status, time_sort] = updateS_acc(arr, tau, k0);
    if (status ~= 0)
        S = nan;
        l1nrm = nan;
        return;
    end

    % soft threshold
    S = sign_s .* max(S - tk, 0.0);
    l1nrm = norm(S(:), 1);
end
