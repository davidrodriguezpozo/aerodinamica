function [alpha_zl] = compute_alphazl (N, X, rho, c)
    for alpha = deg2rad(3) : deg2rad(-.25) : deg2rad(-9)
        cl = compute_cl(alpha, N, X, rho, c);
        if cl<=0
            alpha_zl = alpha;
            break
        end
    end
end