f(x) = sin(x) + atan(x)
x0 = 1

k_max = 100 # Maximum number of iterations.
tol_res = 1e-6 # Tolerance on the residual.
m = 3 # Parameter m.

x = [x0, f(x0)] # Vector of iterates x.
g = f(x) - x # Vector of residuals.

G_k = g[2] - g[1] # Matrix of increments in residuals.
X_k = x[2] - x[1] # Matrix of increments in x.

k = 2
while k < k_max && abs(g[k]) > tol_res
    m_k = min(k, m);
 
    # Solve the optimization problem by QR decomposition.
    [Q, R] = qr(G_k)
    gamma_k = R \ (Q' * g(k))
 
    # Compute new iterate and new residual.
    push!(x, x[k] + g[k] - (X_k + G_k) * gamma_k);
    push!(g, f(x[k + 1]) - x[k + 1]);
 
    # Update increment matrices with new elements.
    X_k = [X_k..., x[k + 1] - x[k]]
    G_k = [G_k..., g[k + 1] - g[k]]
 
    n = size(X_k, 2)
    if n > m_k
        X_k = X_k[:, n - m_k + 1:end];
        G_k = G_k[:, n - m_k + 1:end];
    end
 
    k = k + 1;
end

% Prints result: Computed fixed point 2.013444 after 9 iterations
fprintf("Computed fixed point %f after %d iterations\n", x(end), k);
