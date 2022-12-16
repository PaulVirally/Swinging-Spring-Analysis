function euler_lagrange_ode!(du, u, p, t)
    x, y, z, vx, vy, vz = u;
    wr, wz, lambda = p;
    
    du[1] = vx;
    du[2] = vy;
    du[3] = vz;
    du[4] = lambda * x * z - wr^2 * x;
    du[5] = lambda * y * z - wr^2 * y;
    du[6] = 1/2 * lambda * (x^2 + y^2) - wz^2 * z;
    nothing
end

function euler_lagrange_u_dyn(du, v, u, p, t)
    du .= v
    nothing
end

function euler_lagrange_v_dyn(dv, v, u, p, t)
    x, y, z = u;
    wr, wz, lambda = p;
    
    dv[1] = lambda * x * z - wr^2 * x;
    dv[2] = lambda * y * z - wr^2 * y;
    dv[3] = 1/2 * lambda * (x^2 + y^2) - wz^2 * z;
    nothing
end

function averaged_lagrangian_ode!(du, u, p, t)
    a, b, c, va, vb, vc = u;
    wr, wz, lambda = p;
    kappa = lambda / (4 * wr);
    
    du[1] = va;
    du[2] = vb;
    du[3] = vc;
    du[4] = -1im * kappa * conj(a) * c;
    du[5] = -1im * kappa * conj(b) * c;
    du[6] = 1/4 * kappa * (a^2 + b^2);
    nothing
end

function averaged_lagrangian_u_dyn(du, v, u, p, t)
    du .= v
    nothing
end

function averaged_lagrangian_v_dyn(dv, v, u, p, t)
    a, b, c = u;
    wr, wz, lambda = p;
    kappa = lambda / (4 * wr);
    
    dv[1] = -1im * kappa * conj(a) * c;
    dv[2] = -1im * kappa * conj(b) * c;
    dv[3] = -1im * 1/4 * kappa * (a^2 + b^2);
    nothing
end

function avglag2xyz(sol, wr)
    sol[1, :] = real.(sol[1, :] .* exp.(1im .* wr * sol.t))
    sol[2, :] = real.(sol[2, :] .* exp.(1im .* wr * sol.t))
    sol[3, :] = real.(sol[3, :] .* exp.(2im .* wr * sol.t))
    return real(sol)
end

function hamiltonian(sol, wr, wz, lambda)
    T = sol[4, :].^2 + sol[5, :].^2 + sol[6, :].^2;
    U = wr^2 * (sol[1, :].^2 + sol[2, :].^2) + wz^2 * sol[3, :].^2 - lambda * (sol[1, :].^2 + sol[2, :].^2) .* sol[3, :];
    return 1/2 .* (T .+ U)
end