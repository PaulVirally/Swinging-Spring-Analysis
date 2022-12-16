println("Importing...")
using DifferentialEquations
using LaTeXStrings
include("equations.jl")
include("plot.jl")
println("Done importing")

g = π^2;
l = 1;
wr = π;
wz = 2π;
lambda = l*g^2*wz^2/wr^4;

p = (wr, wz, lambda);
u0 = [0.006, 0.0, 0.012];
v0 = [0.0, 0.00489, 0.0];
tspan = (0.0, 1000.0);

reltol = 1e-5;
abstol = 1e-6;
num_steps = 50000; #1000000;
dt = (tspan[2] - tspan[1])/num_steps;

solver_ode = RK4(); # Fourth order solver
solver_dyn = CandyRoz4(); # Fourth order symplectic solver

function solve_ode_dyn(ode, solver_ode, dyn_v, dyn_u, solver_dyn, u0, v0, tspan, p, reltol, abstol, dt)
    prob_ode = ODEProblem(ode, vcat(u0, v0), tspan, p);
    prob_dyn = DynamicalODEProblem(dyn_v, dyn_u, v0, u0, tspan, p);

    sol_ode = solve(prob_ode, solver_ode, saveat=dt, reltol=reltol, abstol=abstol);
    sol_dyn = solve(prob_dyn, solver_dyn, saveat=dt, dt=dt, reltol=reltol, abstol=abstol);

    return sol_ode, sol_dyn
end

println("Solving the ODEs...")
sol_ode_el, sol_dyn_el = solve_ode_dyn(euler_lagrange_ode!, solver_ode, euler_lagrange_v_dyn, euler_lagrange_u_dyn, solver_dyn, u0, v0, tspan, p, reltol, abstol, dt);
# sol_ode_al, sol_dyn_al = solve_ode_dyn(averaged_lagrangian_ode!, solver_ode, averaged_lagrangian_v_dyn, averaged_lagrangian_u_dyn, solver_dyn, complex(u0), complex(v0), tspan, p, reltol, abstol, dt);

# The symplectic integrator swaps the velocities with the positions, so we need to swap them back
sol_dyn_el[[1, 2, 3, 4, 5, 6], :] = sol_dyn_el[[4, 5, 6, 1, 2, 3], :]
# sol_dyn_al[[1, 2, 3, 4, 5, 6], :] = sol_dyn_al[[4, 5, 6, 1, 2, 3], :]

# Convert from abc space to xyz space
# sol_ode_al = avglag2xyz(sol_dyn_al, wr);
# sol_dyn_al = avglag2xyz(sol_dyn_al, wr);
println("Done solving the ODEs")


println("Plotting...")
colors = palette(:default);

H_ode_el = hamiltonian(sol_ode_el, wr, wz, lambda);
H_dyn_el = hamiltonian(sol_dyn_el, wr, wz, lambda);
H_ode_el_plt = Plots.plot(sol_ode_el.t, H_ode_el, title="RK4", xlabel=L"t", ylabel="Hamiltonian", legend=false, color=colors[1], fontfamily="Computer Modern");
H_dyn_el_plt = Plots.plot(sol_dyn_el.t, H_dyn_el, title="CandyRoz4", xlabel=L"t", ylabel="Hamiltonian", legend=false, color=colors[2], fontfamily="Computer Modern");
H_el_plt = Plots.plot(H_ode_el_plt, H_dyn_el_plt, layout=(2,1), fontfamily="Computer Modern");
H_both_el_plt = Plots.plot(sol_ode_el.t, [H_ode_el, H_dyn_el], label=["RK4" "CandyRoz4"], xlabel=L"t", ylabel="Hamiltonian", fontfamily="Computer Modern");
H_plt = Plots.plot(H_both_el_plt, H_el_plt, layout=(1,2));
savefig(H_plt, "figures/hamiltonian.pdf")

ode_el_xy_plt, ode_el_xyz_plt = plot_still(sol_ode_el, colors[1])
dyn_el_xy_plt, dyn_el_xyz_plt = plot_still(sol_dyn_el, colors[2])
savefig(ode_el_xy_plt, "figures/rk4_el_xy.pdf")
savefig(ode_el_xyz_plt, "figures/rk4_el_xyz.pdf")
savefig(dyn_el_xy_plt, "figures/candy4roz4.pdf")
savefig(dyn_el_xyz_plt, "figures/candyroz4_el_xyz.pdf")

plt_2d, plt_3d = plot_comparison(sol_ode_el, sol_dyn_el, "RK4", "CandyRoz4", "Difference", colors[3]);
savefig(plt_2d, "figures/comparison_2d.pdf")
savefig(plt_3d, "figures/comparison_3d.pdf")

plot_animation(sol_ode_el, "figures/rk4.mp4")
plot_animation(sol_dyn_el, "figures/candyroz4.mp4")
println("Done plotting")