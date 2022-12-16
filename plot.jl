using Plots
using GLMakie
using DataStructures: CircularBuffer
using LaTeXStrings

function plot_animation(sol, fname)
    soln = mapreduce(permutedims, vcat, sol.u);
    x = soln[:, 1];
    y = soln[:, 2];
    z = soln[:, 3];
    t = sol.t;
    mins = minimum(soln, dims=1)
    maxs = maximum(soln, dims=1)
    limits = vcat(mins[1:3], maxs[1:3]).*1.1

    trail_len = length(t) #1000
    # println("Trail length: $trail_len")

    points = CircularBuffer{Point2f}(trail_len)
    fill!(points, Point3f(x[1], y[1], z[1]))
    points = Observable(points)

    colors = CircularBuffer{Int}(trail_len)
    fill!(colors, 0)
    colors = Observable(colors)

    fig, ax, l = lines(points, color=colors, colormap=:greens, transparency=true,
        axis=(; type=Axis3, protrusions=(0, 0, 0, 0), viewmode=:fit, limits=(limits[1], limits[4], limits[2], limits[5], limits[3], limits[6])))

    updates_per_frame = 100;
    num_frames = fld(length(t), updates_per_frame)
    record(fig, fname, 0:(num_frames-1)) do frame
        for i in 1:updates_per_frame
            idx = frame * updates_per_frame + i
            pushfirst!(points[], Point3f(x[idx], y[idx], z[idx]))
            pushfirst!(colors[], frame)
        end
        ax.azimuth[] = π/4 + 0.3 * sin(3π * frame/num_frames)
        ax.elevation[] = 0.2π - 0.3 * sin(6π * frame/num_frames)
        notify.((points, colors))
        l.colorrange = (frame, 0)
    end
end

function plot_still(sol, color)
    xy = Plots.plot(sol, idxs=(1, 2), linealpha=0.3, xlab=L"x", ylab=L"y", dpi=150, legend=false, aspect_ratio=:equal, fontfamily="Computer Modern", color=color)
    xyz = Plots.plot(sol, idxs=(1, 2, 3), linealpha=0.3, xlab=L"x", ylab=L"y", zlab=L"z", camera=[45, 65], dpi=150, legend=false, fontfamily="Computer Modern", color=color)
    return xy, xyz
end

function plot_comparison(sol1, sol2, l1, l2, l3, c3)
    tx = Plots.plot(sol1.t, [sol1[1, :], sol2[1, :]], linealpha=0.1, label=[l1 l2], legend=:topleft, ylab=L"x", fontfamily="Computer Modern");
    tx_delta = twinx(tx);
    plt_tx = Plots.plot!(tx_delta, sol1.t, sol1[1, :] - sol2[1, :], linealpha=0.1, color=c3, label=l3, legend=:bottomleft, ytickfontcolor=c3, fontfamily="Computer Modern");

    ty = Plots.plot(sol1.t, [sol1[2, :], sol2[2, :]], linealpha=0.1, label=[l1 l2], legend=:topleft, ylab=L"y", fontfamily="Computer Modern");
    ty_delta = twinx(ty);
    plt_ty = Plots.plot!(ty_delta, sol1.t, sol1[2, :] - sol2[2, :], linealpha=0.1, color=c3, label=l3, legend=:bottomleft, ytickfontcolor=c3, fontfamily="Computer Modern");

    tz = Plots.plot(sol1.t, [sol1[3, :], sol2[3, :]], linealpha=0.1, label=[l1 l2], legend=:topleft, ylab=L"z", xlab=L"t", fontfamily="Computer Modern");
    tz_delta = twinx(tz);
    plt_tz = Plots.plot!(tz_delta, sol1.t, sol1[3, :] - sol2[3, :], linealpha=0.1, color=c3, label=l3, legend=:bottomleft, ytickfontcolor=c3, fontfamily="Computer Modern");

    plt_txyz = Plots.plot(plt_tx, plt_ty, plt_tz, layout=(3, 1));
    plt_xy = Plots.plot([sol1[1, :], sol2[1, :]], [sol1[2, :], sol2[2, :]], linealpha=0.1, xlab=L"x", ylab=L"y", label=[l1 l2], aspect_ratio=:equal, fontfamily="Computer Modern");
    plt_2d = Plots.plot(plt_txyz, plt_xy, layout=(1, 2), w=1)
    plt_3d = Plots.plot([sol1[1, :], sol2[1, :]], [sol1[2, :], sol2[2, :]], [sol1[3, :], sol2[3, :]], linealpha=0.1, xlab=L"x", ylab=L"y", zlabel=L"z", camera=[45, 65], label=[l1 l2], fontfamily="Computer Modern");
    return plt_2d, plt_3d
end