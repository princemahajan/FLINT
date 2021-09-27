##

using DifferentialEquations, LinearAlgebra, Printf, StaticArrays, CPUTime;

function lorenz(u, p, t)
    @inbounds begin
        dx = p[1]*(u[2]-u[1])
        dy = u[1]*(p[2]-u[3]) - u[2]
        dz = u[1]*u[2] - p[3]*u[3]
    end
    SA_F64[dx,dy,dz]
end


function FireODEIntTest(X0, ODE, tspan, OdeMethod,  atol, rtol, DenseOn, Runs, par)
    
    # run solver bunch of times for benchmarking
    sol = 0
    rtime = 0.0
    X0new = copy(X0)

    # Define ODE problem
    prob = ODEProblem(ODE, X0new, tspan, par);                

    rtime = @elapsed for i in 1:Runs        
        # randomsize initial conditions
        X0new[1] = X0new[1] + 0.000000000001*X0new[1].*rand()
    
        # # Define ODE problem
        # prob = ODEProblem(ODE, X0new, tspan, par);                

        #GC.gc()
        sol = solve(prob, OdeMethod,reltol=rtol,abstol=atol, dense=DenseOn)
    end

    return rtime, sol
end


##

function main_Lorenz()


    print("FLINT benchmarks against Julia DiffEq package (Lorenz w/o events)\n")
    print("\n")


    u0 = MVector{3,Float64}([1.0,0.0,0.0])
    p = SA_F64[10.0,28.0,8/3]
    tspan = (0.0,100.0)
    rtol = 1.0e-11::Float64
    atol = 1.0e-14::Float64
    loops = 100

    # dense off
    rt_Tsit5, s_Tsit5 = FireODEIntTest(u0, lorenz, tspan, Tsit5(),  atol, rtol, false, loops, p)
    rt_DP5, s_DP5 = FireODEIntTest(u0, lorenz, tspan, DP5(),  atol, rtol, false, loops, p)
    rt_DP8, s_DP8 = FireODEIntTest(u0, lorenz, tspan, DP8(),  atol, rtol, false, loops, p)
    rt_Vern6, s_Vern6 = FireODEIntTest(u0, lorenz, tspan, Vern6(),  atol, rtol, false, loops, p)
    rt_Vern9, s_Vern9 = FireODEIntTest(u0, lorenz, tspan, Vern9(),  atol, rtol, false, loops, p)

    # #dense on
    # rt_Tsit5, s_Tsit5 = FireODEIntTest(u0, lorenz, tspan, Tsit5(),  atol, rtol, false, loops, p)
    # rt_DP5, s_DP5 = FireODEIntTest(u0, lorenz, tspan, DP5(),  atol, rtol, false, loops, p)
    # rt_DP8, s_DP8 = FireODEIntTest(u0, lorenz, tspan, DP8(),  atol, rtol, false, loops, p)
    # rt_Vern6, s_Vern6 = FireODEIntTest(u0, lorenz, tspan, Vern6(),  atol, rtol, false, loops, p)
    # rt_Vern9, s_Vern9 = FireODEIntTest(u0, lorenz, tspan, Vern9(),  atol, rtol, false, loops, p)

    allrt = [rt_Tsit5,rt_DP5,rt_DP8,rt_Vern6,rt_Vern9]
    solvers = ["Tsit5","DP5","DP8","Vern6","Vern9"]

    cd(@__DIR__)
    DOP54_t, DOP853_t,Vern65E_t,Vern98R_t = open("results_Lorenz.txt") do f
        dop54 = 0.0
        dop853 = 0.0
        vern65e = 0.0
        vern98r = 0.0
        for l in eachline(f)
            if occursin("DOP54", l) == true
                strs = split(l,' ',keepempty=false)
                dop54 = parse(Float32,strs[2])
            end
            if occursin("DOP853", l) == true
                strs = split(l,' ',keepempty=false)
                dop853 = parse(Float32,strs[2])
            end
            if occursin("Verner65E", l) == true
                strs = split(l,' ',keepempty=false)
                vern65e = parse(Float32,strs[2])
            end
            if occursin("Verner98R", l) == true
                strs = split(l,' ',keepempty=false)
                vern98r = parse(Float32,strs[2])
            end
            if occursin("Interpolated", l) == true
                break
            end
        end
        (dop54,dop853,vern65e,vern98r)
    end

    FLINT_res = [DOP54_t,DOP54_t,DOP853_t,Vern65E_t,Vern98R_t]
    FLINT_alg = ["DOP54","DOP54","DOP853","Vern65E","Vern98R"]

    print("\n")
    print(string("Solver: Julia time (s) [scaled time], FLINT time (s)\n",
                "=====================================================\n"))
    for ctr in 1:length(solvers)
        sctime = allrt[ctr]/FLINT_res[ctr]
        @Printf.printf("%7s: %.3e [%.1f], \t\t %.3e (FLINT %7s)\n",
        solvers[ctr],allrt[ctr],sctime,FLINT_res[ctr],
        FLINT_alg[ctr])
    end



    # function test
    t0 = 0.0::Float64
    y = Array{Float64}([0.0,0.0,0.0])
    y0 = MVector{3,Float64}(0.0,0.0,0.0)
    su0 = SA_F64[1.0,0.0,0.0]
    floops = 100_000_000
    rtime = @CPUelapsed for i in 1:floops
        #su0 = su0 + 0.0001*su0
        y = lorenz(su0,p,t0)
        y0 = y + y0
    end
    println(floops, " Func calls: ", rtime*1e3, " ms")


end
