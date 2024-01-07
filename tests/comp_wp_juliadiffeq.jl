## Julia DifferentialEquations ODE Integrators Work-Precision Script
# Author: Bharat Mahajan

using DifferentialEquations, LinearAlgebra, Printf, StaticArrays;
using DelimitedFiles, Plots;

function cr3bpeom(X, (mu), t)
    @inbounds begin
        r1 = sqrt((X[1] + mu)^2 + X[2]^2 + X[3]^2)
        r2 = sqrt((X[1] - 1.0 + mu)^2 + X[2]^2 + X[3]^2)

        # dX[1:3] .= @view X[4:6] 
        # dX[4] = X[1] + 2*X[5] - (1 - mu)*(X[1] + mu)/r1^3 - mu*(X[1] - 1 + mu)/r2^3
        # dX[5] = X[2] - 2*X[4] - (1 - mu)*X[2]/r1^3 - mu*X[2]/r2^3
        # dX[6] = 0.0

        dX4 = X[1] + 2*X[5] - (1 - mu)*(X[1] + mu)/r1^3 - mu*(X[1] - 1 + mu)/r2^3
        dX5 = X[2] - 2*X[4] - (1 - mu)*X[2]/r1^3 - mu*X[2]/r2^3
        dX6 = - (1 - mu)*X[3]/r1^3 - mu*X[3]/r2^3
    end
    SA_F64[X[4],X[5],X[6],dX4,dX5,dX6]
end


function JacobiConstant(odesol, (mu))
    
    # mass-ratio
    #mu = 0.012277471::Float64

    X = odesol[1:3,:]
    Xdot = odesol[4:6,:]

    ~,n = size(X)
    
    # distances
    R1 = copy(X)
    R2 = copy(X)
    R1[1,:] = R1[1,:] .+ mu
    R2[1,:] = R2[1,:] .- 1 .+ mu
    r1 = [norm(R1[:,i]) for i = 1:n]
    r2 = [norm(R2[:,i]) for i = 1:n]
    
    # radius and velocity nagnitudes
    R =  [sqrt(sum((X.^2)[1:2,i])) for i = 1:n]
    V = [sqrt(sum((Xdot.^2)[:,i])) for i = 1:n]

    # pseudo-potential
    Omega = 1/2*(R.^2) .+ (1-mu)./r1 .+ mu./r2
    
    # Jacobi's Constant
    C = (V.^2)/2 .- Omega
    
    return C
end;


# test result Struct
struct ODETestResult
    FinalStateErr
    IOMErr
    MeanTime
    Solver
end

# Integrate ODE
function FireODE(atol, rtol, method, prob, DenseOn)
    solve(prob,method,abstol=atol, reltol=rtol,timeseries_errors = false, dense_errors=false,dense=DenseOn)  
end
            
# Main ODE Integration Method Test function
function FireODEIntTest(X0, ODE, tspan, OdeMethod,  atol, rtol, 
                FinalState, IOM, IOMFunc::Function, DenseOn, Runs, par)
    
    # run solver bunch of times for benchmarking
    sol = 0
    rtime = 0.0
    X0new = copy(X0)

    # Define ODE problem
    prob = ODEProblem(ODE, X0new, tspan, par);                

    rtime = @elapsed for i in 1:Runs        
        #GC.gc()
        sol = FireODE(atol, rtol, OdeMethod.solver, prob, DenseOn)
    end

    # error metrics
    serr = norm(sol[1:3,end] - FinalState[1:3])
    
    # compute IOM variations
    IOMsol = IOMFunc(sol, par)
    ierr = maximum(abs.(IOMsol .- IOM))
    
    ODETestResult(serr, ierr, rtime, OdeMethod.name)
end

##

struct Solver
    solver
    dense::Bool
    name::String
end


function main_wp_diffeq()

    Runs = 100

    ###### @show rtol = 1.0./10.0.^(6:1:13);
    rtol = 1.0./10.0.^(1:15);
    atol = 1e-80;

    mu::Float64 = 1.0::Float64/(81.30059::Float64 + 1.0)

    # Initial conditions
    R0 = [parse(Float64,"1.02202151273581740824714855590570360"), 0.0, 
            parse(Float64,"-0.182096761524240501132977765539282777")]
    V0 = [0.0, parse(Float64,"-0.103256341062793815791764364248006121"), 0.0]

    # Time period
    Period = parse(Float64,"1.5111111111111111111111111111111111111111")

    # propagate for bunch of orbits
    tf = Period;
    tspan = [0.0, tf];

    # X0 = [R0;V0]
    X0 = SVector{6,Float64}(vcat(R0, V0))

    IOM0 = JacobiConstant(X0, mu)

    # The 2nd element of each solver sets the dense=true option
    # The 3rd element is the name used in results

    ODESolvers = [
                    #Solver(dopri5(),false,"DOPRI5"),                              # Hairer's DOPRI5
                    Solver(DP5(),true,"DP5"),                                     # Julia's implementation of DOPRI5
                    Solver(Tsit5(),true,"TSIT5"),                                 # Tsitouras 5/4 RK
                    #Solver(dop853(),false,"DOP853"),                              # Hairer's DOP853
                    Solver(DP8(),true,"DP8"),                                     # Julia's implementation of DOP853
                    #Solver(ARKODE(Sundials.Explicit(), order = 8),true,"ARKODE8"),# Sundial's explicit RK
                    Solver(Vern6(),true,"Vern6"),                                 # Verner's 7/6 RK
                    #Solver(Vern7(),true,"Vern7"),                                 # Verner's 7/6 RK
                    #Solver(Vern8(),true,"Vern8"),                                 # Verner's 7/6 RK    
                    Solver(Vern9(),true,"Vern9"),                                 # Verner's 9/8 RK
                    #Solver(ddeabm(),false,"DDEABM"),                              # Shampine's DDEABM
                    #Solver(VCABM(),true,"JABM"),      # forcing dt_min           # Julia's implementation of Shampine's DDEABM
                    #Solver(CVODE_Adams(max_order = 12),true,"SunAM12") # Too much error!   # Sundial's Adams-Moulton solver
                ]


    nsol = length(ODESolvers);
    ntol = length(rtol);

    SpRes = Matrix(undef, nsol, ntol);
    DnRes = Matrix(undef, nsol, ntol);

    for itr in 1:ntol
        print(string("\n Running Sparse with rtol: ", rtol[itr],"\n"));
        for ctr in 1:nsol
            print(string(ODESolvers[ctr].name, ", "));
            SpRes[ctr,itr] = FireODEIntTest(X0, cr3bpeom, tspan, ODESolvers[ctr],  
                atol, rtol[itr], [R0; V0], IOM0, JacobiConstant, false,Runs,(mu));
        end
    end

    for itr in 1:ntol
        print(string("\n Running Dense with rtol: ", rtol[itr],"\n"));
        for ctr in 1:nsol
            print(string(ODESolvers[ctr].name, ", "));
            DnRes[ctr,itr] = FireODEIntTest(X0, cr3bpeom, tspan, ODESolvers[ctr],  
                atol, rtol[itr], [R0; V0], IOM0, JacobiConstant, true,Runs,(mu));
        end
    end

    # Extract FLINT computation results from txt files

    cd(@__DIR__)
    FLINT_solvers = ["DOP54","DOP853","Verner65E","Verner98R"];

    FLINT_data = []
    for FSol in FLINT_solvers
        fname = "wp_" * FSol * ".txt";
        fsdata, hdr = readdlm(fname; header = true);
        fsdata = Float64.(fsdata[:,1:end-1]);
        fsdata = replace(fsdata, 0 => NaN);
        push!(FLINT_data, fsdata);
    end  

    # make work-precision Plots
    gr()

    lw = 2
    ms = 3
    ls1 = :solid
    ls2 = :dot

    # Sparse time - rtol plots
    plt1 = plot(FLINT_data[1][:,1], FLINT_data[1][:,3], xaxis=:log10, yaxis=:log10, 
                label=string("FLINT ",FLINT_solvers[1]), ls = ls1, legend=:top,
                lw=lw, ms=ms,markershape = :auto, minorgrid=true)

    for ctr = 2:size(FLINT_data)[1]
        plot!(FLINT_data[ctr][:,1], FLINT_data[ctr][:,3], 
                    label=string("FLINT ",FLINT_solvers[ctr]), ls=ls1,
                    ms=ms, markershape = :auto, lw=lw, legend=:top,)
    end

    # Julia DiffEq results
    for ctr in 1:nsol
        tsdata = Float64.([SpRes[ctr,i].MeanTime for i =1:ntol])
        plot!(rtol, tsdata, label=string("Julia ", SpRes[ctr,1].Solver), 
                    ls=ls2,ms=ms, markershape = :auto, lw=lw,legend=:top)
    end
   
    xlabel!("Rel Tol");
    ylabel!("Time (s)");
    title!("Sparse Output Work-Precision")



    # Sparse time - Error plots
    plt2 = plot(FLINT_data[1][:,4], FLINT_data[1][:,3], xaxis=:log10, yaxis=:log10, 
            label=string("FLINT ",FLINT_solvers[1]), 
            lw=lw, ms=ms,markershape = :auto, minorgrid=true)

    for ctr = 2:size(FLINT_data)[1]
        plot!(FLINT_data[ctr][:,4], FLINT_data[ctr][:,3], 
                    label=string("FLINT ",FLINT_solvers[ctr]),
                    ms=ms, markershape = :auto, lw=lw)
    end

    # Julia DiffEq results
    for ctr in 1:nsol
        tsdata1 = Float64.([SpRes[ctr,i].IOMErr for i =1:ntol])
        tsdata2 = Float64.([SpRes[ctr,i].MeanTime for i =1:ntol])
        plot!(tsdata1, tsdata2, label=string("Julia ", SpRes[ctr,1].Solver), 
                    ls=ls2,ms=ms, markershape = :auto, lw=lw,legend=:top)
    end


    xlabel!("Error");
    ylabel!("Time (s)");
    title!("Sparse Output Work-Precision")


    # Dense time - rtol plots
    plt3 = plot(FLINT_data[1][:,1], FLINT_data[1][:,5], xaxis=:log10, yaxis=:log10, 
                label=string("FLINT ",FLINT_solvers[1]),
                lw=lw, ms=ms,markershape = :auto, minorgrid=true)

    for ctr = 2:size(FLINT_data)[1]
        plot!(FLINT_data[ctr][:,1], FLINT_data[ctr][:,5], 
                    label=string("FLINT ",FLINT_solvers[ctr]),
                    ms=ms, markershape = :auto, lw=lw)
    end

    # Julia DiffEq results
    for ctr in 1:nsol
        tddata = Float64.([DnRes[ctr,i].MeanTime for i =1:ntol])
        plot!(rtol, tddata, label=string("Julia ", DnRes[ctr,1].Solver), 
                    ls=ls2,ms=ms, markershape = :auto, lw=lw,legend=:top)
    end


    xlabel!("Rel Tol");
    ylabel!("Time (s)");
    title!("Dense Output Work-Precision")



    # Dense time - Error plots
    plt4 = plot(FLINT_data[1][:,6], FLINT_data[1][:,5], xaxis=:log10, yaxis=:log10, 
            label=string("FLINT ",FLINT_solvers[1]), 
            lw=lw, ms=ms,markershape = :auto, minorgrid=true)

    for ctr = 2:size(FLINT_data)[1]
        plot!(FLINT_data[ctr][:,6], FLINT_data[ctr][:,5], 
                    label=string("FLINT ",FLINT_solvers[ctr]),
                    ms=ms, markershape = :auto, lw=lw)
    end

    # Julia DiffEq results
    for ctr in 1:nsol
        tddata1 = Float64.([DnRes[ctr,i].IOMErr for i =1:ntol])
        tddata2 = Float64.([DnRes[ctr,i].MeanTime for i =1:ntol])
        plot!(tddata1, tddata2, label=string("Julia ", DnRes[ctr,1].Solver), 
                    ls=ls2,ms=ms, markershape = :auto, lw=lw,legend=:top)
    end



    xlabel!("Error");
    ylabel!("Time (s)");
    title!("Dense Output Work-Precision")


    # for fdata in FLINT_data
    #     sctime = SpRes[ctr,end].MeanTime/FLINT_res[ctr]
    #     @Printf.printf("%7s: %.3e [%.1f], \t\t %.3e (FLINT %7s)\n",
    #     ODESolvers[ctr].name,SpRes[ctr,end].MeanTime,sctime,FLINT_res[ctr],
    #     FLINT_alg[ctr])
    # end

    return (plt1, plt2, plt3, plt4)



end


