##

using DifferentialEquations, LinearAlgebra, Printf, StaticArrays, CPUTime;


function cr3bpeom(X, (mu), t)
    @inbounds begin
        r1 = sqrt((X[1] + mu)^2 + X[2]^2)
        r2 = sqrt((X[1] - 1.0 + mu)^2 + X[2]^2)

        # dX[1:3] .= @view X[4:6] 
        # dX[4] = X[1] + 2*X[5] - (1 - mu)*(X[1] + mu)/r1^3 - mu*(X[1] - 1 + mu)/r2^3
        # dX[5] = X[2] - 2*X[4] - (1 - mu)*X[2]/r1^3 - mu*X[2]/r2^3
        # dX[6] = 0.0

        dX4 = X[1] + 2*X[5] - (1 - mu)*(X[1] + mu)/r1^3 - mu*(X[1] - 1 + mu)/r2^3
        dX5 = X[2] - 2*X[4] - (1 - mu)*X[2]/r1^3 - mu*X[2]/r2^3
    end
    SA_F64[X[4],X[5],X[6],dX4,dX5,0.0]
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
    R =  [sqrt(sum((X.^2)[:,i])) for i = 1:n]
    V = [sqrt(sum((Xdot.^2)[:,i])) for i = 1:n]

    # pseudo-potential
    Omega = 1/2*(R.^2) .+ (1-mu)./r1 .+ mu./r2
    
    # Jacobi's Constant
    C = (V.^2)/2 .- Omega
    
    return C
end;


# Event-1: detect a narrow square pulse
# Note that discrete callback wont be able to precisely locate the
# 0-crossings for this scenario
condition1 = function (u, t, integrator)
    if t > 0.5 && t < 0.5005
       1.0
    else
       -1.0
    end
end


# event times and solutions 
Te = Float64[]
#Xe = Array{Float64}(undef,0,2)

function affect1!(integrator)
    push!(Te, integrator.t)
end

# Event function-2: continuous callback
const Mask2 = Ref(false)
condition2 = function (u, t, integrator)
    if Mask2[] == false
        if u[1] < 0
            # TBD: how to safely ignore the event in this region?
            -1.0
        else
            u[2]
        end
    else
        -1.0
    end
end

function affect2!(integrator)
    Mask2[] = true
    push!(Te, integrator.t)
 end


# Event function-3: continuous callback
const Mask3 = Ref(false)
condition3 = function (u, t, integrator)
    if Mask3[] == false
        if u[1] > 0
            # TBD: how to safely ignore the event in this region?
            -1.0
        else
            u[2]
        end
    else
        -1.0
    end
end

function affect3!(integrator)
    # modify solution to set it back to IC
    integrator.u[1:3] .= @SVector [0.994::Float64, 0.0, 0.0]
    integrator.u[4:6] .= @SVector [0.0, parse(Float64,"-2.00158510637908252240537862224"), 0.0]
    Mask3[] = true
    push!(Te, integrator.t)
end


# test result Struct
struct ODETestResult
    FinalStateErr
    IOMErr
    MeanTime
    OdeSol
end

# Integrate ODE
function FireODE(atol, rtol, method, prob, DenseOn, cb)
#  if DenseOn == true
         #solve(prob,method,abstol=atol, reltol=rtol,timeseries_errors = false, dense_errors=false, dense=true,callback=cb)  
#         solve(prob,method,abstol=atol, reltol=rtol,timeseries_errors = false, dense_errors=false, dense=false)  
#     else 
        solve(prob,method,abstol=atol, reltol=rtol,timeseries_errors = false, dense_errors=false,callback=cb,dense=false)  
        # solve(prob,method,abstol=atol, reltol=rtol,timeseries_errors = false, dense_errors=false, dense=false)  
    # end
end
            
# Main ODE Integration Method Test function
function FireODEIntTest(X0, ODE, tspan, OdeMethod,  atol, rtol, FinalState, IOM, IOMFunc::Function, DenseOn, Runs, par)
    
    n = length(FinalState)

    # run solver bunch of times for benchmarking
    sol = 0
    rtime = 0.0
    X0new = copy(X0)

    cb1 = ContinuousCallback(condition1, nothing;
            rootfind=true,save_positions=(true,true),abstol=0.001)
    cb2 = ContinuousCallback(condition2, affect2!;
            rootfind=true,save_positions=(true,true),abstol=1e-9)
    cb3 = ContinuousCallback(condition3, affect3!;
            rootfind=true,save_positions=(true,true),abstol=1e-9)
    # callback set
    #cb = CallbackSet(cb1,cb2,cb3)
    cb = CallbackSet(cb2)

    # Define ODE problem
    prob = ODEProblem(ODE, X0new, tspan, par);                

    rtime = @elapsed for i in 1:Runs
        # randomsize initial conditions
        X0new[1] = X0new[1] + 0.000000000000*X0new[1].*rand()
        
        # Define ODE problem
        prob = ODEProblem(ODE, X0new, tspan, par);                
        
        #GC.gc()
        sol = FireODE(atol, rtol, OdeMethod, prob, DenseOn, cb)
    end


    # final state error metric
    serr = norm(sol[1:3,end] - FinalState[1:3])
    
    # compute IOM variations
    IOMsol = IOMFunc(sol, par)
    ierr = abs.(IOMsol .- IOM)
    
    ODETestResult(serr, ierr, rtime, sol)
end

##

struct Solver
    solver
    dense::Bool
    name::String
end


function main_CR3BP()

    print("FLINT benchmarks against Julia DiffEq package (CR3BP w/ events)\n")
    print("\n")

    Runs = 1000

    ###### @show rtol = 1.0./10.0.^(6:1:13);
    rtol = 1.0./10.0.^(11);
    atol = rtol*1e-3;

    mu = 0.012277471::Float64

    # Initial conditions
    R0 = [0.994::Float64, 0.0, 0.0]
    V0 = [0.0, parse(Float64,"-2.00158510637908252240537862224"), 0.0]

    # Time period
    Period = parse(Float64,"17.0652165601579625588917206249")

    # propagate for bunch of orbits
    t0 = 0.0::Float64;
    tf = 2.5*Period;
    tspan = [t0, tf];

    # X0 = [R0;V0]
    X0 = MVector{6,Float64}([0.994,0,0,0,parse(Float64,"-2.00158510637908252240537862224"),0.0])

    #op1 = ODEProblem(cr3bpeom,X0,tspan);

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


    nsol = length(ODESolvers)
    ntol = length(rtol)

    SolResults3B = Matrix(undef, nsol, ntol)

    for itr in 1:ntol
        for ctr in 1:nsol
            print(string("Running: ",ODESolvers[ctr].name, "\n"))
            SolResults3B[ctr,itr] = FireODEIntTest(X0, cr3bpeom, tspan, ODESolvers[ctr].solver,  atol[itr], rtol[itr], [R0; V0], 0, JacobiConstant, ODESolvers[ctr].dense,Runs,(mu));
        end
        #print(string("tol: ", rtol[itr],"\n"))
    end

    # Extract FLINT computation results from results.txt file

    cd(@__DIR__)
    DOP54_t, DOP853_t,Vern65E_t,Vern98R_t = open("results_CR3BP.txt") do f
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
    for ctr in 1:nsol
        sctime = SolResults3B[ctr,end].MeanTime/FLINT_res[ctr]
        @Printf.printf("%7s: %.3e [%.1f], \t\t %.3e (FLINT %7s)\n",
        ODESolvers[ctr].name,SolResults3B[ctr,end].MeanTime,sctime,FLINT_res[ctr],
        FLINT_alg[ctr])
    end


    # function test
    y0r = MVector{6,Float64}([R0[1],R0[2],R0[3],V0[1],V0[2],V0[3]])
    yfr = MVector{6,Float64}([0.0,0.0,0.0,0.0,0.0,0.0])
    dy = copy(yfr)
    t0 = 0.0
    floops = 10_000_000
    rtime = @CPUelapsed for i in 1:floops
        y0r = y0r + 0.0001*y0r
        dy = cr3bpeom(y0r, (mu), t0)
        yfr = yfr + dy
    end
    println(floops, " Func calls: ", rtime*1e3, " ms")


    #return SolResults3B



end


