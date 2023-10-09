#using Pkg
#Pkg.add(["Plots", "DifferentialEquations"])

using Plots, DifferentialEquations


function model!(du, u, p, t)
    n, γ, μ = p
    du[1] = -( n * γ + μ) * u[1]
    if n > 1
        for i in 2:Int(n)
            du[i] = n * γ * u[i-1] - (n * γ + μ) * u[i]
        end
    end
    du[Int(n)+1] = n * γ * u[Int(n)] - μ * u[Int(n)+1]
end

u0 = [1,0,0,0,0,0,0] # I1 I2 I3 ... R
tspan = (0.0, 30.0)
n, γ, μ = p = [6, 0.1, 0]

prob = ODEProblem(model!, u0, tspan, p)
sol = solve(prob) # Tsit5

plot(sol, xlabel="Time", ylabel="Proporion", legend=(0.9, 0.9)) # idxs = (0, 4) = Time + Recovered