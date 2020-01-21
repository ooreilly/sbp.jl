module StaggeredRK4
"""
Staggered Time Integrators for Wave Equations
Michelle Ghrist, Bengt Fornberg, and Tobin A. Driscoll
https://doi.org/10.1137/S0036142999351777
"""

function staggered_rk4(
                f::Function, 
                g::Function, 
               u::AbstractArray, 
               v::AbstractArray, 
               t::Float64, dt::Float64)

        uout = step(f, g, u, v, t, dt)
        vout = step(g, f, v, uout, t + 0.5 * dt, dt)


        return uout, vout
        
end

function step(
                f::Function, 
                g::Function, 
               u::AbstractArray, 
               v::AbstractArray, 
               t::Float64, dt::Float64)

        z1 = dt * g(u, t)
        z1 = v - z1
        z1 = dt * f(z1, t + 0.5 * dt - dt)
        z2 = dt * f(v, t + 0.5 * dt)
        z3 = u + z2
        z3 = dt * g(z3, t + dt)
        z3 = v + z3
        z3 = dt * f(z3, t + 0.5 * dt + dt)
        u = u + z1 / 24 + 11 / 12 * z2  + z3 / 24 

        #z1 = dt * f(v - dt * g(u, t), t + 0.5 * dt - dt)
        #z2 = dt * f(v, t + 0.5 * dt)
        #z3 = dt * f(v + dt * g(u + z2, t + dt) , t + 0.5 * dt + dt)
        #u = u + z1 / 24 + 11 / 12 * z2  + z3 / 24 

        return u

end

function staggered_rk4_linear(
                              Ap, Av, p, v, dt::Float64)
        p = step_linear(Ap, Av, p, v, dt)
        v = step_linear(Av, Ap, v, p, dt)
        return p, v

end

function step_linear(Ap, Av, p, v, dt)
        r1 = dt * Ap * v
        p = p + r1
        r2 = dt * Av * r1
        r1 = dt * Ap * r2
        p = p + r1 / 24
        return p
end


end
