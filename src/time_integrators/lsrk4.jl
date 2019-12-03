module LSRK4
"""
 Low storage Runge-Kutta
 M.H. Carpenter and C.A. Kennedy. 
 Fourth-order 2N-storage Runge-Kutta schemes. 
 TechnicalReport NASA TM-109112,
 4th order 5-4 solution 3
 y = lsrk4(g,y,t,dt)
"""


# Butcher Tableaus 
a = [0
     -567301805773/1357537059087
     -2404267990393/2016746695238 
     -3550918686646/2091501179385
     -1275806237668/842570457699]
b = [1432997174477/9575080441755
    5161836677717/13612068292357
    1720146321549/2090206949498
    3134564353537/4481467310338
    2277821191437/14882151754819]
c = [0
    1432997174477/9575080441755
    2526269341429/6820363962896
    2006345519317/3224310063776
    2802321613138/2924317926251]


"""
 Input:
     g: g(y,t) function to integrate in time
     y: Solution vector at current time t
    dy: Temporary storage of the same size as y
     t: Current time
    dt: Time step

 Output:
     y: Solution at time t + dt
"""
function lsrk4(g::Function, 
               y::AbstractArray, 
               dy::AbstractArray, 
               t::Float64, dt::Float64)

        dy = 0 * dy
        
        for k in 1:size(a, 1)
        
            # Set rates
            @inbounds dy = a[k] .* dy + g(y, t + c[k] * dt)
        
            # Update fields
            @inbounds y = y + dt .* b[k] .* dy
        
        end  
        return y
end

end

