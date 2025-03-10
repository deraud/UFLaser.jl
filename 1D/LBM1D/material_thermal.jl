module Thermal

using ..Structs
using ..Constants

export calculateϕₑ!
export calculateϕₗ!

function calculateCe!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        Cen = Be*Tf/π^2 + ((3*N*kb/2 - Be.*Tf/π^2)/(Tf-Tf/π^2))*(lb.Tₑ[i]-Tf/π^2)
        if lb.Tₑ[i] < Tf/π^2
            t.Cₑ[i] = Be*lb.Tₑ[i]
        elseif Tf/pi^2 < lb.Tₑ[i] && lb.Tₑ[i] < 3*Tf/π^2
            t.Cₑ[i] = 2*Be*lb.Tₑ[i]/3 + Cen/3
        elseif 3*Tf/pi^2 < lb.Tₑ[i] && lb.Tₑ[i] < Tf
            t.Cₑ[i] = N*kb + Cen/3
        else
            t.Cₑ[i] = 3*N*kb/2
        end
    end
end

function calculateKbe!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        t.υₑ[i] = lb.Tₑ[i]/Tf
        t.υₗ[i] = lb.Tₗ[i]/Tf
        t.kₑ[i] = X * ((t.υₑ[i]^2 + 0.16).^(5/4)*(t.υₑ[i]^2 + 0.44)*t.υₑ[i]) / ((t.υₑ[i]^2 + 0.092)^(1/2)*(t.υₑ[i]^2 + η*t.υₗ[i]))
    end
end

function calculateCl!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        t.Cₚₗ[i] = 105.1 + 0.294*lb.Tₗ[i]- 8.731e-4*lb.Tₗ[i]^2 + 1.787e-6*lb.Tₗ[i]^3- 7.051e-10*lb.Tₗ[i]^4 + 1.538e-13*lb.Tₗ[i]^5
        t.Cₗ[i] = t.Cₚₗ[i]*ρ        
    end
end

function calculatekl!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        t.keq[i] = 320.973 - 0.0111lb.Tₗ[i] - 2.747e-5*lb.Tₗ[i]^2 - 4.048e-9*lb.Tₗ[i]^3; 
        t.kₗ[i] = 0.01*t.keq[i]
    end
end

function calculateα!(lb::Lattice, t::ThermalProperties)
    calculateCe!(lb, t)
    calculateKbe!(lb, t)
    Threads.@threads for i in 1:Nx
        t.αₑ[i] = t.kₑ[i]/t.Cₑ[i]
    end
    return t.αₑ
end

function calculateαₗ!(lb::Lattice, t::ThermalProperties)
    calculateCl!(lb, t)
    calculatekl!(lb, t)
    Threads.@threads for i in 1:Nx
        t.αₗ[i] = t.kₗ[i]/t.Cₗ[i]
    end
    return t.αₗ
end

function calculateϕₑ!(lb::Lattice, t::ThermalProperties)
    calculateα!(lb, t)
    Threads.@threads for i in 1:Nx
        t.αₙ[i] = t.αₑ[i]*dt/dx^2
        lb.ϕₑ[i] = 1/(2*t.αₙ[i] + 1)
    end
end

function calculateϕₗ!(lb::Lattice, t::ThermalProperties)
    calculateαₗ!(lb, t)
    Threads.@threads for i in 1:Nx
        t.αₙₗ[i] = t.αₗ[i]*dt/dx^2
        lb.ϕₗ[i] = 1/(2*t.αₙₗ[i] + 1)
    end
end

end