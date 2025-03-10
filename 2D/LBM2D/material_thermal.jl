module Thermal

using ..Structs
using ..Constants

export calculateCe!
export calculateKbe!
export calculateCl!
export calculatekl!
export calculateα!
export calculateϕₑ!
export calculateϕₗ!

using Einsum

function calculateCe!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            Cen = Be*Tf/π^2 + ((3*N*kb/2 - Be.*Tf/π^2)/(Tf-Tf/π^2))*(lb.Tₑ[i,j]-Tf/π^2)
            if lb.Tₑ[i,j] < Tf/π^2
                t.Cₑ[i,j] = Be*lb.Tₑ[i,j]
            elseif Tf/pi^2 < lb.Tₑ[i,j] && lb.Tₑ[i,j] < 3*Tf/π^2
                t.Cₑ[i,j] = 2*Be*lb.Tₑ[i,j]/3 + Cen/3
            elseif 3*Tf/pi^2 < lb.Tₑ[i,j] && lb.Tₑ[i,j] < Tf
                t.Cₑ[i,j] = N*kb + Cen/3
            else
                t.Cₑ[i,j] = 3*N*kb/2
            end
        end
    end
end

function calculateKbe!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.υₑ[i,j] = lb.Tₑ[i,j]/Tf
            t.υₗ[i,j] = lb.Tₗ[i,j]/Tf
            t.kₑ[i,j] = X * ((t.υₑ[i,j]^2 + 0.16).^(5/4)*(t.υₑ[i,j]^2 + 0.44)*t.υₑ[i,j]) / ((t.υₑ[i,j]^2 + 0.092)^(1/2)*(t.υₑ[i,j]^2 + η*t.υₗ[i,j]))
        end
    end
end

function calculateCl!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.Cₚₗ[i,j] = 105.1 + 0.294*lb.Tₗ[i,j]- 8.731e-4*lb.Tₗ[i,j]^2 + 1.787e-6*lb.Tₗ[i,j]^3- 7.051e-10*lb.Tₗ[i,j]^4 + 1.538e-13*lb.Tₗ[i,j]^5
            t.Cₗ[i,j] = t.Cₚₗ[i,j]*ρ        
        end
    end
end

function calculatekl!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.keq[i,j] = 320.973 - 0.0111lb.Tₗ[i,j] - 2.747e-5*lb.Tₗ[i,j]^2 - 4.048e-9*lb.Tₗ[i,j]^3; 
            t.kₗ[i,j] = 0.01*t.keq[i,j]
        end
    end
end

function calculateα!(lb::Lattice, t::ThermalProperties)
    calculateCe!(lb, t)
    calculateKbe!(lb, t)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.αₑ[i,j] = t.kₑ[i,j]/t.Cₑ[i,j]
        end
    end
    return t.αₑ
end

function calculateαₗ!(lb::Lattice, t::ThermalProperties)
    calculateCl!(lb, t)
    calculatekl!(lb, t)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.αₗ[i,j] = t.kₗ[i,j]/t.Cₗ[i,j]
        end
    end
    return t.αₗ
end

function calculateϕₑ!(lb::Lattice, t::ThermalProperties)
    calculateα!(lb, t)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.αₙ[i,j] = t.αₑ[i,j]*dt/dx^2
            lb.ϕₑ[i,j] = 1/(2*t.αₙ[i,j] + 1)
        end
    end
end

function calculateϕₗ!(lb::Lattice, t::ThermalProperties)
    calculateαₗ!(lb, t)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.αₙₗ[i,j] = t.αₗ[i,j]*dt/dx^2
            lb.ϕₗ[i,j] = 1/(2*t.αₙₗ[i,j] + 1)
        end
    end
end

end