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
            t.Cₑ[i,j] = γ*lb.Tₑ[i,j]
        end
    end
end

function calculateKbe!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.kₑ[i,j] = 235
        end
    end
end

function calculateCl!(lb::Lattice, t::ThermalProperties)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            t.Cₗ[i,j] = 2.42e6        
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