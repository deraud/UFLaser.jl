module lbElectron

using ..Structs
using ..Constants
using ..Thermal
using ..Coupling

export CollisionElectronLB!
export StreamingElectron!
export UpdateElectron!
export BoundaryElectron!

function CollisionElectronLB!(lb::Lattice, t::ThermalProperties)
    lb.fₑ = collisionElectron(lb, t) + calculateLbSource(lb, t) - calculateLbCoupling(lb, t)
end

function collisionElectron(lb::Lattice, t::ThermalProperties)
    calculateϕₑ!(lb, t)
    Threads.@threads for k in 1:Q
        lb.fₑeq[:,k] = w[k]*lb.Tₑ[:]
        for i in 1:Nx
            lb.fprop[i,k] = lb.fₑ[i,k] - 2*lb.ϕₑ[i]*(lb.fₑ[i,k]-lb.fₑeq[i,k])
        end
    end
    return lb.fprop
end

function calculateLbSource(lb::Lattice, t::ThermalProperties)
    Threads.@threads for k in 1:Q
        for i in 1:Nx
            lb.Source[i,k] = w[k].*dt.*(lb.S[i])/t.Cₑ[i];
        end
    end
    return lb.Source
end

function calculateLbCoupling(lb::Lattice, t::ThermalProperties)
    calculateG!(lb)
    Threads.@threads for k in 1:Q
        for i in 1:Nx
            lb.Coupling[i,k] = dt*lb.G[i]*(lb.fₑ[i,k]-lb.fₗ[i,k])/t.Cₑ[i]
        end
    end
    return lb.Coupling
end

function StreamingElectron!(lb::Lattice)
    for k in 1:Q
        lb.fₑ[:,k] = circshift(lb.fₑ[:,k], [cx[k]])
    end
end

function UpdateElectron!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        lb.Tₑ[i] = 0
        for k in 1:Q
            lb.Tₑ[i] += lb.fₑ[i, k]
        end
    end
end

function BoundaryElectron!(lb::Lattice)
    # Apply boundary conditions (Adiabatic: zero temperature gradient)
    lb.Tₑ[1] = lb.Tₑ[2]           # Left boundary (copy from adjacent node)
    lb.Tₑ[end] = lb.Tₑ[end-1]     # Right boundary
    
    for k in 1:Q
        lb.fₑ[1, k] = w[k] * lb.Tₑ[1]       # Left boundary
        lb.fₑ[end, k] = w[k] * lb.Tₑ[end]   # Right boundary
    end
end

end