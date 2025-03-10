module Coupling

using ..Structs
using ..Constants

using Einsum

export calculateG!

function calculateG!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            lb.G[i,j] = Gᵣ*(Aₑ/Bₗ*(lb.Tₑ[i,j]-lb.Tₗ[i,j])+1)
        end
    end
end


end
