module Coupling

using ..Structs
using ..Constants

export calculateG!

function calculateG!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        lb.G[i] = Gᵣ*(Aₑ/Bₗ*(lb.Tₑ[i]+lb.Tₗ[i])+1)
    end
end


end