module Coupling

using ..Structs
using ..Constants

using Einsum

export calculateG!

function calculateG!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            lb.G[i,j] = 5.69e17
        end
    end
end


end
