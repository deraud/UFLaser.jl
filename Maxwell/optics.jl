module OpticsFunctions

using Interpolations

using ..Structs
using ..Constants

export calculateIndexK!
export calculateIndexN!

function calculateIndexN!(lb::Lattice, op::Optics)
    Threads.@threads for i in trunc(Int, Nx/2):Nx
        for j in 1:Ny
            op.n[i,j] = interpolate_value_n(lb.Tₑ[i,j])
        end
    end
    Threads.@threads for i in 1:trunc(Int, Nx/2)
        for j in 1:Ny
            op.n[i,j] = 1
        end
    end
end

function calculateIndexK!(lb::Lattice, op::Optics)
    Threads.@threads for i in trunc(Int, Nx/2):Nx
        for j in 1:Ny
            op.k[i,j] = interpolate_value_k(lb.Tₑ[i,j])
        end
    end
    Threads.@threads for i in 1:trunc(Int, Nx/2)
        for j in 1:Ny
            op.k[i,j] = 1e-9
        end
    end
end

interpK = LinearInterpolation(T_itp, k_itp, extrapolation_bc=Flat())
interpN = LinearInterpolation(T_itp, n_itp, extrapolation_bc=Flat())

function interpolate_value_k(x_val)
    if x_val < 300
        return k_itp[1]
    elseif x_val > 30000
        return k_itp[end]
    else
        return interpK(x_val)
    end
end

function interpolate_value_n(x_val)
    if x_val < 300
        return n_itp[1]
    elseif x_val > 30000
        return n_itp[end]
    else
        return interpN(x_val)
    end
end

end