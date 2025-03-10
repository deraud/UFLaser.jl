module LaserEnergy

using ..Structs
using ..Constants

export calculateS!

function calculateS!(lb::Lattice, op::Optics, t)
    calculateAbsorption!(op)
    calculateI!(lb, op, t)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            op.S[i,j] = op.I[i,j]*op.α[i,j]
        end
    end
end

function calculateE!(lb::Lattice, op::Optics, t)
    GaussianBeamSpotSize!(lb, op)
    BeamWaveCurvature!(lb, op)
    GouyPhaseShift!(lb, op)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            op.E[i,j] = (w0/op.w[i])*exp(-(t-2tₚ)^2/tₚ^2) *
            exp(-(lb.y[j])^2/op.w[i]^2 - (2π*sqrt(lb.x[i]^2)*op.k[i,j]/(λ)) - π*lb.y[j]^2*im/(op.R[i]/λ) +  op.ς[i]*im)
        end
    end

end

function GaussianBeamSpotSize!(lb::Lattice, op::Optics)
    Threads.@threads for i in 1:Nx
        op.w[i] = w0*sqrt(1 + (lb.x[i]/zᵣ)^2)
    end
end

function BeamWaveCurvature!(lb::Lattice, op::Optics)
    Threads.@threads for i in 1:Nx
        if lb.x[i] == 0
            op.R[i] = Inf 
        else
            op.R[i] = lb.x[i] * (1 + (zᵣ / lb.x[i])^2)
        end
    end
end

function GouyPhaseShift!(lb::Lattice, op::Optics)
    Threads.@threads for i in 1:Nx
        op.ς[i] = atan(lb.x[i]/zᵣ)
    end
end

function calculateI!(lb::Lattice, op::Optics, t)
    calculateE!(lb, op, t)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            op.I[i,j] = op.n[i,j]/2*sqrt(ϵ₀/μ₀)*abs(op.E[i,j])^2
        end
    end
end


function calculateAbsorption!(op::Optics)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            op.α[i,j] = 4π*op.k[i,j]/λ
        end
    end
end


end