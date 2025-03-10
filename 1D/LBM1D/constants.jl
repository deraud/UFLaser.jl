module Constants

using StaticArrays

export Nx, scale, Lx, dx, dt
export J, tₚ, R, δ, δᵦ
export ρ, Gᵣ, Aₑ, Bₗ
export N, kb, Be, Tf
export Q, w, cx
export X, η
export tmax, Nt, timeToFemto

const Nx::Int = 100    # grid points in x
const scale = 10
const Lx::Float64 = 100e-9        # Length of the domain in x (m)
const dx::Float64 = Lx / Nx       # Spatial step size in x (m)
const dt::Float64 = 1e-15

const J::Float64 = 10#*4.66
const tₚ::Float64 = 96e-15
const R::Float64 = 0.926
const δ::Float64 = 20.6e-9 # Penetration Depth
const δᵦ::Float64 = 0e-9 # Ballistic Range

const ρ::Float64 = 19.3e3
const Gᵣ::Float64 = 2.2e16
const Aₑ::Float64 = 1.2e7
const Bₗ::Float64 = 1.23e11

const N::Float64 = 16*10^12;
const kb::Float64 = 5.670374419e-8
const Be::Float64 = 70
const Tf::Float64 = 6.42e4

const X::Int = 353
const η::Float64 = 0.16

## D1Q3
const Q::Int64 = 2 # Number of directions
const w = @SVector[1/2, 1/2] # weights
const cx = @SVector[1, -1] # x-direction

# Time and space arrays
const tmax::Float64 = 4000e-15     # Maximum simulation time (s)
const Nt::Int = Int(round(tmax / dt))    # Number of time steps
const timeToFemto::Array{Float64,1} = collect(0:dt:tmax)

end