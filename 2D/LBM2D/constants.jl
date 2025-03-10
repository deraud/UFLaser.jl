module Constants

using StaticArrays

export Nx, Ny, scale, Lx, Ly, dx, dy, dt
export J, tₚ, R, δ, δᵦ
export ρ, Gᵣ, Aₑ, Bₗ
export N, kb, Be, Tf
export Q, w, cx, cy
export X, η, w0
export tmax, Nt, timeToFemto

Nx::Int = 100    # grid points in x
Ny::Int = 100    # grid points in y
const scale = 10
const Lx::Float64 = 100e-9        # Length of the domain in x (m)
const Ly::Float64 = 100e-9        # Length of the domain in y (m)
const dx::Float64 = Lx / Nx       # Spatial step size in x (m)
const dy::Float64 = Ly / Ny       # Spatial step size in y (m)
const dt::Float64 = 1e-15

const J::Float64 = 10
const tₚ::Float64 = 96e-15
const R::Float64 = 0.93
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

const w0 = 200

const X::Int = 353
const η::Float64 = 0.16

## D2Q9
const Q::Int64 = 9 # Number of directions
const w = @SVector[4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36] # weights
const cx = @SVector[0, 1, 0, -1, 0, 1, -1, -1, 1] # x-direction
const cy = @SVector[0, 0, 1, 0, -1, 1, 1, -1, -1] # y-direction


# Time and space arrays
const tmax::Float64 = 4000e-15     # Maximum simulation time (s)
const Nt::Int = Int(round(tmax / dt))    # Number of time steps
const timeToFemto::Array{Float64,1} = collect(0:dt:tmax)

end