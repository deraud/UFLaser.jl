module Constants

using StaticArrays

export Nx, Ny, Lx, Ly, dx, dy, dt
export λ, c, ω#, α
export ϵ₀, μ₀
export h
export w0, tₚ, zᵣ
export T_itp, n_itp, k_itp
export tmax, Nt, timeToFemto


const Nx::Int = 100    # grid points in x
const Ny::Int = 1000    # grid points in y

const Lx::Float64 = 500e-9        # Length of the domain in x (m)
const Ly::Float64 = 5000e-9        # Length of the domain in y (m)

const dx::Float64 = Lx / Nx       # Spatial step size in x (m)
const dy::Float64 = Ly / Ny       # Spatial step size in y (m)
const dt::Float64 = 1e-15

const λ::Float64 = 800e-9
const c::Float64 = 3e8
const ω::Float64 = 2π * c / λ
#const α::Float64 = 4π*k/λ

const ϵ₀::Float64 = 8.854187817e-12  # permittivity of free space (F/m)
const μ₀::Float64 = 4π * 1e-7        # permeability of free space (H/m)

const h::Float64 = 4.135667696e-15 # planck constant in eV

const tₚ::Float64 = 100e-15
const w0::Float64 = 2500e-9
const zᵣ::Float64 = π*w0^2/λ

const T_itp = @SVector[300, 15000, 30000]
const n_itp = @SVector[3.58, 2.94, 2.55]
const k_itp = @SVector[4.62, 3.31, 2.91]

const tmax::Float64 = 400e-15     # Maximum simulation time (s)
const Nt::Int = Int(round(tmax / dt))    # Number of time steps
const timeToFemto::Array{Float64,1} = collect(0:dt:tmax)


end