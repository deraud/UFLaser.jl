include("./constants.jl")

module Structs

using ..Constants

export Lattice
export ThermalProperties
export LBMDomain
export CellType

module CellType
    const FLUID = 2^0          # 1
    const INTERFACE = 2^1      # 2
    const GAS = 2^2            # 4
    const OBSTACLE = 2^3       # 8
    const INLET = 2^4          # 16
    const OUTLET = 2^5         # 32
    
    const NO_FLUID_NEIGH = 2^6  # 64
    const NO_EMPTY_NEIGH = 2^7  # 128
    const NO_IFACE_NEIGH = 2^8  # 256
    
    const TO_FLUID = 2^9       # 512
    const TO_GAS = 2^10        # 1024
end

Base.@kwdef mutable struct Lattice
    Tₑ::Array{Float64,2} = 300*ones(Nx,Ny) # Electron Temperature
    Tₗ::Array{Float64,2} = 300*ones(Nx,Ny) # Lattice Temperature

    f::Array{Float64,3} = zeros(Nx,Ny,Q) 
    fprop::Array{Float64,3} = zeros(Nx,Ny,Q) 
    fpropl::Array{Float64,3} = zeros(Nx,Ny,Q) 
    feq::Array{Float64,3} = zeros(Nx,Ny,Q) 
    #f₀::Array{Float64,3} = zeros(Nx,Ny,Q) 
    #feq₀::Array{Float64,3} = zeros(Nx,Ny,Q) 
    ϕₑ::Array{Float64,2} = zeros(Nx,Ny)
    ϕₗ::Array{Float64,2} = zeros(Nx,Ny)
    Source::Array{Float64,3} = zeros(Nx,Ny,Q)
    Coupling::Array{Float64,3} = zeros(Nx,Ny,Q)

    fₑ::Array{Float64,3} = zeros(Nx,Ny,Q) 
    fₗ::Array{Float64,3} = zeros(Nx,Ny,Q)
    fₑeq::Array{Float64,3} = zeros(Nx,Ny,Q) 
    fₗeq::Array{Float64,3} = zeros(Nx,Ny,Q)
    #fₑ₀::Array{Float64,3} = zeros(Nx,Ny,Q) 
    #fₗ₀::Array{Float64,3} = zeros(Nx,Ny,Q)
    #fₑeq₀::Array{Float64,} = zeros(Nx,Ny,Q) 
    #fₗeq₀::Array{Float64,3} = zeros(Nx,Ny,Q)
    S::Array{Float64,2} = zeros(Nx,Ny)
    x::Array{Float64,1} = collect(0:dx:Lx)
    y::Array{Float64,1} = collect(0:dy:Ly)
    G::Array{Float64,2} = zeros(Nx,Ny)

    Ablation::Array{Float64,2} = ones(Nx,Ny)
    AblationEdge::Array{Int,1} = zeros(Ny)
    AblationEdgeHorizontal::Array{Int,1} = ones(Nx)

    UpperBoundary::Int64 = Nx
end

Base.@kwdef mutable struct ThermalProperties
    Cₑ::Array{Float64,2} = zeros(Nx,Ny)
    kₑ::Array{Float64,2} = zeros(Nx,Ny)
    υₑ::Array{Float64,2} = zeros(Nx,Ny)
    υₗ::Array{Float64,2} = zeros(Nx,Ny)
    αₑ::Array{Float64,2} = zeros(Nx,Ny)
    αₙ::Array{Float64,2} = zeros(Nx,Ny)

    Cₗ::Array{Float64,2} = zeros(Nx,Ny)
    kₗ::Array{Float64,2} = zeros(Nx,Ny)
    Cₚₗ::Array{Float64,2} = zeros(Nx,Ny)
    αₗ::Array{Float64,2} = zeros(Nx,Ny)
    αₙₗ::Array{Float64,2} = zeros(Nx,Ny)
    keq::Array{Float64,2} = zeros(Nx,Ny)
end

@kwdef mutable struct LBMDomain
    nx::Int
    ny::Int
    omega::Float64 = 1.0
    gravity::Vector{Float64} = [0.0, -0.1]
    
    # Distribution functions and macroscopic quantities
    fin::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    fout::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    equi::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    fdist::Array{Float64, 3} = zeros(Float64, 9, nx, ny)
    inlet::Array{Float64, 3} = zeros(Float64, 2, nx, ny)
    u::Array{Float64, 3} = zeros(Float64, 2, nx, ny)
    rho::Matrix{Float64} = zeros(Float64, nx, ny)
    mass::Matrix{Float64} = zeros(Float64, nx, ny)
    cell_type::Matrix{Int} = fill(CellType.GAS, nx, ny)
    
    # Storage for previous timestep
    rho_prev::Matrix{Float64} = zeros(Float64, nx, ny)
    mass_prev::Matrix{Float64} = zeros(Float64, nx, ny)
    u_prev::Array{Float64, 3} = zeros(Float64, 2, nx, ny)
    cell_type_prev::Matrix{Int} = fill(CellType.GAS, nx, ny)
    
    # History for visualization
    velocity_history::Vector{Array{Float64, 3}} = Array{Float64, 3}[]
    cell_type_history::Vector{Matrix{Int}} = Matrix{Int}[]

    Te::Matrix{Float64} = zeros(Float64, nx, ny)  # Electron temperature
    Tl::Matrix{Float64} = zeros(Float64, nx, ny)  # Lattice temperature
end

end