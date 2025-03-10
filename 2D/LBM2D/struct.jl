include("./constants.jl")

module Structs

using ..Constants

export Lattice
export ThermalProperties

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

end