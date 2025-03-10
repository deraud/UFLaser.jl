include("./constants.jl")

module Structs

using ..Constants

export Lattice
export ThermalProperties

Base.@kwdef mutable struct Lattice
    # Electron and Lattice Temperature (1D Vector)
    Tₑ::Vector{Float64} = 300 * ones(Nx)  
    Tₗ::Vector{Float64} = 300 * ones(Nx)  

    # Distribution Functions (2D Matrices)
    fₑ::Matrix{Float64} = zeros(Nx, Q) 
    fₗ::Matrix{Float64} = zeros(Nx, Q)

    fprop::Matrix{Float64} = zeros(Nx, Q) 
    fpropl::Matrix{Float64} = zeros(Nx, Q) 

    fₑeq::Matrix{Float64} = zeros(Nx, Q) 
    fₗeq::Matrix{Float64} = zeros(Nx, Q)

    # Scalar Fields (1D Vector Instead of Nx×1 Matrix)
    ϕₑ::Vector{Float64} = zeros(Nx)
    ϕₗ::Vector{Float64} = zeros(Nx)

    Source::Matrix{Float64} = zeros(Nx, Q)
    Coupling::Matrix{Float64} = zeros(Nx, Q)

    S::Vector{Float64} = zeros(Nx)
    x::Vector{Float64} = collect(0:dx:Lx)  
    G::Vector{Float64} = zeros(Nx)

    Ablation::Vector{Float64} = ones(Nx)  
    AblationEdge::Vector{Int} = zeros(Nx)

    UpperBoundary::Int64 = Nx
end


Base.@kwdef mutable struct ThermalProperties
    Cₑ::Array{Float64,1} = zeros(Nx)
    kₑ::Array{Float64,1} = zeros(Nx)
    υₑ::Array{Float64,1} = zeros(Nx)
    υₗ::Array{Float64,1} = zeros(Nx)
    αₑ::Array{Float64,1} = zeros(Nx)
    αₙ::Array{Float64,1} = zeros(Nx)

    Cₗ::Array{Float64,1} = zeros(Nx)
    kₗ::Array{Float64,1} = zeros(Nx)
    Cₚₗ::Array{Float64,1} = zeros(Nx)
    αₗ::Array{Float64,1} = zeros(Nx)
    αₙₗ::Array{Float64,1} = zeros(Nx)
    keq::Array{Float64,1} = zeros(Nx)
end

end