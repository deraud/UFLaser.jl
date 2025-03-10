include("./constants.jl")

module Structs

using ..Constants

export Lattice
export Optics

Base.@kwdef mutable struct Lattice


    x::Array{Float64,1} = collect(-Lx/2:dx:Lx/2)
    y::Array{Float64,1} = collect(-Ly/2:dy:Ly/2)
    Tₑ::Array{Float64,2} = ones(Nx,Ny)*300

end

Base.@kwdef mutable struct Optics
    w::Array{Float64,1} = zeros(Nx)
    R::Array{Float64,1} = zeros(Nx)
    ς::Array{Float64,1} = zeros(Nx)
    I::Array{Float64,2} = zeros(Nx,Ny)
    E::Array{ComplexF64,2} = zeros(Nx,Ny)
    S::Array{Float64,2} = zeros(Nx,Ny)
    k::Array{Float64,2} = zeros(Nx,Ny)
    n::Array{Float64,2} = zeros(Nx,Ny)
    α::Array{Float64,2} = zeros(Nx,Ny)
end

end