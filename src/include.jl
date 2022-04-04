using PyPlot
using LinearAlgebra
using LoopVectorization
using DifferentialEquations
using ArgParse
using HDF5
using StructArrays

include("plots.jl")
include("loop.jl")
include("frg.jl")
include("integrators.jl")

main()