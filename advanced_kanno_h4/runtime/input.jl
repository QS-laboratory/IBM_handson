# Path-only adapter. The upstream Gate_generation_from_MPS_MPO.jl is byte-identical.
const H4_ROOT = normpath(joinpath(@__DIR__, ".."))
const geometry_id = parse(Int, get(ENV, "H4_GEOMETRY", "0"))
const save_dir_pickle = joinpath(H4_ROOT, "data", "geometry_$(geometry_id)")
const save_dir_npy = joinpath(H4_ROOT, "results", "geometry_$(geometry_id)")
const pickle_dir = save_dir_pickle
mkpath(save_dir_npy)
