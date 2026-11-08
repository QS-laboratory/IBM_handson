using Serialization
using LinearAlgebra
BLAS.set_num_threads(1)
include(joinpath(@__DIR__, "..", "runtime", "Gate_generation_from_MPS_MPO.jl"))

function state_vector(psi)
    physical = [only(setdiff(inds(psi[j]),
                   vcat(j>1 ? collect(commoninds(psi[j],psi[j-1])) : Index[],
                        j<length(psi) ? collect(commoninds(psi[j],psi[j+1])) : Index[]))) for j in 1:length(psi)]
    return vec(Array(contract(psi), physical...))
end

function run_case(name, num_sweeps; dmrg_sweeps=200, depth=8)
    println("H4 $(geometry_id) $(name); depth=$(depth), sweeps=$(num_sweeps), dmrg_sweeps=$(dmrg_sweeps), Julia $(VERSION)")
    flush(stdout)
    # Upstream implementation is unchanged; increase DMRG sweeps for convergence.
    cache=joinpath(save_dir_pickle,"$(name)_dmrg$(dmrg_sweeps)_mps.dat")
    if isfile(cache)
        psi=deserialize(cache)
    else
        psi=get_psi_sp(name, dmrg_sweeps, vcat(fill(10,3),fill(50,12),fill(1000,5)), 2, 1e-8)
        serialize(cache,psi)
    end
    vector=state_vector(psi)
    source=np.load(joinpath(save_dir_pickle,"$(name).npz"))
    h=convert(Matrix{ComplexF64},source.__getitem__("hamiltonian"))
    spectrum=convert(Vector{Float64},source.__getitem__("full_eigenvalues"))
    ground=vector[1:2:end]*sqrt(2)
    excited=vector[2:2:end]*sqrt(2)
    eg=real(dot(ground,h*ground)); ee=real(dot(excited,h*excited))
    println("DMRG_AUDIT ground=$(eg) excited=$(ee) exact_ground=$(spectrum[1]) exact_excited=$(spectrum[2])")
    flush(stdout)
    @assert abs(eg-spectrum[1]) < 1e-7 "Ground MPS not converged; do not compress"
    @assert abs(ee-spectrum[2]) < 1e-7 "Excited MPS not converged; do not compress"
    serialize(joinpath(save_dir_pickle,"$(name)_mps.dat"),psi)
    np.save(joinpath(save_dir_npy,"$(name)_target_mps_vector.npy"),vector)
    prep, fidelities = get_and_save_prep_MPSBW_2d_list(psi, depth, num_sweeps, name; seed=1)
    np.save(joinpath(save_dir_npy,"$(name)_prep_history.npy"),fidelities)
    ham = hamcall(name)
    ref = get_gate_second_trotter(ham,8,0.1,100)
    serialize(joinpath(save_dir_pickle,"$(name)_evolution_mpo.dat"),ref)
    left,right = get_left_right_list(ref)
    ref_dense = reshape(Array(contract(ref),right...,left...),256,256)
    np.save(joinpath(save_dir_npy,"$(name)_reference_mpo_dense.npy"),ref_dense)
    evol, errors = get_and_save_evol_BW_2d_list(ref,depth,num_sweeps,name; seed=1)
    np.save(joinpath(save_dir_npy,"$(name)_evol_history.npy"),errors)
    println("COMPLETED $(name)");flush(stdout)
end

sweeps=parse(Int,get(ENV,"H4_OPT_SWEEPS","1000"))
for name in (isempty(ARGS) ? ["canonical","boys_ordered","boys"] : ARGS)
    run_case(name,sweeps)
end
