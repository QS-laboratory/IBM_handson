from pathlib import Path
import sys,json,hashlib
import numpy as np
import pytest
from qiskit.quantum_info import Statevector,Operator
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from author_analysis import load_meta,load_arrays,brickwall,qpde_circuit,evaluate,gaussian


@pytest.mark.parametrize('geometry',range(5))
def test_orbital_invariance_and_table_s1(geometry):
    expected=[.0159,.0793,.1283,.0769,.0171][geometry]
    spectra=[]
    for name in ['canonical','boys','boys_ordered']:
        m=load_meta(name,geometry);a=load_arrays(name,geometry)
        assert m['overlap_orthonormality_error']<1e-10
        assert abs(m['gap_exact_hartree']-expected)<.00005
        assert np.allclose(m['global_lowest_states_particle_numbers'],[4,4])
        assert np.isclose(m['triplet_s_squared'],2)
        spectra.append(a['full_eigenvalues'])
        np.testing.assert_allclose(np.linalg.eigvalsh(a['hamiltonian']),a['full_eigenvalues'],atol=1e-10)
    np.testing.assert_allclose(spectra,np.tile(spectra[0],(3,1)),atol=1e-9)
    m=load_meta('boys_ordered',geometry)
    assert abs(m['exchange_cost']-m['ga_exhaustive_min_cost'])<1e-10


def test_upstream_code_and_existing_materials_unchanged():
    assert (ROOT/'runtime/Gate_generation_from_MPS_MPO.jl').read_bytes()==(ROOT/'upstream/src/Gate_generation_from_MPS_MPO.jl').read_bytes()
    old=json.loads((ROOT/'existing_materials_sha256.json').read_text())
    revision=ROOT/'audience_revision_sha256.json'
    if revision.exists():old.update(json.loads(revision.read_text())['sha256'])
    for name,digest in old.items():
        if name in {'README.md','.gitignore'}:continue  # Repository publication metadata.
        assert hashlib.sha256((ROOT.parent/name).read_bytes()).hexdigest()==digest


@pytest.mark.parametrize('name',['canonical','boys','boys_ordered'])
def test_author_export_and_reference(name):
    m=evaluate(name)
    assert m['optimization_sweeps']==1000
    assert m['dmrg_ground_residual']<1e-7
    assert m['dmrg_excited_residual']<1e-7
    np.testing.assert_allclose(m['dmrg_particle_numbers'],[4,4],atol=1e-6)
    np.testing.assert_allclose(m['dmrg_s_squared'],[0,2],atol=1e-6)
    assert m['reference_mpo_adjoint_error_vs_matrix_exponential']<1e-3
    assert m['dense_second_trotter_error_vs_exponential']<1e-6
    assert m['reference_mpo_vs_dense_trotter_error']<1e-3
    assert abs(m['preparation_overlap']-m['author_preparation_final_cost'])<1e-6
    assert abs(m['evolution_delta_recomputed']-m['author_evolution_final_delta'])<1e-6
    for kind,expected in [('prep',32),('evol',28)]:
        qc=brickwall(name,kind)
        assert len(qc.data)==expected
        for gate in qc.data:
            sites=[qc.find_bit(q).index for q in gate.qubits]
            assert sites[1]-sites[0]==1


@pytest.mark.parametrize('name',['canonical','boys_ordered'])
def test_qpde_all_zero_amplitude(name):
    v=Statevector.from_instruction(brickwall(name,'prep')).data
    u=Operator(brickwall(name,'evol')).data
    u=np.linalg.matrix_power(u,2)
    a,b=np.vdot(v[::2],u@v[::2]),np.vdot(v[1::2],u@v[1::2])
    for eps in [0.,.1,.5]:
        expected=abs(a+np.exp(1j*eps*.2)*b)**2
        actual=abs(Statevector.from_instruction(qpde_circuit(name,eps,2)).data[0])**2
        assert abs(expected-actual)<1e-10


def test_gaussian_is_direct_author_function():
    x=np.array([-.3,0.,.5])
    np.testing.assert_allclose(gaussian(x,.7,.1,.2),.7*np.exp(-(x-.1)**2/(2*.2**2)))


def test_ideal_reference_control_reaches_exact_gap():
    p=ROOT/'results/geometry_0/canonical_qpde_reference_control_exact.json'
    control=json.loads(p.read_text())
    assert control['reference_control']
    assert control['status']=='variance_threshold_reached'
    assert abs(control['iterations'][-1]['posterior_mean']-control['exact_gap'])<1e-5
