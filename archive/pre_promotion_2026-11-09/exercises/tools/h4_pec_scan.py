"""No-shot BPE/BPDE-C peak scans; globals come from the teaching runtime."""
def fit_phase_scan(beta, values):
    design=np.c_[np.ones(len(beta)),np.cos(beta),np.sin(beta)]
    coeff=np.linalg.lstsq(design,values,rcond=None)[0]
    amplitude=float(np.hypot(coeff[1],coeff[2]))
    if amplitude<1e-10:raise ValueError('Flat signal: phase peak is undefined')
    return float(np.arctan2(coeff[2],coeff[1])),float(np.max(abs(design@coeff-values)))


def track_survival_phase(phases, weights, L):
    """Track from t=0; certify no hidden winding using a derivative bound.

    Fractional steps interpolate the principal log of the saved one-step U.
    They are simulation-only auxiliary information for branch selection.
    No FCI target energy or Hamiltonian expectation enters branch selection.
    """
    for power in range(2,11):
        step=2**power;previous=np.ones(2,complex);phase=np.zeros(2)
        coarse_previous=np.ones(2,complex);coarse_phase=np.zeros(2)
        minimum=np.ones(2);end=None
        for start in range(0,L*step+1,2048):
            ix=np.arange(start,min(start+2048,L*step+1))
            a=np.exp(1j*(ix/step)[:,None]*phases)@weights
            phase+=np.sum(np.angle(a*np.conj(np.vstack([previous,a[:-1]]))),axis=0)
            previous=a[-1];end=a[-1]
            coarse=a[ix%2==0]
            if len(coarse):
                coarse_phase+=np.sum(np.angle(coarse*np.conj(np.vstack([coarse_previous,coarse[:-1]]))),axis=0)
                coarse_previous=coarse[-1]
            minimum=np.minimum(minimum,abs(a).min(axis=0))
        derivative=np.abs(phases)@weights
        ratio=float(np.max(derivative/(step*np.maximum(minimum,1e-300))))
        agreement=float(np.max(abs(phase-coarse_phase)))
        if ratio<=.2 and agreement<1e-8:
            return end,phase,dict(certified=True,auxiliary_step_L=1/step,
                derivative_bound=derivative.tolist(),minimum_amplitude=minimum.tolist(),
                hD_over_minimum=ratio,fine_coarse_difference=agreement)
    raise ValueError('Phase branch not certified: shorten L or use a better preparation')


def run_pec_scan(prep_m=8,evolution_m=32,L=64,output_path='outputs/exercise4_pec_scan.json'):
    if not isinstance(L,(int,np.integer)) or isinstance(L,(bool,np.bool_)) or not 1<=L<=1024:
        raise ValueError('PEC_L must be an integer in 1..1024')
    snapshot_root()  # all 20 structures are required; no R=2 fallback for this scan
    selected=(SELECTED_GEOMETRY,SELECTED_PREP_M,SELECTED_EVOLUTION_M)
    rows=[];beta=np.arange(9)*2*np.pi/9
    try:
        select_circuits(prep_m,evolution_m)
        for geometry in SNAPSHOT_INDEX['geometries']:
            select_geometry(geometry['id'])
            phases,vectors=SPECTRA['both']
            weights=np.stack([abs(vectors.conj().T@psi)**2 for psi in APPROX_STATES],axis=1)
            a,unwrapped,certificate=track_survival_phase(phases,weights,L)
            # Scan the same auxiliary phase beta in the BPE and BPDE-C circuits.
            bpe=(1+np.real(np.exp(1j*beta)*a[0]))/2
            bpde=abs(a[0]+np.exp(1j*beta)*a[1])**2/4
            peak0,res0=fit_phase_scan(beta,bpe);peakgap,resgap=fit_phase_scan(beta,bpde)
            anchors=np.array([-unwrapped[0],unwrapped[0]-unwrapped[1]])
            principal=np.array([peak0,peakgap])
            lifted=principal+2*np.pi*np.rint((anchors-principal)/(2*np.pi))
            np.testing.assert_allclose(lifted,anchors,atol=1e-8,rtol=0)
            energy0,gap=lifted/(L*DT)
            expected=[float(np.vdot(psi,H@psi).real) for psi in APPROX_STATES]
            enuc=float(MODEL['enuc'])
            rows.append(dict(geometry=geometry['id'],R=geometry['R'],prep_m=prep_m,evolution_m=evolution_m,L=int(L),
                S0_peak_total=float(energy0+enuc),T_peak_total=float(energy0+gap+enuc),gap_peak=float(gap),
                S0_FCI_total=float(MODEL['energy_S0']+enuc),T_FCI_total=float(MODEL['energy_T']+enuc),gap_FCI=GAP,
                S0_expect_total=expected[0]+enuc,T_expect_total=expected[1]+enuc,gap_expect=expected[1]-expected[0],
                beta=beta.tolist(),BPE_probability=bpe.tolist(),BPDE_C_probability=bpde.tolist(),
                principal_peak=principal.tolist(),unwrapped_peak=lifted.tolist(),
                max_fit_residual=max(res0,resgap),phase_certificate=certificate))
    finally:
        select_geometry(selected[0]);select_circuits(selected[1],selected[2])
    if output_path is not None:
        path=Path(output_path);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(rows,indent=2))
    return rows
