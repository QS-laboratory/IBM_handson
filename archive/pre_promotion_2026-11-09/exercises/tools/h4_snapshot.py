"""Read the pinned teaching snapshot. Embedded in 04/05; no research imports."""
def snapshot_root():
    candidates=[]
    if os.environ.get('H4_DATA_DIR'):
        candidates.append(Path(os.environ['H4_DATA_DIR']).expanduser())
    else:
        for parent in [Path.cwd(),*Path.cwd().parents]:
            candidates.extend([parent/'data/h4_pec_v1',parent/'hands_on/exercises/data/h4_pec_v1'])
    for root in candidates:
        if (root/'manifest.json').is_file():
            if hashlib.sha256((root/'manifest.json').read_bytes()).hexdigest()!=SNAPSHOT_DIGEST:
                raise ValueError('Snapshot manifest differs from the pinned h4-pec-v1 manifest')
            return root
    raise FileNotFoundError('全20構造には教材の data/h4_pec_v1 が必要です。リポジトリを取得し、'
                            'H4_DATA_DIRをそのフォルダのパスに設定してください。既定R=2は内蔵データでも実行できます。')


def list_geometries():
    print('Geometry ID          R / Angstrom   FCI S0 total / Ha   FCI T0 total / Ha')
    for row in SNAPSHOT_INDEX['geometries']:
        print(f"{row['id']:20s} {row['R']:12.8f} {row['energy_S0']+row['enuc']:19.10f} {row['energy_T']+row['enuc']:19.10f}")


def snapshot_geometry(geometry):
    row=next((r for r in SNAPSHOT_INDEX['geometries'] if r['id']==geometry),None)
    if row is None:raise ValueError('Unknown GEOMETRY; choose an ID from list_geometries()')
    root=snapshot_root()
    def verified(relative):
        raw=(root/relative).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=SNAPSHOT_INDEX['files'][relative]:
            raise ValueError(f'Snapshot checksum mismatch: {relative}')
        return raw
    assets={k:base64.b64encode(verified(geometry+'/'+v)).decode() for k,v in row['assets'].items()}
    catalog=json.loads(verified(geometry+'/catalog.json'))
    for name in ['preparation','evolution']:
        catalog[name]={int(k) if k!='exact' else k:v for k,v in catalog[name].items()}
    provenance=json.loads(verified(geometry+'/provenance.json'))
    return assets,catalog,provenance


def select_geometry(geometry='R_2p00000000'):
    global ASSETS,CIRCUIT_CATALOG,PROVENANCE,MODEL,H,GAP,EXACT_STATES,U_EXACT,A_EXACT,SELECTED_GEOMETRY
    if geometry==SELECTED_GEOMETRY:return
    assets,catalog,provenance=snapshot_geometry(geometry)
    with np.load(io.BytesIO(base64.b64decode(assets['model'])),allow_pickle=False) as saved:
        model={k:saved[k].copy() for k in saved.files}
    ASSETS,CIRCUIT_CATALOG,PROVENANCE=assets,catalog,provenance
    MODEL=model;H=MODEL['hamiltonian'];GAP=float(MODEL['energy_T']-MODEL['energy_S0'])
    EXACT_STATES=[MODEL['S0'],MODEL['T']];U_EXACT=expm(-1j*H*DT);A_EXACT=prepare_pair(*EXACT_STATES)
    SELECTED_GEOMETRY=geometry
    select_circuits(SELECTED_PREP_M,SELECTED_EVOLUTION_M)
