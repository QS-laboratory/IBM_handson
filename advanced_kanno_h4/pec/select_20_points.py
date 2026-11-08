"""Select 20 existing PEC points; no new FCI or DMRG calculations."""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def main():
    data = np.genfromtxt(ROOT / 'linear_stretch.csv', delimiter=',', names=True)
    x = data['distance_angstrom']
    y = np.array([data['singlet_total_hartree'], data['triplet_total_hartree']])
    # Endpoints, sampled minima, workshop reference, and a long-distance anchor.
    anchors = [0.60, 0.90, 1.15, 2.00, 3.30, 4.00]
    chosen = {int(np.argmin(abs(x - r))) for r in anchors}
    while len(chosen) < 20:
        idx = sorted(chosen)
        residual = np.max(abs(PchipInterpolator(x[idx], y[:, idx], axis=1)(x) - y), axis=0)
        residual[idx] = -1
        chosen.add(int(residual.argmax()))
    idx = sorted(chosen)
    selected = x[idx]
    for source, target in [('linear_stretch.csv', 'linear_stretch_20.csv'),
                           ('linear_dmrg.csv', 'linear_dmrg_20.csv')]:
        with (ROOT / source).open(newline='') as stream:
            reader = csv.DictReader(stream)
            fields = reader.fieldnames
            rows = list(reader)
        rows = [row for row in rows if any(abs(float(row['distance_angstrom']) - r) < 1e-9 for r in selected)]
        assert len(rows) == 20
        with (ROOT / target).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    interpolator = PchipInterpolator(selected, y[:, idx], axis=1)
    error = np.max(abs(interpolator(x) - y), axis=1)
    (ROOT / 'selection_20.json').write_text(json.dumps({
        'distance_angstrom': selected.tolist(),
        'anchors_angstrom': anchors,
        'method': 'Greedy maximum absolute PCHIP residual across both states on the original 69-point grid',
        'max_residual_on_original_grid_hartree': dict(zip(['singlet', 'triplet'], error.tolist())),
        'note': 'Interpolation error against existing samples, not a rigorous continuous error bound. No new electronic-structure calculations.'
    }, indent=2) + '\n')
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(8.7, 5.5), layout='constrained')
    fine = np.linspace(x[0], x[-1], 800)
    curves = interpolator(fine)
    for k, (label, color) in enumerate([(r'$S_0$', '#2563eb'), (r'$T_1$', '#d95f02')]):
        ax.plot(x, y[k], color='#94a3b8', lw=3, alpha=.45,
                label='Original 69-point reference' if k == 0 else None)
        ax.plot(fine, curves[k], color=color, lw=1.8, label=label + ' PCHIP (20 points)')
        ax.plot(selected, y[k, idx], 'o', color=color, markerfacecolor='white', ms=5)
    reference = json.loads((ROOT / 'metadata.json').read_text())
    ax.axhline(reference['four_hydrogen_dissociation_limit_hartree'], ls='--',
               color='#64748b', lw=1, label='4 H limit (STO-3G)')
    ax.set(title='Linear H4 / STO-3G: 20 nonuniform sample points',
           xlabel=r'Adjacent H-H distance $R$ ($\AA$)', ylabel='Total energy (Hartree)')
    ax.grid(alpha=.16)
    ax.legend(frameon=False)
    fig.savefig(ROOT / 'h4_linear_pec_20.png', dpi=180)
    fig.savefig(ROOT / 'h4_linear_pec_20.svg')
    plt.close(fig)
    print('Selected distances (angstrom):', ', '.join(f'{r:.2f}' for r in selected))
    print('Maximum residual on original grid (Hartree):', error)


if __name__ == '__main__':
    main()
