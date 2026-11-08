"""Plot only the uniformly stretched linear H4 PECs."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})


def main():
    data=np.genfromtxt(ROOT/'linear_stretch.csv',delimiter=',',names=True)
    meta=json.loads((ROOT/'metadata.json').read_text())
    at_two=meta['linear_at_2_angstrom']
    fig,ax=plt.subplots(figsize=(8,5.2),layout='constrained')
    for label,title,color in [('singlet',r'$S_0$ (lowest singlet)','#2563eb'),('triplet',r'$T_1$ (lowest triplet)','#d95f02')]:
        ax.plot(data['distance_angstrom'],data[f'{label}_total_hartree'],color=color,lw=2.3,label=title)
        ax.scatter([2.],[at_two[f'{label}_total_hartree']],color=color,s=34,zorder=5)
    ax.axhline(meta['four_hydrogen_dissociation_limit_hartree'],color='#64748b',ls='--',lw=1.2,label='4 H dissociation limit (STO-3G)')
    ax.axvline(2.,color='#94a3b8',ls=':',lw=1.1)
    ax.set(xlabel=r'Adjacent H-H distance $R$ ($\AA$)',ylabel='Total energy (Hartree)',
           title='Linear H4: Full-CI / STO-3G\nAll 4 electrons and 4 orbitals; nuclear repulsion included')
    ax.text(.97,.04,r'$R=2\ \AA$: $E(T_1)-E(S_0)=0.4328$ eV',transform=ax.transAxes,ha='right',fontsize=10)
    ax.grid(alpha=.17)
    ax.legend(frameon=False,fontsize=10)
    fig.savefig(ROOT/'h4_linear_pec.png',dpi=180)
    fig.savefig(ROOT/'h4_linear_pec.svg')
    plt.close(fig)
    print('Saved linear H4 PEC (PNG/SVG)')

if __name__=='__main__':main()
