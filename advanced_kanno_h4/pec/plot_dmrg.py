"""Overlay independently computed spin-adapted DMRG and Full-CI PECs."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})


def main():
    data=np.genfromtxt(ROOT/'linear_dmrg.csv',delimiter=',',names=True)
    meta=json.loads((ROOT/'dmrg_metadata.json').read_text())
    reference=json.loads((ROOT/'metadata.json').read_text())
    fig,(ax,error_ax)=plt.subplots(2,1,figsize=(8.7,7.1),layout='constrained',sharex=True,gridspec_kw={'height_ratios':[3.2,1]})
    colors={'singlet':'#2563eb','triplet':'#d95f02'}
    for label,name in [('singlet',r'$S_0$'),('triplet',r'$T_1$')]:
        color=colors[label]
        ax.plot(data['distance_angstrom'],data[f'{label}_fci_total_hartree'],color=color,lw=2.2,label=name+' Full-CI')
        ax.plot(data['distance_angstrom'][::3],data[f'{label}_dmrg_total_hartree'][::3],color=color,ls='none',marker='o',
                markerfacecolor='white',markersize=4.8,markeredgewidth=1.2,label=name+' DMRG')
        error_ax.plot(data['distance_angstrom'],data[f'{label}_dmrg_minus_fci_hartree']*1e14,color=color,lw=1.2,label=name)
    ax.axhline(reference['four_hydrogen_dissociation_limit_hartree'],ls='--',lw=1,color='#64748b',label='4 H limit (STO-3G)')
    ax.axvline(2.,color='#94a3b8',ls=':',lw=1)
    error_ax.axhline(0,color='#64748b',lw=1)
    ax.set(ylabel='Total energy (Hartree)',title='Linear H4 / STO-3G: Full-CI vs SU(2) DMRG\nAll 4 electrons and 4 orbitals; nuclear repulsion included')
    ax.legend(frameon=False,fontsize=10)
    ax.text(.97,.05,'DMRG: M = 16, two independent seeds\nMarkers overlap Full-CI curves',transform=ax.transAxes,ha='right',fontsize=10)
    error_ax.set(xlabel=r'Adjacent H-H distance $R$ ($\AA$)',ylabel=r'DMRG $-$ FCI'+'\n'+r'($10^{-14}$ Hartree)')
    error_ax.legend(frameon=False,ncol=2,fontsize=9,loc='upper right')
    for axes in [ax,error_ax]:axes.grid(alpha=.16)
    fig.savefig(ROOT/'h4_linear_fci_dmrg_pec.png',dpi=180)
    fig.savefig(ROOT/'h4_linear_fci_dmrg_pec.svg')
    plt.close(fig)
    print('Saved DMRG/Full-CI PEC and error comparison')

if __name__=='__main__':main()
