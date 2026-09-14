from .settings import *
from .titration import _Titration
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from uncertainties import ufloat, UFloat


class Boehm(_Titration):

    def __init__(self, boehm: str, data: dict, vmin: float, vmax: float, vtotal: float, valiquot: float, vhcl: float, chcl: UFloat, cnaoh: UFloat, 
                 cboehm: UFloat=None, mboehm: float=None, purityboehm: float=1, vboehm: float=None, mcarbon: float=None):

        super().__init__(data, vmin, vmax)

        self.vtotal: UFloat = ufloat(vtotal, vtotal/1e-3*U1000)
        self.valiquot: UFloat = ufloat(valiquot, valiquot/1e-3*U1000)
        self.vhcl: UFloat = ufloat(vhcl, vhcl/1e-3*U1000)

        self.chcl: UFloat = chcl
        self.cnaoh: UFloat =  cnaoh

        if cboehm:
            self.cboehm = cboehm
        elif mboehm:
            self.mboehm = ufloat(mboehm, UMASS)
            self.purityboehm = purityboehm
            self.vboehm = ufloat(vboehm, UVOLFLASK)
            self.cboehm = (self.mboehm*self.purityboehm/M[boehm])/self.vboehm

        if boehm == 'Na2CO3':
            mol = 2
        else:
            mol = 1
        if boehm == 'NaOH':
            specify = ' (Boehm)'
        else:
            specify = ''

        self.n = (self.valiquot*self.cboehm*mol + self.veq*self.cnaoh - self.vhcl*self.chcl)*(self.vtotal/self.valiquot)
        
        if mcarbon:
            self.mcarbon = ufloat(mcarbon, UMASS)
            self.nperg = self.n/self.mcarbon
            nperg = f'\n({self._numerr(self.nperg, prefix=1e-6)} μmol/g)'
            blank = ''
        else:
            self.mcarbon = 0
            self.nperg = None
            nperg = ''
            blank = ' (blank)'

        self.desc = f'{X[boehm]} Boehm base{blank} with NaOH after HCl acidification'

        self.info = (f'$V_{{\\mathrm{{total}}}}$ = {self._numerr(self.vtotal, prefix=1e-3)} mL\n'
                     f'$V_{{\\mathrm{{aliquot}}}}$ = {self._numerr(self.valiquot, prefix=1e-3)} mL\n'
                     f'$V_{{\\mathrm{{HCl}}}}$ = {self._numerr(self.vhcl, prefix=1e-3)} mL\n\n'
                     f'[{X[boehm]}{specify}] = {self._numerr(self.cboehm, sci=False)} M\n'
                     f'[{X['NaOH']}] = {self._numerr(self.cnaoh, sci=False)} M\n'
                     f'[{X['HCl']}] = {self._numerr(self.chcl, sci=False)} M\n\n'
                     f'$m_{{\\mathrm{{carbon}}}}$ = {self._numerr(self.mcarbon, prefix=1e-3)} mg{blank}\n'
                     f'$n_{{\\mathrm{{acidic}}}}$ = {self._numerr(self.n, prefix=1e-6)} μmol'
                     f'{nperg}')

    def __str__(self):
        string = self.info + '\n\n' + super().__str__()
        for substr in ['$', '{', '}', '\\mathrm', '_']:
            string = string.replace(substr, '')
        return string

    @staticmethod
    def result(nahco3: Boehm, na2co3: Boehm, naoh: Boehm, nahco3blank: Boehm=None, na2co3blank: Boehm=None, naohblank: Boehm=None) -> None:

        WIDTH = 0.35
        def othersetup(ax):
            fig.canvas.draw() # calculate ticks preliminarily
            ticks = ax.get_yticks()
            yscale = np.diff(ticks)[0]
            ax.set_ylim(ticks[0]-yscale*0.2, ticks[-1]+yscale*0.2)
            ax.set_yticks(np.arange(ticks[0], ticks[-1]+1e-10, yscale))
            ax.spines[:].set_visible(False)
            ax.tick_params(axis='both', length=0)
            ax.grid(axis='y', alpha=0.35)
            ax.axhline(y=0, linewidth=1, color='black')
            ax.set_axisbelow(True)
            ax.set_ylabel('$n_{\\mathrm{acidic}}$ / μmol g$^{-1}$')

        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        data = [(nahco3, nahco3blank), (na2co3, na2co3blank), (naoh, naohblank)]
        npergs = []

        for i, (boehm, blank) in enumerate(data):
            if blank:
                ax.bar(i-WIDTH*1.2/2, boehm.nperg.n*1e6, yerr=boehm.nperg.s*1e6, capsize=2, width=WIDTH, color='slategray')
                ax.bar(i+WIDTH*1.2/2, (blank.n/boehm.mcarbon).n*1e6, yerr=(blank.n/boehm.mcarbon).s*1e6, 
                       capsize=2, width=WIDTH, color='lightsteelblue')
                ax.text(i-WIDTH*1.2/2, -0.1, f'({Boehm._numerr(boehm.nperg, prefix=1e-6)})', ha="center", transform=ax.get_xaxis_transform())
                ax.text(i+WIDTH*1.2/2, -0.1, f'({Boehm._numerr(blank.n/boehm.mcarbon, prefix=1e-6)})', ha="center", transform=ax.get_xaxis_transform())
                npergs.append((boehm.n-blank.n)/boehm.mcarbon)
            else:
                ax.bar(i, boehm.nperg.n*1e6, yerr=boehm.nperg.s*1e6, capsize=2, width=WIDTH, color='slategray')
                ax.text(i, -0.1, f'({Boehm._numerr(boehm.nperg, prefix=1e-6)})', ha="center", transform=ax.get_xaxis_transform())
                npergs.append(boehm.nperg)

        ax.set_xticks(list(range(len(data))))
        ax.set_xticklabels([X['NaHCO3'], X['Na2CO3'], X['NaOH']])
        ax.set_title(f'Amount of acidic species reacted with each Boehm base')
        ax.legend(handles=[Patch(facecolor='slategray', label='Carbon'), Patch(facecolor='lightsteelblue', label='Blank')], 
                  loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
        othersetup(ax)
        plt.show()

        carboxylic = npergs[0] # able to react with nahco3
        lactonic = npergs[1] - npergs[0] # able to react with na2co3 but not nahco3
        phenolic = npergs[2] - npergs[1] # able to react with naoh but not na2co3 (or nahco3)

        results = [carboxylic, lactonic, phenolic]

        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        bars = ax.bar(np.arange(len(results)), [y.n*1e6 for y in results], yerr=[y.s*1e6 for y in results], 
                      capsize=2, color='skyblue')
        for i, y in enumerate(results):
                ax.text(i, -0.1, f'({Boehm._numerr(y, prefix=1e-6)})', ha="center", transform=ax.get_xaxis_transform())
        ax.set_xticks(list(range(len(results))))
        ax.set_xticklabels(['carboxylic', 'lactonic', 'phenolic'])
        ax.set_title(f'Amount of acidic surface functional groups')
        othersetup(ax)
        plt.show()

        return carboxylic, lactonic, phenolic