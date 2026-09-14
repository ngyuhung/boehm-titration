from .settings import *
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from uncertainties import ufloat, UFloat
from numpy.typing import NDArray

class _Titration:

    def __init__(self, data: dict, vmin: float, vmax: float):

        self.data = data
        self.vmin = vmin
        self.vmax = vmax

        self.desc = ''
        self.info = ''

        ph = []
        vol = []
        vprev = 0 # previous volume added

        for v, p in data.items():

            v = v*1e-6 # volumes recorded in uL

            ph.append(ufloat(p, UPH))

            if v == 0: # no volume added yet, normal float with zero uncertainty
                vol.append(0)
                continue

            if isinstance(vol[-1], UFloat):
                vprev = vol[-1].n
            else:
                vprev = vol[-1]

            vdiff = v - vprev
            if vdiff <= 10e-6: # 10 uL addition
                uv = U10
            elif vdiff <= 40e-6: # 40 uL addition
                uv = U40
            elif vdiff <= 1000e-6: # 1000 uL addition
                uv = U1000

            # cumulative uncertainty with more titrant added
            vol.append(vol[-1] + ufloat(vdiff, uv))
            
        self.vol: NDArray = np.array(list(vol)) # L
        self.ph: NDArray = np.array(list(ph))
        self.gra: NDArray = 10**(-self.ph)*self.vol # mol

        # computing veq

        vol = np.array([v.n if isinstance(v, UFloat) else v for v in self.vol])
        gra = np.array([g.n if isinstance(g, UFloat) else g for g in self.gra])
        # volerr = np.array([v.s if isinstance(v, UFloat) else v for v in self.vol])
        # graerr = np.array([g.s if isinstance(g, UFloat) else g for g in self.gra])

        self.imin: int = min(range(len(vol)), key=lambda i: abs(vol[i] - self.vmin)) # index of the closest point to vmin
        self.imax: int = min(range(len(vol)), key=lambda i: abs(vol[i] - self.vmax)) # index of the closest point to vmax

        fit = linregress(vol[self.imin:self.imax+1], gra[self.imin:self.imax+1])
        self.r2: float = fit.rvalue**2 # coefficient of determination
        self.slope: UFloat = ufloat(fit.slope, fit.stderr) # M
        self.intercept: UFloat = ufloat(fit.intercept, fit.intercept_stderr) # mol
        self.veq: UFloat = -self.intercept/self.slope # L

    def __str__(self):
        string = 'Titration data in (vol/uL, pH): '
        for v, p in zip(self.vol*1e6, self.ph):
            if isinstance(v, UFloat):
                string += f'({v:.1u}, '
            else:
                string += f'({v:.1f}, '
            if isinstance(p, UFloat):
                string += f'{p:.1u}), '
            else:
                string += f'{p:.1f}), '
        return string[:-2].replace('+/-', ' ± ')
        
    @staticmethod
    def _numerr(val: UFloat, prefix: int = 1, sci: bool = True) -> str:
        """Returns a formatted string of `num` ± `err` with the correct dp in scientific notation with optional unit prefix specifications."""

        if not isinstance(val, UFloat):
            return f'{val:.1f}'

        num = val.n*(1/prefix) if np.isfinite(val.n) else 0 # scale by the unit prefix
        err = val.s*(1/prefix) if np.isfinite(val.s) else 0

        i = int(f'{err:.0e}'.split('e')[-1])+1 # find the exponent such that the error has 1sf in the first decimal place

        if sci and prefix==1: # scientific notation to 1dp by default
            return rf'${num*(10**-i):.1f}$ ± ${err*(10**-i):.1f}$ $\times$ 10$^{{{i}}}$'
        else: # use more dp (implied default if a unit prefix is specified)
            return f'${num:.{max(-i+1,0)}f}$ ± ${err:.{max(-i+1,0)}f}$'


    def plot(self) -> None:
        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
        vol = np.array([v.n if isinstance(v, UFloat) else v for v in self.vol])*1e6
        ph = np.array([p.n if isinstance(p, UFloat) else p for p in self.ph])

        ax.scatter(vol, ph, marker='x', color='black') # plot data points, error bars are too small to be seen
        ax.plot(vol, ph, color='black', alpha = 0.2) # plot titration curve
        ax.axvline(x=self.veq.n*1e6, color='black', linestyle='--', alpha=0.2) # plot equivalence point
        ax.text(0.02, 0.97, 
                f'$V_{{\\mathrm{{eq}}}}$ = {self._numerr(self.veq, prefix=1e-6)} μL\n\n'
                f'{self.info}',
                transform=ax.transAxes, ha='left', va='top')

        fig.canvas.draw() # calculate ticks preliminarily
        xscale = np.diff(ax.get_xticks())[0]
        yscale = np.diff(ax.get_yticks())[0]
        xmin = np.floor((min(vol)-0.3*xscale)/xscale)*xscale
        xmax = np.ceil((max(vol)+0.3*xscale)/xscale)*xscale
        ymin = np.floor((min(ph)-0.3*yscale)/yscale)*yscale
        ymax = np.ceil((max(ph)+0.3*yscale)/yscale)*yscale
        ax.set_xlim(xmin, xmax)
        ax.set_xticks(np.arange(xmin, xmax+abs(1e-10*xmax), xscale))
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(np.arange(ymin, ymax+abs(1e-10*ymax), yscale))

        ax.set_xlabel('$V_{\\mathrm{NaOH}}$ / μL')
        ax.set_ylabel('pH')
        if self.desc == '':
            conj = ''
        else:
            conj = ' for '
        ax.set_title(f'Titration curve{conj}{self.desc}')
        plt.show()

    def gran(self) -> None:

        vol = np.array([v.n if isinstance(v, UFloat) else v for v in self.vol])
        gra = np.array([g.n if isinstance(g, UFloat) else g for g in self.gra])
        # volerr = np.array([v.s if isinstance(v, UFloat) else v for v in self.vol])
        # graerr = np.array([g.s if isinstance(g, UFloat) else g for g in self.gra])

        for i in range(2): # plot twice, first time with all points, second time zoomed in around regression line
            fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)

            if i == 0:
                mask = slice(0, len(vol)-1)
            if i == 1:
                mask = slice(self.imin-2, self.imax+int((len(vol)-self.imax)/2))
            x = vol[mask]*1e6
            y = gra[mask]*1e9
            ax.scatter(x, y, marker='x', color='black') # plot data points

            # if i == 0:
            #     ax.scatter(x, y, marker='x', color='black') # plot data points
            # if i == 1: # plot error bars only for linear regression
            #     xerr = volerr[mask]*1e6
            #     yerr = graerr[mask]*1e9
            #     ax.errorbar(x, y, xerr=xerr, yerr=yerr, 
            #                 capsize=2, elinewidth=0.5, capthick=0.5, linestyle='none', color='black')
            
            x1 = np.linspace(vol[self.imin], vol[self.imax], 100)
            x2 = np.linspace(vol[self.imin]*0.9, vol[self.imin], 100)
            x3 = np.linspace(vol[self.imax], vol[self.imax]*1.1, 100)
            y1 = self.slope.n*x1 + self.intercept.n
            y2 = self.slope.n*x2 + self.intercept.n
            y3 = self.slope.n*x3 + self.intercept.n
            ax.plot(x1*1e6, y1*1e9, color='black', alpha=0.2)
            ax.plot(x2*1e6, y2*1e9, color='black', linestyle='--', alpha=0.2)
            ax.plot(x3*1e6, y3*1e9, color='black', linestyle='--', alpha=0.2)
            ax.axhline(y=0, color='black', linestyle='--', alpha = 0.2)

            fig.canvas.draw() # calculate ticks preliminarily
            xscale = np.diff(ax.get_xticks())[0]
            yscale = np.diff(ax.get_yticks())[0]
            xmin = np.floor((min(x)-0.3*xscale)/xscale)*xscale
            xmax = np.ceil((max(x)+0.3*xscale)/xscale)*xscale
            ymin = np.floor((min(y)-0.5*yscale)/yscale)*yscale
            ymax = np.ceil((max(y)+0.5*yscale)/yscale)*yscale
            ax.set_xlim(xmin, xmax)
            ax.set_xticks(np.arange(xmin, xmax+abs(1e-10*xmax), xscale))
            ax.set_ylim(ymin, ymax)
            ax.set_yticks(np.arange(ymin, ymax+abs(1e-10*ymax), yscale))

            ax.text(0.98, 0.97, 
                    f'$V_{{\\mathrm{{eq}}}}$ = x-intercept = {self._numerr(self.veq, prefix=1e-6)} μL\n'
                    f'y-intercept = {self._numerr(self.intercept, prefix=1e-9)} nmol\n'
                    f'slope = {self._numerr(self.slope, prefix=1e-6)} μM\n'
                    f'$r^2$ = {self.r2:.3f}\n\n'
                    f'{self.info}',
                    transform=ax.transAxes, ha='right', va='top')

            ax.set_xlabel('$V_{\\mathrm{NaOH}}$ / μL')
            ax.set_ylabel('[H$\\!^+\\!] \\times V_{\\mathrm{NaOH}}$ / nmol')
            if self.desc == '':
                conj = ''
            else:
                conj = ' for '
            ax.set_title(f'Gran plot{conj}{self.desc}')
            plt.show()