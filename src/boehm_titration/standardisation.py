from .settings import *
from .titration import _Titration
from uncertainties import ufloat, UFloat


class Standardisation(_Titration):

    def __init__(self, unknown: str, standard: str, data: dict, vmin: float, vmax: float, vanalyte: float, vtitrant: float, 
                 cstandard: UFloat=None, mstandard: float=None, vstandard: float=None):

        super().__init__(data, vmin, vmax)

        self.desc = f'{X[unknown]} standardisation with {X[standard]}'

        self.vanalyte = ufloat(vanalyte, vanalyte/1e-3*U1000)
        self.vtitrant = ufloat(vtitrant, vtitrant/1e-3*U1000)

        if cstandard:
            self.cstandard = cstandard
        elif mstandard:
            self.mstandard = ufloat(mstandard, UMASS)
            self.vstandard = ufloat(vstandard, UVOLFLASK)
            self.cstandard = (self.mstandard/M[standard])/self.vstandard

        # NaOH is always the titrant
        if standard == "NaOH":
            self.c = self.cstandard*(self.veq+self.vtitrant)/self.vanalyte
            self.analyte = unknown
        elif unknown == "NaOH":
            self.c = self.cstandard*self.vanalyte/(self.veq+self.vtitrant)
            self.analyte = standard

        self.info = (f'$V_{{\\mathrm{{{self.analyte}}}}}$ = {self._numerr(self.vanalyte, prefix=1e-3)} mL\n'
                     f'$V_{{\\mathrm{{NaOH}}}}$ = {self._numerr(self.vtitrant, prefix=1e-3)} mL\n\n'
                     f'[{X[standard]}] = {self._numerr(self.cstandard, sci=False)} M\n'
                     f'[{X[unknown]}] = {self._numerr(self.c, sci=False)} M')

    def __str__(self):
        string = self.info + '\n\n' + super().__str__()
        for substr in ['$', '{', '}', '\\mathrm', '_']:
            string = string.replace(substr, '')
        return string

