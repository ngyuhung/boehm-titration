from molmass import Formula

FIGSIZE = (8, 5) # figure size for all plots
DPI = 400 # figure dpi resolution

"""
Dictionary of LaTeX text formatting for Matplotlib rendering
"""

X = {'NaHCO3': 'NaHCO$_3$',
     'Na2CO3': 'Na$_2$CO$_3$',
     'NaOH': 'NaOH',
     'HCl': 'HCl',
     'KHP': 'KHP',}

"""
Dictionary of molar masses
"""

M = {'NaHCO3': Formula('NaHCO3').mass,
     'Na2CO3': Formula('Na2CO3').mass,
     'NaOH': Formula('NaOH').mass,
     'HCl': Formula('HCl').mass,
     'KHP': Formula('C8H5O4K').mass,}

"""
Uncertainties
"""

U10 = 0.3e-6 # uncertainty/uL of 10 uL additions using a 10-100 uL Eppendorf Research plus micropipette
U40 = 0.5e-6 # uncertainty/uL of 40 uL additions using a 10-100 uL Eppendorf Research plus micropipette
U1000 = 6e-6 # uncertainty/uL of 1000 uL additions using a 100-1000 uL Eppendorf Research plus micropipette

UPH = 0.01 # uncertainty of MageTech pHTemp2000 meter
UMASS = 0.1e-3 # uncertainty/g of electronic mass balance
UVOLFLASK = 0.1e-3 # uncertainty/L of volumetric flask