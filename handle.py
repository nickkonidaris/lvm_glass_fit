
import numpy as np
from astropy.modeling import models, fitting
from pylab import *
rcParams.update({"figure.autolayout": True})


# Lines from https://www.us.schott.com/d/advanced_optics/38cfacda-3e20-49d5-a432-d9066f56c05b/1.2/tie_29_refractive_index_and_dispersion_us.pdf

n_bak2_catalog = np.array([[1529.6, 1060.0, 1013.98, 852.11, 706.5188, 656.2725, 643.8469, 632.8, 589.2938, 587.5618, 546.074, 486.1327, 479.9914, 435.8343, 404.6561, 365.0146], [1.52385,1.52919,1.52980,1.53234,1.53564,1.53721,1.53765,1.53806,1.53988,1.53996,1.54212,1.54625,1.54677,1.55117,1.55525,1.56221]])

n_bak2_catalog[0] /= 1000 # To micron

n_bak2_coeff = np.array([1.01662154E+00,3.19903051E-01,9.37232995E-01,5.92383763E-03,2.03828415E-02,1.13118417E+02])


# from https://www.us.schott.com/advanced_optics/english/download/index.html https://www.us.schott.com/d/advanced_optics/0dccafb2-f1b8-499e-98b1-b6c3d09ee9bd/1.45/schott-optical-glass-overview-excel-table-english-january-2019.xls
def sellmeiern(lam_micron, B1=n_bak2_coeff[0],B2=n_bak2_coeff[1],B3=n_bak2_coeff[2],C1=n_bak2_coeff[3],C2=n_bak2_coeff[4],C3=n_bak2_coeff[5]):
    " Returns n for Sellmeier "
    l = lam_micron**2

    if np.any((300 < lam_micron) & (lam_micron < 1000)):
        print("LIKELY YOUR WAVELENGTH IS INCORRECT")

    return np.sqrt(B1*l/(l-C1) + B2*l/(l-C2) + B3*l/(l-C3)+1)

#def simple_sellmeiern(lam_micron, B1=n_bak2_coeff[0],B2=n_bak2_coeff[1],C1=n_bak2_coeff[3],C2=n_bak2_coeff[4]):
def simple_sellmeiern(lam_micron, B1=n_bak2_coeff[0],B2=n_bak2_coeff[1],C1=n_bak2_coeff[3],C2=n_bak2_coeff[4]):
    " Returns n for Sellmeier "
    l = lam_micron**2

    if np.any((300 < lam_micron) & (lam_micron < 1000)):
        print("LIKELY YOUR WAVELENGTH IS INCORRECT")

    return np.sqrt(B1*l/(l-C1) + B2*l/(l-C2) + 1)
    #return np.sqrt(B1*l/(l-C1) + B2*l/(l-C2) + n_bak2_coeff[2]*l/(l-n_bak2_coeff[5]) + 1)



model_sellmeiern = models.custom_model(simple_sellmeiern)
s_init = model_sellmeiern()
fit_t = fitting.LevMarLSQFitter()

def check_catalog():
    " Compare N-BAK2 catalog values against the fitted values. "
    ls,ns = n_bak2_catalog
    ff = fit_t(s_init, *n_bak2_catalog, maxiter=10000)
    print("Fit Sellmeier coefficients to catalog N-BAK2 data")
    print("%10s %10s %10s %15s" % ("Wave", "Cat.", "Fit", "Delta"))
    for ix, l in enumerate(ls):
        n = ns[ix]
        fn = ff(l)

        print("%10.4f %10.6f %10.6f %15.2e" % (l,n,fn, n-fn))
        #assert(abs(n-fn) < 1e-5)

check_catalog()


def dict_to_list(gd):
    " Take a glass dictionary (wave: index) and convert to wavelength, index array "

    l = np.array(list(gd.keys()))
    n = np.array(list(gd.values()))

    return l,n
    

def fit_glass(gd, meltname="", compare=None, Plot=True, model="two"):
    " Fit a glass dictionary and return a fitting function "

    if model == "two":
        model_sellmeiern = models.custom_model(simple_sellmeiern)
    if model == "three":
        model_sellmeiern = models.custom_model(sellmeiern)
    s_init = model_sellmeiern()
    fit_t = fitting.LevMarLSQFitter()

    ls, ns = dict_to_list(gd)
    rs = np.zeros_like(ns)
    offsets = np.zeros_like(ns)
    ff = fit_t(s_init, ls, ns, maxiter=10000)

    print("Fit Sellmeier coefficients to melt: %s" % meltname)
    print("%10s %10s %10s %15s %15s" % ("Wave", "Melt", "Fit", "Melt-Fit", "Melt-Cat."))
    for ix, l in enumerate(ls):
        n = ns[ix]
        fn = ff(l)
        rs[ix] = n-fn

        if compare is None: d = None
        else: d = n-compare[l]
        offsets[ix] = d
        print("%10.4f %10.6f %10.6f %15.2e %15.2e" % (l, n, fn, n-fn, d))

    if Plot:
        figure(1, figsize=(9,9))
        clf()
        subplot(311)
        ll = np.arange(.35,1,.01)
        plot(ll, ff(ll) - sellmeiern(ll))
        grid(True)
        title("Fit Function - N-BAK2 function: %s" % (ff.parameters))
        xlabel("Wavelength [micron]")
        
        subplot(312)
        plot(ls, rs, 'o')
        rms= np.std(rs)
        title("Residual of Melt-Fit: %s (SD: %5.1e)" % (meltname, rms))
        xlabel("Wavelength [micron]")
        ylabel("Delta Melt - Fit")
        ylim(-5e-6,5e-6)
        grid(True)
        tight_layout()

        subplot(313)
        plot(ls, offsets, 'o')
        title("Difference of Melt-Catalog: %s" % meltname)
        xlabel("Wavelength [micron]")
        ylabel("Delta Melt - Catalog")
        grid(True)
        tight_layout()

        
        savefig("%s_%s.pdf" % (meltname, model))

    return ff

d = 0.5875618 ; e = 0.546074 ; r = 0.7065188 ; Fp = 0.4799914 ; Cp = 0.6438469 
g = 0.4358343 ; h = 0.4046561 ; i = 0.3650146 ; s = 0.85211 ; t = 1.01398

N_BAK2 = {d: 1.53996, e: 1.54212, r: 1.53564, Fp: 1.54677, Cp: 1.53765,
        g: 1.55117, h: 1.55525, i: 1.56221, s: 1.53234, t: 1.52980}

B1236231A = {d: 1.5399, e: 1.542057, r: 1.53558, Fp: 1.546715, Cp: 1.537585,
        g: 1.551121, h: 1.555210, s: 1.532267, t: 1.529729, i: 1.56216}
B1236215A = {d: 1.539893, e: 1.542050, r: 1.535573, Fp: 1.546708, Cp: 1.537578,
        g: 1.551114, h: 1.555203, s: 1.532260, t: 1.529722, i: 1.56215}
B1236123A = {d: 1.539818, e: 1.541977, r: 1.535499, Fp: 1.546628, Cp: 1.537505,
        g: 1.551033, h: 1.555114, s: 1.532191, t: 1.529649, i: 1.56205}

ff=fit_glass(N_BAK2, meltname="Catalog", compare=N_BAK2)
#ff=fit_glass(B1236231A, meltname="B1236231A", compare=N_BAK2)
#ff=fit_glass(B1236215A, meltname="B1236215A", compare=N_BAK2)
#ff=fit_glass(B1236123A, meltname="B1236123A", compare=N_BAK2)


pars = []
meltnames = ["B1236231A", "B1236215A", "B1236123A"]
strs = ["%10s %15s %15s %15s %15s %15s %15s" % ("melt", "B1", "B2", "B3", "C1", "C2", "C3")]
for melt in meltnames:
    ff = fit_glass(eval(melt), meltname=melt, compare=N_BAK2, model="three")
    strs.append("%10s %s" % (melt, np.array_str(ff.parameters, max_line_width=50000)))

print("\n".join(strs))
