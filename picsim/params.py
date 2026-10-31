"""Central parameter registry. Every number used by the models is defined here once.

Edit values here, not inside scripts. Each entry carries a provenance tag (see provenance.py).
"""
from .provenance import P

PROPOSAL = "ANRF/ARG/2026/009937 draft v97"
ZHAO = "Zhao et al., arXiv:2508.09920 (abstract and Sec. 2, read 2026-10-31)"
MONTIFIORE = "Montifiore et al., arXiv:2601.15695 (abstract, read 2026-10-31)"
TRC = "arXiv:2606.05673, introduction (values at 193 THz; ~7% dispersion to visible), read 2026-10-31"

# ---------------------------------------------------------------- wavelengths
LAM_TRANSPORT = P(1.0644350, "um", "USER", f"{PROPOSAL}, Table 2: side-B tweezers at 1064.4350 nm (45S1/2 double magic)")
LAM_STATIC = P(1.075, "um", "USER", f"{PROPOSAL}, Table 2: static array at 1075 nm")
LAM_ADDRESS = P(0.8946, "um", "USER", f"{PROPOSAL}: 895 nm Raman addressing (Cs D1 line, 894.6 nm)")
LAM_IMAGE = P(0.852, "um", "USER", f"{PROPOSAL}: 852 nm imaging (Cs D2)")

# ---------------------------------------------------------------- platform
T_SIN = P(0.400, "um", "USER", f"{PROPOSAL}, Sec. 7.4: CoE-CPPICS stoichiometric LPCVD SiN up to 400 nm")
T_SIN_SWEEP = P((0.200, 0.300, 0.400), "um", "CHOICE", "thickness sweep agreed with PI; 200 nm covers foundry runs")
T_BOX = P(4.0, "um", "CHOICE", "buried thermal oxide; checked against substrate leakage in scripts/m1_platform.py")
T_TOPCLAD = P(3.0, "um", "CHOICE", "upper oxide cladding under heater")
MIN_FEATURE = P(0.100, "um", "CHOICE", "minimum e-beam feature assumed for taper tips and gaps")

# ---------------------------------------------------------------- materials
N_SIN_MODEL = P(0.0, "-", "SOURCE", "Luke et al., Opt. Lett. 40, 4823 (2015), doi:10.1364/OL.40.004823; "
                "Sellmeier via refractiveindex.info Si3N4/Luke (340 nm LPCVD film), read 2026-10-31")
N_SIO2_MODEL = P(0.0, "-", "SOURCE", "Malitson, JOSA 55, 1205 (1965), doi:10.1364/JOSA.55.001205; "
                 "via refractiveindex.info SiO2/Malitson, read 2026-10-31")
N_ALN_MODEL = P(0.0, "-", "SOURCE", "Pastrnak & Roskovcova, Phys. Stat. Sol. 14, K5 (1966), doi:10.1002/pssb.19660140127; "
                "ordinary index via refractiveindex.info AlN/Pastrnak-o, read 2026-10-31")
N_SI_1064 = P(3.55, "-", "RECALLED", "crystalline Si near 1.06 um, approx. 3.55 + 0.0001i; verify (Green 2008)")
DNDT_SIN = P(2.45e-5, "1/K", "SOURCE", TRC)
DNDT_SIO2 = P(1.0e-5, "1/K", "SOURCE", TRC)
N2_SIN = P(2.5e-19, "m^2/W", "RECALLED", "Kerr index of LPCVD Si3N4, approx. 2.4-2.5e-19 m^2/W (Ikeda 2008); verify")

K_SIO2 = P(1.38, "W/(m K)", "RECALLED", "fused/thermal silica thermal conductivity; verify")
K_SI = P(148.0, "W/(m K)", "RECALLED", "silicon thermal conductivity at 300 K; verify")
K_SIN = P(20.0, "W/(m K)", "RECALLED", "LPCVD SiN film, literature spread ~10-30 W/(m K); verify")
K_AIR = P(0.026, "W/(m K)", "RECALLED", "air at 300 K")
K_METAL = P(20.0, "W/(m K)", "RECALLED", "TiN/NiCr heater film order of magnitude; verify")
RHO_CP_SIO2 = P(2200 * 740, "J/(m^3 K)", "RECALLED", "density 2200 kg/m^3 x c_p 740 J/(kg K)")
RHO_CP_SI = P(2330 * 700, "J/(m^3 K)", "RECALLED", "density 2330 kg/m^3 x c_p 700 J/(kg K)")
RHO_CP_SIN = P(3100 * 700, "J/(m^3 K)", "RECALLED", "density 3100 kg/m^3 x c_p 700 J/(kg K)")

# ---------------------------------------------------------------- fibre interfaces
MFD_LENSED = P(2.5, "um", "CHOICE", "lensed/tapered fibre spot, typical catalogue range 2-3 um")
MFD_HI1060 = P(6.2, "um", "RECALLED", "Corning HI1060 MFD at 1060 nm approx. 6.2 um; verify datasheet")
MFD_PM780 = P(5.3, "um", "RECALLED", "PM780-HP MFD at 850 nm approx. 5.3 um; verify datasheet")

# ---------------------------------------------------------------- system targets (Table 1, Sec. 4)
NA_OBJ = P(0.5, "-", "USER", f"{PROPOSAL}, Table 1: meta-objective NA >= 0.5")
FIELD_HALF = P(100.0, "um", "USER", f"{PROPOSAL}, Table 1: Strehl >= 0.8 over >= +/-100 um")
WD = P(15.0, "mm", "USER", f"{PROPOSAL}, Sec. 4(i): similar working distance to 1.5 cm")
SITE_PITCH = P(4.0, "um", "USER", f"{PROPOSAL}, Sec. 3: 'about five sites at 4 um pitch'")
N_OPA = P((32, 64), "-", "USER", f"{PROPOSAL}, Sec. 4 Phase C: 32-64 channels")
K_TWEEZERS = P(21, "-", "USER", f"{PROPOSAL}, Table 1: >20 simultaneous tweezers; Sec. 4: twenty-one tweezers")
UPDATE_RATE = P(100e3, "Hz", "USER", f"{PROPOSAL}, Table 1: >= 100 kHz")
POS_NOISE = P(10e-3, "um", "USER", f"{PROPOSAL}, Sec. 4 error budget: <= 10 nm rms")
P_TWEEZER_ATOMS = P(6.6e-3, "W", "USER", f"{PROPOSAL}, Sec. 4(ii): ~6.6 mW per tweezer at the atoms")
LOSS_CHIP_OPA = P(10.0, "dB", "USER", f"{PROPOSAL}, Sec. 4(ii): ~10 dB of chip and OPA loss")
P_CHIP = P(1.4, "W", "USER", f"{PROPOSAL}, Sec. 4(ii): ~1.4 W on chip for 21 tweezers")
W0_TWEEZER = P(1.04, "um", "USER", f"{PROPOSAL}, Sec. 4 error budget: 1.04 um waist")
ER_TARGET = P(50.0, "dB", "USER", f"{PROPOSAL}, Table 1: >= 50 dB at 895 nm")
T_SWITCH = P(100e-9, "s", "USER", f"{PROPOSAL}, Table 1: <= 100 ns")
N_ROUTER = P(16, "-", "USER", f"{PROPOSAL}, Sec. 4: 16 channels on the 895 nm router")
P_SITE_ADDR = P((1e-6, 20e-6), "W", "USER", f"{PROPOSAL}, Sec. 4: ~1-20 uW per site at the atoms")
FACET_LOSS = P((1.0, 2.0), "dB", "USER", f"{PROPOSAL}, Sec. 4: ~1-2 dB per facet")
TREE_LOSS = P(1.0, "dB", "USER", f"{PROPOSAL}, Sec. 4: 5-6 MMI stages, ~1 dB")

# ---------------------------------------------------------------- actuators
ALN_VPI_795 = P(74.7, "V", "SOURCE", f"{ZHAO}: single-ended V_pi 74.7 +/- 3.7 V at 795 nm")
ALN_VPI_1013 = P(200.0, "V", "SOURCE", f"{ZHAO}: V_pi 200 V at 1013 nm")
ALN_RISE = P(26e-9, "s", "SOURCE", f"{ZHAO}: 10-90% rise 26 +/- 7 ns")
ALN_ER_795 = P(71.4, "dB", "SOURCE", f"{ZHAO}: mean ER 71.4 +/- 1.1 dB, cascaded MZI per channel")
ALN_XTALK = P(-68.0, "dB", "SOURCE", f"{ZHAO}: nearest-neighbour on-chip crosstalk -68.0 +/- 1.0 dB")
PZT_VPI_532 = P(2.8, "V", "SOURCE", f"{MONTIFIORE}: coil MZM V_pi 2.8 V at 532 nm")
PZT_BW = P(0.4e6, "Hz", "SOURCE", f"{MONTIFIORE}: coil MZM 3-dB bandwidth 0.4 MHz")
PZT_ER = P(21.5, "dB", "SOURCE", f"{MONTIFIORE}: coil MZM ER 21.5 dB")

# ---------------------------------------------------------------- atom
CS_MASS = P(132.905451933 * 1.66053906660e-27, "kg", "RECALLED", "133Cs atomic mass (Steck, Cesium D line data)")
ALPHA_CS_1064 = P(1160.0, "a.u.", "USER", f"{PROPOSAL}, Sec. 3 citing Le Kien et al. EPJD 67, 92 (2013)")
ALPHA_CS_1075 = P(1120.0, "a.u.", "USER", f"{PROPOSAL}, Sec. 4(ii)")
DEPTH_TRANSPORT = P(1.2e-3, "K", "USER", f"{PROPOSAL}, Sec. 4(ii): transport tweezers ~1.2 mK")
SURVIVAL_MOVE = P(0.99, "-", "USER", f"{PROPOSAL}, Table 1: >= 99% per move")
