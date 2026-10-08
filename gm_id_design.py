"""
Robust charge amplifier design using gm/ID lookup tables.
SNU Advanced Analog CMOS Design -- Design Problem 1.

Common-source NMOS M1, ideal current-source load, capacitive feedback (CF)
around a sensor capacitance (CS), driving CB at the output.  All device
parameters come from the supplied illustrative gm/ID lookup table.

This script performs:
  Part A: closed-form derivation check (symbolic results are in the report;
          this file evaluates them numerically).
  Part B: feasible-region sweep in (x, ID) for each L, at all three corners.
  Part C: minimum-power (= minimum ID) design search per L, with grid-
          convergence demonstration.
  Part D: parasitic-aware re-check / re-optimization, and the exact
          transient (impulse) response to a 20 fC charge step.
  Part E: fixed-ID vs constant-gm bias-strategy comparison.

Run:  python3 gm_id_design.py
Outputs: PNG plots + a JSON results dump, printed summary tables.
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. Lookup table (illustrative data given in the problem, Section 2)
# ----------------------------------------------------------------------
X_PTS = np.array([5.0, 10.0, 15.0, 20.0, 25.0])            # gm/ID  [S/A]

LUT = {
    0.10e-6: dict(fT=np.array([35, 28, 18, 8, 1.2]) * 1e9,   # Hz
                  A=np.array([15, 20, 27, 35, 42]),          # gm/gds
                  J=np.array([180, 70, 20, 4, 0.30])),       # A/m  (uA/um numerically equal)
    0.20e-6: dict(fT=np.array([22, 17, 10, 4, 0.60]) * 1e9,
                  A=np.array([28, 40, 56, 72, 90]),
                  J=np.array([110, 42, 12, 2.2, 0.16])),
    0.40e-6: dict(fT=np.array([10, 8, 5, 2, 0.30]) * 1e9,
                  A=np.array([55, 80, 110, 145, 180]),
                  J=np.array([55, 21, 6, 1.1, 0.08])),
}
L_LIST = [0.10e-6, 0.20e-6, 0.40e-6]

def interp_A(L, x):
    return np.interp(x, X_PTS, LUT[L]["A"])

def interp_log(L, x, key):
    vals = LUT[L][key]
    return np.exp(np.interp(x, X_PTS, np.log(vals)))

def interp_fT(L, x):
    return interp_log(L, x, "fT")

def interp_J(L, x):
    return interp_log(L, x, "J")

# ----------------------------------------------------------------------
# 2. Fixed specs
# ----------------------------------------------------------------------
CS = 1.00e-12
CF = 0.250e-12
CB = 0.500e-12
VDD = 1.2
VO_LO, VO_HI = 0.35, 0.85
LOAD_COMPLIANCE = 0.15
kB = 1.380649e-23
F_EVAL = 100e6
W_EVAL = 2 * np.pi * F_EVAL

CORNERS = {
    "nominal":   dict(T=300.0, a=1.00, b=1.00, d=1.00, gamma=1.10),
    "slow_hot":  dict(T=398.0, a=0.75, b=1.10, d=1.20, gamma=1.30),
    "fast_cold": dict(T=233.0, a=1.20, b=0.95, d=1.35, gamma=1.00),
}

ID_MAX = 1.00e-3
W_MAX = 300e-6
WL_MAX = 80e-12
VDSAT_MAX = 0.35
H0_MIN = 0.85
CAL_ERR_MAX = 0.05
F3DB_MIN = 100e6
NOISE_MAX = 3.00e-12

# ----------------------------------------------------------------------
# 3. Nominal device sizing (Section 3)
# ----------------------------------------------------------------------
def device_nominal(L, x, ID):
    g = x * ID
    fT = interp_fT(L, x)
    A = interp_A(L, x)
    J = interp_J(L, x)
    W = ID / J
    Cgg = g / (2 * np.pi * fT)
    Cgs, Cgd, Cdb = 0.80 * Cgg, 0.10 * Cgg, 0.30 * Cgg
    gds = g / A
    go = 1.60 * gds
    return dict(g=g, W=W, fT=fT, A=A, J=J, Cgg=Cgg, Cgs=Cgs, Cgd=Cgd, Cdb=Cdb,
                gds=gds, go=go)

# ----------------------------------------------------------------------
# 4. Corner evaluation (Section 4-6): ZT(s), H0, pole/zero, f3dB, noise
# ----------------------------------------------------------------------
def corner_eval(L, x, ID, corner_name, Cpar_in=0.0, Cpar_out=0.0, nom=None):
    if nom is None:
        nom = device_nominal(L, x, ID)
    c = CORNERS[corner_name]
    g_c = c["a"] * nom["g"]
    go_c = c["d"] * nom["go"]
    Cgs_c = c["b"] * nom["Cgs"]
    Cgd_c = c["b"] * nom["Cgd"]
    Cdb_c = c["b"] * nom["Cdb"]

    Cx = CF + Cgd_c
    Ci = CS + Cpar_in + Cgs_c + Cx
    Co = CB + Cpar_out + Cdb_c + Cx

    A1 = Cx * g_c + Ci * go_c
    A2 = Ci * Co - Cx ** 2

    H0 = CF * g_c / A1
    wp = A1 / A2           # rad/s, LHP pole
    wz = g_c / Cx           # rad/s, RHP zero

    term = 1.0 / wp ** 2 - 2.0 / wz ** 2
    if term > 0:
        w3 = 1.0 / np.sqrt(term)
        f3dB = w3 / (2 * np.pi)
        f3dB_exists = True
    else:
        f3dB = None
        f3dB_exists = False

    Hinf = H0 * wp / wz   # high-frequency plateau magnitude (always finite, see derivation)

    Sio = 4 * kB * c["T"] * c["gamma"] * g_c
    Sin = Sio * (W_EVAL ** 2 * Ci ** 2) / (g_c ** 2 + (W_EVAL * Cx) ** 2)
    in_noise = np.sqrt(Sin)

    xc = c["a"] * x
    VDsat_c = 2.0 / xc

    return dict(g=g_c, go=go_c, Cgs=Cgs_c, Cgd=Cgd_c, Cdb=Cdb_c, Cx=Cx, Ci=Ci, Co=Co,
                A1=A1, A2=A2, H0=H0, wp=wp, wz=wz, f3dB=f3dB, f3dB_exists=f3dB_exists,
                Hinf=Hinf, in_noise=in_noise, VDsat=VDsat_c, xc=xc, T=c["T"])

# ----------------------------------------------------------------------
# 5. Full design evaluation / feasibility check
# ----------------------------------------------------------------------
def evaluate_design(L, x, ID, Cpar_in=0.0, Cpar_out=0.0):
    nom = device_nominal(L, x, ID)
    W = nom["W"]
    results = {}
    for cname in CORNERS:
        results[cname] = corner_eval(L, x, ID, cname, Cpar_in, Cpar_out, nom=nom)
    H0_nom = results["nominal"]["H0"]
    Kcal = 1.0 / H0_nom if H0_nom > 0 else np.inf

    ok = True
    reasons = []

    if ID > ID_MAX:
        ok, _ = False, reasons.append(f"ID={ID*1e3:.4f}mA > {ID_MAX*1e3:.2f}mA")
    if W > W_MAX:
        ok, _ = False, reasons.append(f"W={W*1e6:.1f}um > {W_MAX*1e6:.0f}um")
    if W * L > WL_MAX:
        ok, _ = False, reasons.append(f"W*L={W*L*1e12:.2f}um^2 > {WL_MAX*1e12:.0f}um^2")

    for cname, r in results.items():
        if (not r["f3dB_exists"]) or (r["f3dB"] < F3DB_MIN):
            ok = False
            reasons.append(f"{cname}: f3dB fail ({r['f3dB']/1e6 if r['f3dB_exists'] else float('nan'):.1f} MHz)")
        if r["in_noise"] > NOISE_MAX:
            ok = False
            reasons.append(f"{cname}: noise fail ({r['in_noise']*1e12:.3f} pA/rtHz)")
        if r["H0"] < H0_MIN:
            ok = False
            reasons.append(f"{cname}: H0 fail ({r['H0']:.3f})")
        cal_err = abs(r["H0"] * Kcal - 1.0)
        if cal_err > CAL_ERR_MAX:
            ok = False
            reasons.append(f"{cname}: cal_err fail ({cal_err*100:.2f}%)")
        if r["VDsat"] > VDSAT_MAX:
            ok = False
            reasons.append(f"{cname}: VDsat fail ({r['VDsat']*1e3:.1f} mV)")

    load_ok = (VDD - VO_HI) >= LOAD_COMPLIANCE
    return dict(ok=ok, reasons=reasons, results=results, Kcal=Kcal, W=W, nom=nom,
                load_ok=load_ok, L=L, x=x, ID=ID)

# ----------------------------------------------------------------------
# 6. Part B: feasible-region sweep
# ----------------------------------------------------------------------
def sweep_feasible(L, nx=161, nid=161, id_lo=1e-7, id_hi=ID_MAX):
    xs = np.linspace(5.0, 25.0, nx)
    ids = np.geomspace(id_lo, id_hi, nid)
    feas = np.zeros((nid, nx), dtype=bool)
    for i, ID in enumerate(ids):
        for j, x in enumerate(xs):
            d = evaluate_design(L, x, ID)
            feas[i, j] = d["ok"]
    return xs, ids, feas

# ----------------------------------------------------------------------
# 7. Part C: minimum-power design search (grid + local refine)
# ----------------------------------------------------------------------
def min_power_design(L, nx_coarse=81, nid_coarse=400, refine=True):
    xs, ids, feas = sweep_feasible(L, nx=nx_coarse, nid=nid_coarse)
    best = None
    for i, ID in enumerate(ids):
        row_feas = feas[i, :]
        if row_feas.any():
            j = np.argmax(row_feas)
            best = (xs[j], ID)
            break
    if best is None:
        return None
    x0, ID0 = best
    if refine:
        # local 2D refine around (x0, ID0): shrink window, finer grid
        x_win = max(0.5, (xs[1] - xs[0]) * 3)
        id_lo = max(1e-7, ID0 / 3)
        id_hi = min(ID_MAX, ID0 * 1.2)
        xs_f = np.linspace(max(5.0, x0 - x_win), min(25.0, x0 + x_win), 61)
        ids_f = np.geomspace(id_lo, id_hi, 400)
        best_f = None
        for ID in ids_f:
            feasible_x = [x for x in xs_f if evaluate_design(L, x, ID)["ok"]]
            if feasible_x:
                best_f = (feasible_x[0], ID)
                # want smallest feasible x too but ID already minimal for this row
                break
        if best_f is not None:
            x0, ID0 = best_f
    return x0, ID0

def converge_check(L):
    """Demonstrate numerical convergence by halving grid steps."""
    r1 = min_power_design(L, nx_coarse=81, nid_coarse=400)
    r2 = min_power_design(L, nx_coarse=161, nid_coarse=800)
    if r1 is None or r2 is None:
        return r1, r2, None
    dP = abs(r2[1] - r1[1]) / r1[1]
    return r1, r2, dP

# ----------------------------------------------------------------------
# 8. Part D: parasitic-aware re-check
# ----------------------------------------------------------------------
CPAR_IN = 50e-15
CPAR_OUT = 75e-15

def check_with_parasitics(L, x, ID):
    return evaluate_design(L, x, ID, Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT)

def min_power_design_parasitic(L, nx_coarse=81, nid_coarse=400, refine=True):
    xs = np.linspace(5.0, 25.0, nx_coarse)
    ids = np.geomspace(1e-7, ID_MAX, nid_coarse)
    best = None
    for ID in ids:
        row_ok = [evaluate_design(L, x, ID, Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT)["ok"] for x in xs]
        if any(row_ok):
            j = row_ok.index(True)
            best = (xs[j], ID)
            break
    if best is None:
        return None
    x0, ID0 = best
    if refine:
        x_win = max(0.5, (xs[1] - xs[0]) * 3)
        id_lo = max(1e-7, ID0 / 3)
        id_hi = min(ID_MAX, ID0 * 1.2)
        xs_f = np.linspace(max(5.0, x0 - x_win), min(25.0, x0 + x_win), 61)
        ids_f = np.geomspace(id_lo, id_hi, 400)
        for ID in ids_f:
            feasible_x = [x for x in xs_f if evaluate_design(L, x, ID, Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT)["ok"]]
            if feasible_x:
                x0, ID0 = feasible_x[0], ID
                break
    return x0, ID0

# ----------------------------------------------------------------------
# 8b. Part D: quantify minimal spec relaxation when infeasible
# ----------------------------------------------------------------------
def worst_case_metrics(L, x, ID, Cpar_in, Cpar_out):
    """Return worst-corner f3dB (min, Hz) and worst-corner noise (max, A/rtHz),
    plus a bool for whether every OTHER constraint (H0, cal_err, VDsat, ID/W/WL)
    is satisfied."""
    nom = device_nominal(L, x, ID)
    W = nom["W"]
    other_ok = True
    if ID > ID_MAX or W > W_MAX or W * L > WL_MAX:
        other_ok = False
    worst_f3dB = np.inf
    worst_noise = -np.inf
    results = {}
    for cname in CORNERS:
        r = corner_eval(L, x, ID, cname, Cpar_in, Cpar_out, nom=nom)
        results[cname] = r
        f3 = r["f3dB"] if r["f3dB_exists"] else 0.0
        worst_f3dB = min(worst_f3dB, f3)
        worst_noise = max(worst_noise, r["in_noise"])
        if r["H0"] < H0_MIN or r["VDsat"] > VDSAT_MAX:
            other_ok = False
    H0_nom = results["nominal"]["H0"]
    Kcal = 1.0 / H0_nom if H0_nom > 0 else np.inf
    for cname, r in results.items():
        if abs(r["H0"] * Kcal - 1.0) > CAL_ERR_MAX:
            other_ok = False
    return worst_f3dB, worst_noise, other_ok


def fractional_violation(L, x, ID, Cpar_in, Cpar_out):
    """Worst-case fractional shortfall across every spec (f3dB, noise, H0,
    VDsat, calibration error) over all three corners, plus hard ID/W/WL
    limits. 0.0 means every spec is met (feasible); positive means the
    worst-violated spec misses by that fraction."""
    nom = device_nominal(L, x, ID)
    W = nom["W"]
    hard_ok = (ID <= ID_MAX) and (W <= W_MAX) and (W * L <= WL_MAX)
    worst = 0.0
    results = {}
    for cname in CORNERS:
        r = corner_eval(L, x, ID, cname, Cpar_in, Cpar_out, nom=nom)
        results[cname] = r
        f3 = r["f3dB"] if r["f3dB_exists"] else 0.0
        worst = max(worst,
                     max(0.0, (F3DB_MIN - f3) / F3DB_MIN),
                     max(0.0, (r["in_noise"] - NOISE_MAX) / NOISE_MAX),
                     max(0.0, (H0_MIN - r["H0"]) / H0_MIN),
                     max(0.0, (r["VDsat"] - VDSAT_MAX) / VDSAT_MAX))
    H0_nom = results["nominal"]["H0"]
    Kcal = 1.0 / H0_nom if H0_nom > 0 else np.inf
    for r in results.values():
        cal_err = abs(r["H0"] * Kcal - 1.0)
        worst = max(worst, max(0.0, (cal_err - CAL_ERR_MAX) / CAL_ERR_MAX))
    if not hard_ok:
        worst = max(worst, 10.0)
    return worst, results, Kcal


def find_closest_to_feasible(L_list=L_LIST, Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT,
                              nx=161, nid=240, refine=True):
    """Two-stage (coarse grid + local refine) search for the (L, x, ID) point
    that minimizes the worst fractional spec violation -- i.e. the design
    requiring the SMALLEST possible relaxation of whichever spec(s) are
    binding, holding every other spec at its original value."""
    best = None
    for L in L_list:
        for x in np.linspace(5.0, 25.0, nx):
            for ID in np.geomspace(1e-6, ID_MAX, nid):
                w, _, _ = fractional_violation(L, x, ID, Cpar_in, Cpar_out)
                if best is None or w < best[0]:
                    best = (w, L, x, ID)
    if refine and best is not None:
        _, L0, x0, ID0 = best
        xs_f = np.linspace(max(5.0, x0 - 1.0), min(25.0, x0 + 1.0), 161)
        ids_f = np.linspace(max(1e-6, ID0 * 0.85), min(ID_MAX, ID0 * 1.15), 251)
        for x in xs_f:
            for ID in ids_f:
                w, _, _ = fractional_violation(L0, x, ID, Cpar_in, Cpar_out)
                if w < best[0]:
                    best = (w, L0, x, ID)
    return best


# ----------------------------------------------------------------------
# 9. Part D: exact transient response to a charge impulse
# ----------------------------------------------------------------------
def transient_response(L, x, ID, Q=20e-15, corner="nominal", Cpar_in=0.0, Cpar_out=0.0):
    r = corner_eval(L, x, ID, corner, Cpar_in, Cpar_out)
    g, go = r["g"], r["go"]
    Cx, Ci, Co = r["Cx"], r["Ci"], r["Co"]
    A1, A2 = r["A1"], r["A2"]
    wp = A1 / A2
    tau = 1.0 / wp

    vo_inf = -Q * g / A1
    vo_0p = Q * Cx / A2

    def vo_t(t):
        return vo_inf + (vo_0p - vo_inf) * np.exp(-t / tau)

    if abs(vo_inf) > 1e-30:
        ratio = abs((vo_0p - vo_inf) / (0.01 * vo_inf))
        t1pct = tau * np.log(ratio) if ratio > 1 else 0.0
    else:
        t1pct = None

    return dict(vo_inf=vo_inf, vo_0p=vo_0p, tau=tau, t1pct=t1pct, H0=r["H0"], vo_t=vo_t)

# ----------------------------------------------------------------------
# 10. Part E: bias-strategy comparison (fixed-ID vs constant-gm)
# ----------------------------------------------------------------------
def bias_strategy_compare(L, x, ID_nom, Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT):
    nom = device_nominal(L, x, ID_nom)
    out = {}
    for cname, c in CORNERS.items():
        # --- fixed-ID strategy ---
        r_fixed = corner_eval(L, x, ID_nom, cname, Cpar_in, Cpar_out, nom=nom)
        P_fixed = VDD * (ID_nom + 50e-6)

        # --- constant-gm strategy ---
        ID_c = ID_nom / c["a"]
        g_c = c["a"] * nom["g"]        # = nom["g"]  (restored)
        go_c = c["d"] * nom["go"]
        Cgs_c = c["b"] * nom["Cgs"]
        Cgd_c = c["b"] * nom["Cgd"]
        Cdb_c = c["b"] * nom["Cdb"]
        Cx = CF + Cgd_c
        Ci = CS + Cpar_in + Cgs_c + Cx
        Co = CB + Cpar_out + Cdb_c + Cx
        A1 = Cx * g_c + Ci * go_c
        A2 = Ci * Co - Cx ** 2
        H0 = CF * g_c / A1
        wp = A1 / A2
        wz = g_c / Cx
        term = 1.0 / wp ** 2 - 2.0 / wz ** 2
        f3dB = (1.0 / np.sqrt(term)) / (2 * np.pi) if term > 0 else None
        Sio = 4 * kB * c["T"] * c["gamma"] * g_c
        Sin = Sio * (W_EVAL ** 2 * Ci ** 2) / (g_c ** 2 + (W_EVAL * Cx) ** 2)
        in_noise = np.sqrt(Sin)
        xc = c["a"] * x
        VDsat_c = 2.0 / xc
        P_gm = VDD * (ID_c + 50e-6)

        out[cname] = dict(
            fixed=dict(ID=ID_nom, P=P_fixed, H0=r_fixed["H0"], f3dB=r_fixed["f3dB"],
                       noise=r_fixed["in_noise"], VDsat=r_fixed["VDsat"]),
            gm=dict(ID=ID_c, P=P_gm, H0=H0, f3dB=f3dB, noise=in_noise, VDsat=VDsat_c),
        )
    return out

# ----------------------------------------------------------------------
# 11. Rejected-example illustrations for Part B write-up
# ----------------------------------------------------------------------
def find_rejected_examples(L):
    # high-efficiency (weak inversion, low current) candidate
    hi_eff = None
    for x in np.linspace(20, 25, 11):
        for ID in np.geomspace(1e-7, 50e-6, 30):
            d = evaluate_design(L, x, ID)
            if not d["ok"]:
                hi_eff = (x, ID, d["reasons"])
                break
        if hi_eff:
            break
    # high-speed (strong inversion, high current) candidate
    hi_speed = None
    for x in np.linspace(5, 7, 5):
        for ID in np.geomspace(500e-6, ID_MAX, 20):
            d = evaluate_design(L, x, ID)
            if not d["ok"]:
                hi_speed = (x, ID, d["reasons"])
                break
        if hi_speed:
            break
    return hi_eff, hi_speed

# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
if __name__ == "__main__":
    summary = {}

    # ---- Part B & C per L ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    best_designs = {}
    for ax, L in zip(axes, L_LIST):
        xs, ids, feas = sweep_feasible(L, nx=121, nid=240)
        ax.pcolormesh(xs, ids * 1e3, feas, shading="auto", cmap="Greens")
        ax.set_yscale("log")
        ax.set_xlabel("x = gm/ID [S/A]")
        ax.set_ylabel("ID [mA]")
        ax.set_title(f"L = {L*1e6:.2f} um  feasible region")
        r1, r2, dP = converge_check(L)
        best_designs[L] = r2 if r2 is not None else r1
        if r2 is not None:
            ax.plot(r2[0], r2[1] * 1e3, "r*", markersize=14)
        summary[f"L={L*1e6:.2f}um"] = dict(
            coarse=r1, fine=r2,
            dP_convergence_pct=None if dP is None else dP * 100,
        )
    plt.tight_layout()
    plt.savefig("feasible_regions.png", dpi=150)
    plt.close()

    # current-vs-x tradeoff plot at min feasible ID per x, nominal L shown as example per all L
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for L in L_LIST:
        xs = np.linspace(5, 25, 81)
        min_ids = []
        for x in xs:
            ids = np.geomspace(1e-7, ID_MAX, 200)
            feasible_ids = [ID for ID in ids if evaluate_design(L, x, ID)["ok"]]
            min_ids.append(feasible_ids[0] * 1e3 if feasible_ids else np.nan)
        ax.plot(xs, min_ids, label=f"L={L*1e6:.2f}um")
    ax.set_xlabel("x = gm/ID [S/A]")
    ax.set_ylabel("minimum feasible ID [mA]")
    ax.set_yscale("log")
    ax.legend()
    ax.set_title("Current vs. inversion-level tradeoff")
    plt.tight_layout()
    plt.savefig("current_vs_x_tradeoff.png", dpi=150)
    plt.close()

    print("=" * 70)
    print("PART C: minimum-power design per L (coarse -> fine grid)")
    print("=" * 70)
    overall_best = None
    for L in L_LIST:
        r1, r2, dP = converge_check(L)
        print(f"L={L*1e6:.2f}um  coarse: x,ID={r1}  fine: x,ID={r2}  "
              f"convergence dP={dP*100 if dP is not None else float('nan'):.3f}%")
        if r2 is not None:
            P = VDD * (r2[1] + 50e-6)
            print(f"   -> P = {P*1e3:.4f} mW")
            if overall_best is None or P < overall_best[2]:
                overall_best = (L, r2, P)

    print("\nOverall minimum-power design:", overall_best)

    L_opt, (x_opt, ID_opt), P_opt = overall_best
    d_opt = evaluate_design(L_opt, x_opt, ID_opt)
    print(f"\nOptimal device: L={L_opt*1e6:.3f}um  x={x_opt:.3f} S/A  ID={ID_opt*1e6:.3f} uA  "
          f"W={d_opt['W']*1e6:.3f} um  P={P_opt*1e3:.4f} mW")
    print("Corner table (optimal, no extra parasitics):")
    for cname, r in d_opt["results"].items():
        print(f"  {cname:10s}: f3dB={r['f3dB']/1e6 if r['f3dB_exists'] else float('nan'):8.2f} MHz  "
              f"noise={r['in_noise']*1e12:6.3f} pA/rtHz  H0={r['H0']:.4f}  "
              f"VDsat={r['VDsat']*1e3:6.2f} mV  "
              f"cal_err={abs(r['H0']*d_opt['Kcal']-1)*100:5.2f}%")

    print("\n" + "=" * 70)
    print("PART B: rejected example designs")
    print("=" * 70)
    for L in L_LIST:
        hi_eff, hi_speed = find_rejected_examples(L)
        print(f"L={L*1e6:.2f}um")
        print(f"  high-efficiency rejected: x,ID,reasons = {hi_eff}")
        print(f"  high-speed rejected:      x,ID,reasons = {hi_speed}")

    print("\n" + "=" * 70)
    print("PART D: parasitic-aware re-check of optimum design")
    print("=" * 70)
    d_par = check_with_parasitics(L_opt, x_opt, ID_opt)
    print("Feasible with parasitics (same device)?", d_par["ok"])
    if not d_par["ok"]:
        print("Reasons:", d_par["reasons"])
        r_new = min_power_design_parasitic(L_opt)
        print("Re-optimized (same L):", r_new)
        if r_new is not None:
            x2, ID2 = r_new
            P2 = VDD * (ID2 + 50e-6)
            print(f"New P = {P2*1e3:.4f} mW  (penalty = {(P2-P_opt)*1e3:.4f} mW, "
                  f"{(P2/P_opt-1)*100:.2f}%)")
            # Try other Ls too for best parasitic-aware design
            best_par = (L_opt, r_new, P2)
            for L in L_LIST:
                if L == L_opt:
                    continue
                rL = min_power_design_parasitic(L)
                if rL is not None:
                    PL = VDD * (rL[1] + 50e-6)
                    print(f"  alt L={L*1e6:.2f}um -> x,ID={rL}  P={PL*1e3:.4f} mW")
                    if PL < best_par[2]:
                        best_par = (L, rL, PL)
            print("Best parasitic-aware design overall:", best_par)
            L_fin, (x_fin, ID_fin), P_fin = best_par
        else:
            print("NO FEASIBLE DESIGN FOUND at this L even after re-optimization.")
            print("Scanning ALL L for any parasitic-feasible point ...")
            found_any = False
            for L in L_LIST:
                rL = min_power_design_parasitic(L)
                if rL is not None:
                    found_any = True
                    PL = VDD * (rL[1] + 50e-6)
                    print(f"  L={L*1e6:.2f}um -> x,ID={rL}  P={PL*1e3:.4f} mW")
            if not found_any:
                print("CONFIRMED: no feasible design exists anywhere in the "
                      "(L,x,ID) space with the added parasitics (exhaustive grid).")
                print("Finding the design that minimizes the worst fractional "
                      "spec violation (smallest possible relaxation, any axis) ...")
                best = find_closest_to_feasible()
                w_best, L_fin, x_fin, ID_fin = best
                wv, res, Kcal = fractional_violation(L_fin, x_fin, ID_fin, CPAR_IN, CPAR_OUT)
                print(f"  Closest-to-feasible point: L={L_fin*1e6:.3f}um x={x_fin:.4f} "
                      f"ID={ID_fin*1e6:.4f}uA")
                print(f"  Worst fractional spec violation = {wv*100:.4f}%  "
                      "(this is the SMALLEST relaxation, on whichever spec(s) bind, "
                      "that restores feasibility)")
                binding = []
                for cname, r in res.items():
                    f3 = r["f3dB"] / 1e6 if r["f3dB_exists"] else float("nan")
                    cal_err = abs(r["H0"] * Kcal - 1.0) * 100
                    print(f"    {cname:10s}: f3dB={f3:8.3f}MHz  noise={r['in_noise']*1e12:7.4f}pA/rtHz  "
                          f"H0={r['H0']:.5f}  VDsat={r['VDsat']*1e3:6.2f}mV  cal_err={cal_err:5.3f}%")
                    if (not r["f3dB_exists"]) or r["f3dB"] < F3DB_MIN:
                        binding.append(f"{cname}: bandwidth short by "
                                        f"{(F3DB_MIN-f3*1e6)/1e6:.4f} MHz")
                    if r["in_noise"] > NOISE_MAX:
                        binding.append(f"{cname}: noise over by "
                                        f"{(r['in_noise']-NOISE_MAX)*1e12:.5f} pA/rtHz")
                    if r["H0"] < H0_MIN:
                        binding.append(f"{cname}: H0 short by {(H0_MIN-r['H0']):.5f}")
                print("  Binding (violated) spec(s):", binding if binding else "none -- fully feasible")
                relaxed_spec = "; ".join(binding) if binding else "fully feasible (within numerical tolerance)"
                P_fin = VDD * (ID_fin + 50e-6)
                print(f"  LABEL: 'Relaxed design' -- requires relaxing: {relaxed_spec}")
                print(f"  Relaxed design: L={L_fin*1e6:.3f}um x={x_fin:.3f} "
                      f"ID={ID_fin*1e6:.3f}uA  P={P_fin*1e3:.4f}mW "
                      f"(vs. original parasitic-free optimum P={P_opt*1e3:.4f}mW, "
                      f"penalty={(P_fin-P_opt)*1e3:.4f}mW / {(P_fin/P_opt-1)*100:.2f}%)")
            else:
                print("A feasible point exists at another L after all -- see above; "
                      "re-run min_power_design_parasitic logic picked the best one.")
                best_par = None
                for L in L_LIST:
                    rL = min_power_design_parasitic(L)
                    if rL is not None:
                        PL = VDD * (rL[1] + 50e-6)
                        if best_par is None or PL < best_par[2]:
                            best_par = (L, rL, PL)
                L_fin, (x_fin, ID_fin), P_fin = best_par
    else:
        L_fin, x_fin, ID_fin, P_fin = L_opt, x_opt, ID_opt, P_opt
        print("Design remains feasible with parasitics added; no re-optimization needed.")

    d_fin = evaluate_design(L_fin, x_fin, ID_fin, Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT)
    print(f"\nFinal parasitic-aware design: L={L_fin*1e6:.3f}um x={x_fin:.3f} ID={ID_fin*1e6:.3f}uA "
          f"W={d_fin['W']*1e6:.3f}um P={P_fin*1e3:.4f}mW  feasible={d_fin['ok']}")
    for cname, r in d_fin["results"].items():
        print(f"  {cname:10s}: f3dB={r['f3dB']/1e6 if r['f3dB_exists'] else float('nan'):8.2f} MHz  "
              f"noise={r['in_noise']*1e12:6.3f} pA/rtHz  H0={r['H0']:.4f}  "
              f"VDsat={r['VDsat']*1e3:6.2f} mV")

    print("\n" + "=" * 70)
    print("PART D (cont.): transient response to 20 fC charge impulse (nominal corner)")
    print("=" * 70)
    tr = transient_response(L_fin, x_fin, ID_fin, Q=20e-15, corner="nominal",
                             Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT)
    print(f"vo(0+) initial jump = {tr['vo_0p']*1e3:.4f} mV")
    print(f"vo(inf) final value = {tr['vo_inf']*1e3:.4f} mV")
    print(f"tau = {tr['tau']*1e9:.4f} ns")
    print(f"t_1% settle = {tr['t1pct']*1e9:.4f} ns" if tr['t1pct'] is not None else "t_1% undefined")
    print(f"H0 (nominal, parasitic-aware) = {tr['H0']:.4f}")
    print(f"Swing check: |vo(inf)| = {abs(tr['vo_inf'])*1e3:.4f} mV vs 250 mV budget "
          f"-> {'OK' if abs(tr['vo_inf'])<=0.25 else 'VIOLATION'}")

    # transient plot
    fig, ax = plt.subplots(figsize=(6, 4.5))
    tt = np.linspace(0, 8 * tr["tau"], 400)
    ax.plot(tt * 1e9, tr["vo_t"](tt) * 1e3)
    ax.axhline(tr["vo_inf"] * 1e3, ls="--", c="gray", lw=0.8)
    ax.set_xlabel("time [ns]")
    ax.set_ylabel("vo(t) [mV]")
    ax.set_title("Charge-impulse acquisition response (20 fC, nominal corner)")
    plt.tight_layout()
    plt.savefig("transient_response.png", dpi=150)
    plt.close()

    print("\n" + "=" * 70)
    print("PART E: fixed-ID vs constant-gm bias strategy")
    print("=" * 70)
    cmp_ = bias_strategy_compare(L_fin, x_fin, ID_fin, Cpar_in=CPAR_IN, Cpar_out=CPAR_OUT)
    for cname, d in cmp_.items():
        print(f"\n{cname}:")
        for strat in ("fixed", "gm"):
            s = d[strat]
            f3 = s["f3dB"] / 1e6 if s["f3dB"] is not None else float("nan")
            print(f"  {strat:6s}: ID={s['ID']*1e6:7.3f}uA  P={s['P']*1e3:7.4f}mW  "
                  f"f3dB={f3:8.2f}MHz  noise={s['noise']*1e12:6.3f}pA/rtHz  "
                  f"H0={s['H0']:.4f}  VDsat={s['VDsat']*1e3:6.2f}mV")

    # Save everything to JSON for the report
    dump = dict(
        optimal_nominal=dict(L=L_opt, x=x_opt, ID=ID_opt, P=P_opt, W=d_opt["W"]),
        optimal_parasitic=dict(L=L_fin, x=x_fin, ID=ID_fin, P=P_fin, W=d_fin["W"],
                                feasible=d_fin["ok"]),
        transient=dict(vo_0p=tr["vo_0p"], vo_inf=tr["vo_inf"], tau=tr["tau"],
                        t1pct=tr["t1pct"], H0=tr["H0"]),
    )
    with open("results_summary.json", "w") as f:
        json.dump(dump, f, indent=2, default=str)

    print("\nDone. Plots: feasible_regions.png, current_vs_x_tradeoff.png, "
          "transient_response.png. Summary: results_summary.json")
