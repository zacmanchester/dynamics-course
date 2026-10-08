"""16.07 Problem Set 1 -- spring-loaded inverted pendulum.
Generates all figures and numbers used in the solutions document.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3,
                     "pdf.fonttype": 42})

# ---------------------------------------------------------------- parameters
m, l, g = 1.0, 1.0, 9.8
k_cr = m * g * l
ratios = [0.5, 1.0, 1.5]
ks = [r * k_cr for r in ratios]
labels = [r"$k=0.5\,k_{cr}$", r"$k=k_{cr}$", r"$k=1.5\,k_{cr}$"]

GREEN, RED, BLUE = "#1a9850", "#d73027", "#2166ac"


def U(th, k):
    return m * g * l * np.cos(th) + 0.5 * k * th**2


def Upp(th, k):
    return k - m * g * l * np.cos(th)


def energy(th, thd, k):
    return 0.5 * m * l**2 * thd**2 + U(th, k)


def rhs(t, y, k):
    th, thd = y
    return [thd, (g / l) * np.sin(th) - k / (m * l**2) * th]


def equilibria(k):
    """Return list of (theta_eq, kind) in [-3,3]."""
    eq = [(0.0, None)]
    r = k / k_cr
    if r < 1:
        ts = brentq(lambda x: np.sin(x) / x - r, 1e-3, np.pi)
        eq += [(ts, None), (-ts, None)]
    out = []
    for th, _ in eq:
        u2 = Upp(th, k)
        if abs(u2) < 1e-9:
            kind = "degenerate"
        elif u2 > 0:
            kind = "stable"
        else:
            kind = "unstable"
        out.append((float(th), kind, float(u2)))
    return out


results = {"k_cr": k_cr}

# ------------------------------------------------------- Fig 1: U(theta), 3(b-d)
th = np.linspace(-3, 3, 1201)
fig, axs = plt.subplots(3, 1, figsize=(6.4, 8.2), sharex=True)
for ax, k, lab, r in zip(axs, ks, labels, ratios):
    ax.plot(th, U(th, k), "k", lw=1.8)
    for te, kind, u2 in equilibria(k):
        if kind == "stable":
            ax.plot(te, U(te, k), "o", ms=9, color=GREEN, mec="k", zorder=5)
        elif kind == "unstable":
            ax.plot(te, U(te, k), "X", ms=10, color=RED, mec="k", zorder=5)
        else:
            ax.plot(te, U(te, k), "D", ms=9, color=BLUE, mec="k", zorder=5)
    ax.set_ylabel(r"$U(\theta)$")
    ax.set_title(lab)
axs[-1].set_xlabel(r"$\theta$ [rad]")
from matplotlib.lines import Line2D
handles = [
    Line2D([], [], marker="o", ls="", ms=9, color=GREEN, mec="k", label="stable (local min, $U''>0$)"),
    Line2D([], [], marker="X", ls="", ms=10, color=RED, mec="k", label="unstable (local max, $U''<0$)"),
    Line2D([], [], marker="D", ls="", ms=9, color=BLUE, mec="k", label="stable ($U''=0$, but $U$ convex)"),
]
fig.legend(handles=handles, loc="lower center", ncol=1, frameon=True, bbox_to_anchor=(0.5, -0.0))
fig.tight_layout(rect=(0, 0.09, 1, 1))
fig.savefig("fig_potential.pdf")
plt.close(fig)

results["equilibria"] = {f"{r}": equilibria(k) for r, k in zip(ratios, ks)}

# separatrix quantities for k = 0.5 k_cr
k05 = 0.5 * k_cr
th_star = [e for e in equilibria(k05) if e[0] > 0][0][0]
E_sep = U(0.0, k05)                           # = m g l
th_max = brentq(lambda x: 2 * m * g * l * (1 - np.cos(x)) - k05 * x**2, 2.0, 3.0)
results.update({"theta_star": th_star, "E_sep": E_sep, "theta_max": th_max,
                "Upp_star": float(Upp(th_star, k05))})

# ------------------------------------------------------------ simulations
N = 20
# Draw theta0 ~ U(-3,3); use the first seed whose batch contains at least one
# above-barrier run (|theta0| > theta_max ~ 2.78, only ~7% of the range) so that
# all three qualitative orbit types appear in the k = 0.5 k_cr plot.
for SEED in range(200):
    th0s = np.random.default_rng(SEED).uniform(-3, 3, N)
    if (np.sum(np.abs(th0s) > th_max) >= 1 and np.sum(th0s > 0) >= 6 and np.sum(th0s < 0) >= 6
            and np.sum(np.abs(th0s) < 0.3) == 0):
        break
results["seed"] = SEED
t_end = 60.0
t_eval = np.linspace(0, t_end, 12001)
sols = {}
stats = {}
for r, k in zip(ratios, ks):
    trajs = []
    for th0 in th0s:
        s = solve_ivp(rhs, (0, t_end), [th0, 0.0], args=(k,), method="RK45",
                      t_eval=t_eval, rtol=1e-10, atol=1e-12)
        trajs.append(s.y)
    sols[r] = np.array(trajs)                 # (N, 2, T)
    Y = sols[r]
    E = energy(Y[:, 0], Y[:, 1], k)
    drift = np.max(np.abs(E - E[:, :1]) / np.maximum(np.abs(E[:, :1]), 1.0))
    stats[r] = {"max_rel_energy_drift": float(drift),
                "max_abs_theta": float(np.max(np.abs(Y[:, 0])))}
results["sim"] = stats
results["theta0"] = [float(x) for x in th0s]

# ------------------------------------------------------ Fig 2: phase planes, 4(b)
cmap = plt.get_cmap("viridis")
def eq_markers(ax, k):
    for te, kind, _ in equilibria(k):
        if kind == "stable":
            ax.plot(te, 0, "o", ms=9, color=GREEN, mec="k", zorder=6)
        elif kind == "unstable":
            ax.plot(te, 0, "X", ms=10, color=RED, mec="k", zorder=6)
        else:
            ax.plot(te, 0, "D", ms=9, color=BLUE, mec="k", zorder=6)

for r, k, lab, name in zip(ratios, ks, labels, ["05", "10", "15"]):
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    Y = sols[r]
    for i in range(N):
        ax.plot(Y[i, 0], Y[i, 1], lw=0.9, color=cmap(i / (N - 1)))
        ax.plot(Y[i, 0, 0], Y[i, 1, 0], ".", color="k", ms=4)
    eq_markers(ax, k)
    ax.set_xlabel(r"$\theta$ [rad]")
    ax.set_ylabel(r"$\dot\theta$ [rad/s]")
    ax.set_title(f"Phase plane, {lab}  ({N} runs, $t\\in[0,{t_end:.0f}]$ s)")
    fig.tight_layout()
    fig.savefig(f"fig_phase_{name}.pdf")
    plt.close(fig)

# -------------------------- Fig 3: k = k_cr zoom near origin with small ICs
fig, ax = plt.subplots(figsize=(6.2, 4.4))
small = [0.1, 0.2, 0.4, 0.8]
for j, a in enumerate(small):
    for sgn in (+1, -1):
        # integrate for ~1.5 quartic-well periods
        c = g / (6 * l)
        T_est = 7.416 / (a * np.sqrt(c))
        s = solve_ivp(rhs, (0, 1.5 * T_est), [sgn * a, 0.0], args=(k_cr,),
                      method="RK45", rtol=1e-11, atol=1e-13,
                      t_eval=np.linspace(0, 1.5 * T_est, 4000))
        ax.plot(s.y[0], s.y[1], lw=1.2, color=cmap(j / (len(small) - 1)),
                label=rf"$\theta(0)=\pm{a}$" if sgn == 1 else None)
ax.plot(0, 0, "D", ms=9, color=BLUE, mec="k", zorder=6)
ax.set_xlabel(r"$\theta$ [rad]"); ax.set_ylabel(r"$\dot\theta$ [rad/s]")
ax.set_title(r"Zoom near upright equilibrium, $k=k_{cr}$")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout(); fig.savefig("fig_phase_cr_zoom.pdf"); plt.close(fig)

# ---------------------------------------- separatrix as a simulated trajectory
# (Sec. 4(c) hint): the separatrix passes through the saddle at theta=0 and has
# a turning point (theta_dot=0) at theta = +-theta_max, since U(theta_max)=U(0).
# Starting at rest exactly at theta_max puts the IC exactly on the separatrix;
# integrating forward/backward traces the loop (it only reaches the saddle
# asymptotically as t -> +-infinity, so a finite-time run traces the loop
# almost all the way around, slowing to a crawl near the saddle).
t_sep = 400.0
s_f = solve_ivp(rhs, (0, t_sep), [th_max, 0.0], args=(k05,), method="RK45",
                rtol=1e-12, atol=1e-14, t_eval=np.linspace(0, t_sep, 20000))
s_b = solve_ivp(rhs, (0, -t_sep), [th_max, 0.0], args=(k05,), method="RK45",
                rtol=1e-12, atol=1e-14, t_eval=np.linspace(0, -t_sep, 20000))
sep_right = np.concatenate([s_b.y[:, ::-1], s_f.y], axis=1)
sep_left = np.array([-sep_right[0], sep_right[1]])  # mirror symmetry theta -> -theta
results["separatrix_sim_check"] = {
    "closest_rad_to_saddle_fwd": float(np.min(np.hypot(s_f.y[0], s_f.y[1]))),
    "max_abs_E_minus_Esep_along_sim": float(np.max(np.abs(
        energy(s_f.y[0], s_f.y[1], k05) - E_sep)))}

# ------------------------------------------ Fig 4: separatrix, k = 0.5 k_cr, 4(c)
def sep_speed(x):
    val = 2 * (g / l) * (1 - np.cos(x)) - k05 / (m * l**2) * x**2
    return np.sqrt(np.clip(val, 0, None))

xs = np.linspace(-th_max, th_max, 2001)
Y = sols[0.5]
E = energy(Y[:, 0], Y[:, 1], k05)
sign_E = np.sign(E[:, 0] - E_sep)
# does sign(E - E_sep) ever change along a trajectory? (should not)
sign_changes = int(np.sum([np.any(np.sign(E[i] - E_sep) != sign_E[i]) for i in range(N)]))
# geometric check: inside the loop <=> |thd| < sep_speed(theta) for |theta| < th_max
inside_left = inside_right = outer = 0
geom_violations = 0
margin_min = np.inf
for i in range(N):
    th_i, thd_i = Y[i, 0], Y[i, 1]
    m_in = np.abs(th_i) < th_max
    if E[i, 0] < E_sep:
        if th_i[0] < 0: inside_left += 1
        else: inside_right += 1
        # must stay strictly inside loop and on the same side of theta=0
        geom_violations += int(np.any(np.abs(thd_i) >= sep_speed(th_i) + 1e-9) or np.any(np.sign(th_i) != np.sign(th_i[0])))
    else:
        outer += 1
        # must stay strictly outside: |thd| > sep_speed wherever the separatrix exists
        geom_violations += int(np.any(np.abs(thd_i[m_in]) <= sep_speed(th_i[m_in]) - 1e-9))
    margin_min = min(margin_min, float(np.min(np.abs(E[i] - E_sep))))
results["separatrix_check"] = {"left_well": inside_left, "right_well": inside_right,
                               "outer": outer, "sign_changes": sign_changes,
                               "geom_violations": int(geom_violations),
                               "min_abs_E_minus_Esep": margin_min}

fig, ax = plt.subplots(figsize=(6.4, 4.8))
for i in range(N):
    ax.plot(Y[i, 0], Y[i, 1], lw=0.9, color=cmap(i / (N - 1)), alpha=0.9)
    ax.plot(Y[i, 0, 0], Y[i, 1, 0], ".", color="k", ms=4)
# separatrix obtained by directly simulating a trajectory released from rest
# at the turning point theta = +- theta_max (Sec. 4(c) hint)
ax.plot(sep_right[0], sep_right[1], "k--", lw=2.2, label="separatrix (simulated)")
ax.plot(sep_left[0], sep_left[1], "k--", lw=2.2)
eq_markers(ax, k05)
ax.set_xlabel(r"$\theta$ [rad]"); ax.set_ylabel(r"$\dot\theta$ [rad/s]")
ax.set_title(r"Separatrix, $k=0.5\,k_{cr}$")
ax.legend(handles=[Line2D([], [], color="k", ls="--", lw=2.2, label=r"separatrix (simulated, $E=mg\ell$)")] + handles[:2],
          loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8, framealpha=0.95)
ax.set_ylim(-3.6, 3.6)
fig.tight_layout(); fig.savefig("fig_separatrix.pdf"); plt.close(fig)

with open("results.json", "w") as f:
    json.dump(results, f, indent=2)
print(json.dumps(results, indent=2))
