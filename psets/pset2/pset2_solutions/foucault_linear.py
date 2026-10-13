"""16.07 Problem Set 2 -- Foucault pendulum simulation (solution code).

Integrates the linearized equations of motion from Section 1(c),
    M qdd + C qd + K q = 0,   q = [x, y]  (east, north),
with M = m I, C = 2 m Omega sin(lat) J, K = (m g / ell) I, and compares the
precession of the swing plane with the prediction psi_dot = -Omega sin(lat).
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ---------------- Parameters ----------------
m = 1.0                        # kg (cancels out of the dynamics)
g = 9.81                       # m/s^2
ell = 10.0                     # m
Omega = 7.2921e-5              # rad/s (sidereal)
lat = np.deg2rad(42.36)        # MIT
Ou = Omega * np.sin(lat)       # vertical component of Earth's rotation
w0sq = g / ell

psi_dot_pred = -Ou
T_prec = 2 * np.pi / Ou        # time for the swing plane to turn 360 deg
T_final = 36 * 3600.0          # > one full precession period (35.5 h)

J = np.array([[0.0, -1.0], [1.0, 0.0]])
M = m * np.eye(2)
C = 2 * m * Ou * J
K = m * w0sq * np.eye(2)
Minv = np.linalg.inv(M)


def dynamics(t, X):
    q, qd = X[:2], X[2:]
    qdd = -Minv @ (C @ qd + K @ q)
    return np.concatenate([qd, qdd])


def energy(X):
    """Earth-frame energy E = 1/2 qd^T M qd + 1/2 q^T K q (conserved since C = -C^T)."""
    q, qd = X[:2], X[2:]
    return 0.5 * np.einsum("i...,ij,j...->...", qd, M, qd) + \
           0.5 * np.einsum("i...,ij,j...->...", q, K, q)


def swing_angle(t, x, y):
    """Swing-plane angle psi at the turning points (local maxima of r)."""
    r = np.hypot(x, y)
    k = np.where((r[1:-1] > r[:-2]) & (r[1:-1] >= r[2:]))[0] + 1
    psi = 0.5 * np.arctan2(2 * x[k] * y[k], x[k]**2 - y[k]**2)  # angle of the line, mod pi
    psi = 0.5 * np.unwrap(2 * psi)                             # remove the pi jumps
    return t[k], psi


if __name__ == "__main__":
    theta0 = np.deg2rad(5.0)
    X0 = [ell * np.sin(theta0), 0.0, 0.0, 0.0]   # released from rest, displaced east
    t_eval = np.arange(0.0, T_final, 0.05)
    sol = solve_ivp(dynamics, [0, T_final], X0, method="DOP853",
                    t_eval=t_eval, rtol=1e-10, atol=1e-10)

    # Energy check
    E = energy(sol.y)
    print(f"max relative energy drift: {np.max(np.abs(E - E[0])) / E[0]:.2e}")

    # Precession rate
    tk, psi = swing_angle(sol.t, sol.y[0], sol.y[1])
    rate = np.polyfit(tk, psi, 1)[0]
    err = 100 * (rate - psi_dot_pred) / psi_dot_pred
    print(f"predicted: {np.rad2deg(psi_dot_pred) * 3600:.6f} deg/h, "
          f"period {T_prec / 3600:.3f} h")
    print(f"measured:  {np.rad2deg(rate) * 3600:.6f} deg/h, error {err:.2e} %")

    # Plot
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(tk / 3600, np.rad2deg(psi), "C0", lw=2.5, label="simulation")
    ax.plot(tk / 3600, np.rad2deg(psi_dot_pred * tk), "k--", lw=1.2,
            label=r"prediction $-\Omega\sin\lambda\,t$")
    ax.axvline(T_prec / 3600, color="gray", lw=0.8, ls=":")
    ax.text(T_prec / 3600 - 0.4, -40, "one precession\nperiod (35.5 h)",
            ha="right", va="top", fontsize=9, color="gray")
    ax.set_xlabel("time [h]")
    ax.set_ylabel(r"swing-plane angle $\psi$ [deg]")
    ax.set_yticks(np.arange(0, -361, -60))
    ax.set_title(r"Linearized Foucault pendulum, $\lambda = 42.36^\circ$N, $\ell = 10$ m")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig("psi_linear.png", dpi=200)
