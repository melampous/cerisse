# Revised section: Nozzle exit profiles (matches the computed matrix)

A prescribed finite-thickness velocity profile at the nozzle exit is adopted as
the primary inlet condition. Within the circular orifice the axial velocity is
defined by the shifted hyperbolic-tangent function

$$
f(r) \;=\; \tfrac12\!\left[1-\tanh\!\left(\frac{r-r_0}{\delta}\right)\right],
\qquad r_0 = R_e - 3\delta,
\tag{3.13}
$$

where $R_e = D_e/2$ and, for the baseline case,

$$
\delta = 0.01\,D_e = 254\ \mu\mathrm{m}.
$$

The velocity components are prescribed as $u_x(r)=u_e f(r)$, $u_\perp=0$.
Locating the centre of the transition at $R_e-3\delta$ reduces the axial
velocity to approximately $0.0025\,u_e$ at the orifice edge.

The static pressure is uniform across the orifice. Temperature and density
follow from a constant total temperature:

$$
p(r)=p_e,\qquad
T(r)=T_0-\frac{u_x^2(r)}{2c_p},\qquad
\rho(r)=\frac{p_e}{R\,T(r)}.
$$

Thus only the axial velocity is prescribed through $f(r)$; density,
temperature and pressure are not blended with the ambient state. The
profile-integrated effective discharge coefficient relative to an ideal
top-hat exit is

$$
C_{d,\mathrm{eff}}=\frac{\dot m}{\rho_e u_e \pi R_e^2}\simeq 0.879 .
$$

This value is a consequence of the prescribed profile and is not an
experimentally measured discharge coefficient. On the upstream boundary the
profile is applied only for $r<R_e$; the region $R_e\le r\le 6D_e$ is an
adiabatic slip wall representing the simplified flush baffle, and $r>6D_e$
carries the quiescent ambient state.

## Exit-profile sensitivity matrix

The sensitivity of the computed near field to the prescribed exit profile is
assessed with three additional cases, chosen to separate the two candidate
control parameters of the shock-cell structure — the integral mass/momentum
deficit of the profile and its local shear-layer scale.

**(i) Top hat.** The limiting case of a discontinuous velocity change at the
lip,

$$
f_{\mathrm{TH}}(r)=\begin{cases}1, & r<R_e,\\ 0, & r\ge R_e,\end{cases}
$$

with $C_{d,\mathrm{eff}}=1$ by construction.

**(ii) Thin shifted-tangent profile.** Equation (3.13) with
$\delta=100\ \mu\mathrm{m}$ ($\delta/D_e\simeq 0.004$), i.e. the same profile
family as the baseline at reduced thickness.

**(iii) Wall-attached tangent profile.** A profile whose maximum velocity
gradient lies at the lip itself,

$$
f_{\mathrm{WT}}(r)=\tanh\!\left(\frac{R_e-r}{a}\right),\qquad
a=200\ \mu\mathrm{m},
$$

which resembles an attached exit boundary layer and, unlike Eq. (3.13),
produces no near-stagnant annulus adjacent to the lip.

Cases (ii) and (iii) are constructed as a matched pair: with the vorticity
thickness defined as
$\delta_\omega=\Delta U/\max|\partial u/\partial r|$, one obtains
$\delta_\omega=2\delta=200\ \mu\mathrm{m}$ for case (ii) and
$\delta_\omega=a=200\ \mu\mathrm{m}$ for case (iii); both profiles also share
$\delta_{99}\simeq 530\ \mu\mathrm{m}$. The two cases nevertheless differ in
half-velocity location, displacement thickness, shape factor and integral
fluxes (table below), so that agreement between them would indicate control
by the local shear scale, whereas ordering by $C_{d,\mathrm{eff}}$ would
indicate control by the integral deficit.

| case | profile | $\delta$ or $a$ [µm] | $\delta_\omega$ [µm] | $y_{50}$ [µm] | $\delta^*$ [µm] | $\theta$ [µm] | $H$ | $C_{d,\mathrm{eff}}$ | $C_{M}$ |
|---|---|---|---|---|---|---|---|---|---|
| top hat | $f_{\mathrm{TH}}$ | — | — | 0 | 0 | 0 | — | 1.000 | 1.000 |
| wall-tanh | $f_{\mathrm{WT}}$ | $a=200$ | 200 | 110 | 139 | 61 | 2.26 | 0.9755 | 0.9669 |
| shifted-tanh | Eq. (3.13) | $\delta=100$ | 200 | 300 | 300 | 50 | 6.02 | 0.9513 | 0.9445 |
| baseline | Eq. (3.13) | $\delta=254$ | 508 | 762 | 762 | 127 | 6.02 | 0.8789 | 0.8623 |

For the shifted-tangent family the nominal momentum thickness is

$$
\theta_0=\int_0^{R_e} f(r)\,[1-f(r)]\,\mathrm{d}r \simeq \frac{\delta}{2},
\tag{3.14}
$$

when the transition is contained within the orifice; for the wall-attached
profile $\theta_0=a(1-\ln 2)\simeq 0.31\,a$.

A discrete-quadrature audit of the inflow boundary confirms that the mass and
momentum fluxes of all four profiles are represented to within 0.01 % on the
production grid (cell-centred sampling at $\Delta x = 99\ \mu\mathrm{m}$),
while the maximum velocity gradient of the 200-µm-scale profiles is resolved
by only $\approx 2$ cells; conclusions drawn from the matched pair are
therefore restricted to integral and mean-cell metrics (§X).

The current profile and the computed variations are shown in Fig. 3.3. The
Pohlhausen–Blasius and one-seventh-power profiles are retained in the figure
for visual comparison only and are not part of the simulation matrix.
