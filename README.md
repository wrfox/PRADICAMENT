# PRADICAMENT
Data analysis tools for electromagnetic field measurement by particle deflectometry in plasmas.


**PRADICAMENT** [1] is a set of  [2,3] tools for analyzing proton radiography (PRAD) fluence images.
The present scope of the project is a fast algorithm to 
process 1-D fluence profiles to obtain the line-integrated magnetic or electric field.
A companion manuscript [1] describes the theory, setup, and notation used here in detail.
The routines are currently written in Matlab.  The 1-D algorithm is complementary to existing 2-D algorithms [4].



References:

[1] W. Fox et al, Proton deflectometry analysis in magnetized plasmas: magnetic field reconstruction in 1-D.  Submitted (2023)

[2] [N. L. Kugland, et al. Invited Article: Relation between electric and magnetic field structures and their proton-beam images. Review of Scientific Instruments 83, 101301 (2012).](https:/doi.org/10.1063/1.4750234)

[3] [D. B. Schaeffer, et al,  Proton imaging of high-energy-density laboratory plasmas. (2023).](https://arxiv.org/abs/2212.08252)

[4] [Bott, A. F. A. et al. Proton imaging of stochastic magnetic fields. Journal of Plasma Physics 83, (2017).](https://doi.org/10.1017/S0022377817000939);  Implementation: https://github.com/flash-center/PROBLEM


# Installation

`git clone https://github.com/wrfox/PRADICAMENT.git`

Clone the package, and add the directory to your Matlab path.


# Problem statement:

Proton radiography (or deflectometry) is commonly used in
plasma physics to observe electromagnetic fields.
The principle of the measurement is that a point source 
of high-energy protons or other charged particles is generated 
and sent through an experimental volume,
where they pick up deflections from the 
electromagnetic fields.   For high energy particles,
the deflections are small-angle and proportional to the line integral
of the electromagnetic fields along the particle trajectories.
The deflections produce a non-uniform proton fluence pattern
on the detector.  The goal of fluence
analysis is to recover the electromagnetic fields which
produce an observed fluence pattern.

**Notation**  The plasma plane is assumed to lie
between the proton source and detector, 
a distance $L_s$ from the proton point source and
$L_d$ from the detector, giving a magnification $M = (1 + L_d/L_s)$.
The routines here work in the coordinate system of 
the plasma-plane, mapping coordinates from the detector back 
to the plasma plane using $1/M$.
In this coordinate system, the effect of the deflection by the electromagnetic fields 
can be regarded as a mapping from $x \to x'$, where 
$x$ is the initial proton position (where a given proton crosses the plasma),
and $x'$ the final proton position, defined by 

$x' = x + \xi(x)$.

Following the equations in [1], $\xi$ can be related to plasma
electromagnetic fields according to 

$\xi(x)  = (1/K_B) \int d\ell \times \mathbf{B} + (1/K_E) \int \mathbf{E} _\perp d\ell $ 

Here $K_B$ and $K_E$ are constants that relate the proton deflection (in plasma-plane coordinates) to the magnetic and electric fields.  In SI units they are given by 

$K_B  = (m_p V_p / e) (L_s + L_d) / L_s  L_d$

and 

$K_E  = (m_p V_p^2 / e) (L_s + L_d) / L_s  L_d$

where $V_p$ is the proton speed, $m_p$ the proton mass, and $e$ the fundamental charge unit.  $K_B$ has units of `T`,
and $K_E$ has units of `V/m`.

For a given set of electromagnetic fields, it is easy to imagine calculating the forward model
proton fluence image:  we just numerically calculate many protons mapped according to $x \to x + \xi(x)$, and 
bin the results.  (This is the routine `prad_fwd`).  
 
**The core routines `prad_inv` and `prad_inv_bc` do the inverse problem (in 1-D), which is to recover $\int d\ell \times B$ given a proton image $I(x)$**

A final required input for all routines is the "initial" proton fluence $I_0(x)$, which is the proton fluence *before* it reaches the plasma.  This important quantity is discussed in publication references.  For the purposes of this documentation
it is taken as a known quantity.


# Notation in routines

$x$ = spatial coordinate, using the coordinates of the plasma plane

$b(x)$ = line integrated magnetic field variable $= \int d\ell \times B$.

$I(x)$ = observed proton fluence (reaching the detector, in coordinates of the plasma plane)

$I_0(x)$ = initial proton fluence (before reaching plasma plane).

$K_B$ = deflection parameter.

In the routines, $x$ is a vector, and $b$, $I$, and $I_0$ are expected to be vectors of the same shape.  Optionally,
$I_0$ can be a scalar, in which case it is taken as uniform across the region.

The routines work in any consistent unit scheme.  This means, given a spatial unit for positions $x$ (i.e. `m`), and $K_B$ calculated in `T`, then the units of $b$ will be `T-m`.  However, $K_B$ could also be converted e.g. to Gauss, and other spatial units can be used.  Be careful with the above equation for $K_B$ however, as that is specifically in SI units.

In general an experimental analysis must discriminate between electric and magnetic deflections, however this must be done
by a higher-level analysis outside of the these routines, which just do the reconstruction. 
(Examples include comparing reconstructions at multiple energies or date taken from various orientations.)
The routines below default to determining the line-integrated magnetic field $b(x)$, *assuming the electric deflection
is negligible*.  They can equally determine the electric deflection (assuming negligible B field)
by passing $K_E$ instead of $K_B$ as a parameter.
Finally, if unity ("1") is passed for $K_B$, then the routines return the deflection field $\xi$,
in the same unit as the spatial coordinate $x$, which
could then be used for further analysis separating E and B deflections.

# Routines

Additional documentation is available via inline Matlab help and comments in the files.

```
[x,I] = prad_fwd(x, b, KB, I0)
```

* `prad_fwd` is the proton forward model in 1-D.  Given a magnetic field profile `b(x)`, and input proton fluence `I0(x)`, generate
 the forward model observed fluence `I(x)`.

```
[x,b] = prad_inv(x, I, I0, KB, [optional x0,b0] )
```

* `prad_inv.m` is the proton inverse solver:  Given proton fluence data, obtain the magnetic fields.  Optionally specify a known
magnetic field boundary conditions `b=b0` at `x=x0`.
  

```
[x,b] = prad_inv_bc(x, I, I0, KB, x_bc, b_bc)
```
* `prad_inv_bc.m` is the proton inverse solver with boundary conditions.   It uses *pair* of boundary conditions on the magnetic field, `(x1, b1), (x2,b2)` where `x_bc = [x1,x2]` and `b_bc = [b1,b2]` and generates a solution which crosses through these points.  To do so, it renormalizes the input fluence $I_0$ to achieve the boundary condition


# Contributions

The algorithm was developed by G. Fiksel and W. Fox.  

Community contributions are welcome.
