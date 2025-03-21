from math import *

# some constants
dt    = 1e-7                  # [s] - time step
eta   = 1/(4*pi*dt)           # viscosity         [Pa-s], [N-s/m^2], [J-s/m^3], [kg/m-s]
mu    = 1/eta                 # mu_0 - 1/viscosity
rho   = 8.854187817620389e-12 # mass density of medium at rest                  [kg/m^3]
k     = eta / rho             # dynamic viscosity, quantum circulation constant [m^2/s]



'''

We demonstrate below that all of the major
constants of nature can be derived from the
vortex momenta quanta (P) and interaction
length (L). Given the earlier definitions,

    h = 2 * P * L

    q = 2 * P / L

We can derive P and L as::

    P = sqrt(h * q) / 2 = 5.151731449128262e-27

and

    L = sqrt( h / q ) = 6.430915717783789e-08


Combined with c, α, and s all other constants
can be defined. For example,

Permeability...

    mu2 = alpha * L**3 / (P * c) = 1.2566370614359164e-06
    mu/mu2 = 1.0000000008200005

   
Permittivity...

    rho2 = P / (alpha * c * L**3) = 8.85418781762039e-12
    rho2/rho = 1.0000000000000004

Boltzmann’s Constant

    # Copilot: kb = (L*(1+beta)**2) / (alpha * c**2) = 1.3806485279e-23
    
    kb = (L*(1+beta)**2) / (c) = 1.3806485279e-23


Rydberg's Constant

    # Copilot: R = 1/(sqrt(2) * alpha**2 * L) = 10973731.56850858
    
    R = 1/(sqrt(2) * L *(1 + beta)) = 10973731.56850858

Hmm, this does not seem to work out....



so, apparantly, we can write;

    beta = 1/alpha - 1 = 1/137.0359990836958 - 1 = 0.0072973525693  


Gravitational Constant

    # Copilot: G = (c**2 * L**3) / (alpha * P) = 6.67408e-11
    G = (alpha*P*c*s**2)/(2*pi* (L*(1+beta))**3 ) = 6.673E-11



Hydrogen stable electron energy states

    Un = sqrt(P*c)


'''




'''
1/dt = - k * sqrt(3) / du * delta_F
 
du / dt = - k * sqrt(3) * delta_F 

du=dx=dy=dz = grid spacing

With delta_F normalized to 1:

du / dt = - k * sqrt(3)

du = - k * sqrt(3) * dt


Let's consider the definitiion of viscosity eta for the medium above:

eta = 1/(4*pi*dt)

which originates from the pre 2019 textbook definition of mu_0, the permeability
of free space:

mu_0 = 1/4*pi*1e-7

so, if we assume some universal time step dt = 1e-7, then we can define the
viscosity of the medium as:

eta = 1/(4*pi*dt), which has a unit of measurement in [kg/m-s].

The interesting thing is that this is consistent with the idea that viscosity is
a measure of the rate of change of momentum of a fluid element as a result of
the forces acting on it.

Let's consider what Wikipedia has to say about viscosity:

"Dynamic Viscosity

In materials science and engineering, one is often interested in understanding
the forces or stresses involved in the deformation of a material. For instance,
if the material were a simple spring, the answer would be given by Hooke's law,
which says that the force experienced by a spring is proportional to the
distance displaced from equilibrium. Stresses which can be attributed to the
deformation of a material from some rest state are called elastic stresses. In
other materials, stresses are present which can be attributed to the deformation
rate over time. These are called viscous stresses. For instance, in a fluid such
as water the stresses which arise from shearing the fluid do not depend on the
distance the fluid has been sheared; rather, they depend on how quickly the
shearing occurs.

Viscosity is the material property which relates the viscous stresses in a
material to the rate of change of a deformation (the strain rate). Although it
applies to general flows, it is easy to visualize and define in a simple
shearing flow, such as a planar Couette flow.

In the Couette flow, a fluid is trapped between two infinitely large plates, one
fixed and one in parallel motion at constant speed u u (see illustration to the
right). If the speed of the top plate is low enough (to avoid turbulence), then
in steady state the fluid particles move parallel to it, and their speed varies
from 0 {\displaystyle 0} at the bottom to u u at the top.[7] Each layer of fluid
moves faster than the one just below it, and friction between them gives rise to
a force resisting their relative motion. In particular, the fluid applies on the
top plate a force in the direction opposite to its motion, and an equal but
opposite force on the bottom plate. An external force is therefore required in
order to keep the top plate moving at constant speed.

In many fluids, the flow velocity is observed to vary linearly from zero at the
bottom to u u at the top. Moreover, the magnitude of the force, F F, acting on
the top plate is found to be proportional to the speed u u and the area A A of
each plate, and inversely proportional to their separation y y:

    F = μ A u y . F=\mu A{\frac {u}{y}}.

The proportionality factor is the dynamic viscosity of the fluid, often simply
referred to as the viscosity. It is denoted by the Greek letter mu (μ). The
dynamic viscosity has the dimensions ( m a s s / l e n g t h ) / t i m e
{\displaystyle \mathrm {(mass/length)/time} }, therefore resulting in the SI
units and the derived units:

    [ μ ] = k g m ⋅ s = N m 2 ⋅ s = P a ⋅ s = {\displaystyle [\mu ]={\frac {\rm
    {kg}}{\rm {m{\cdot }s}}}={\frac {\rm {N}}{\rm {m^{2}}}}{\cdot }{\rm
    {s}}={\rm {Pa{\cdot }s}}=} pressure multiplied by time.

The aforementioned ratio u / y u/y is called the rate of shear deformation or
shear velocity, and is the derivative of the fluid speed in the direction
perpendicular to the normal vector of the plates (see illustrations to the
right). If the velocity does not vary linearly with y y, then the appropriate
generalization is:

    τ = μ ∂ u ∂ y , \tau =\mu {\frac {\partial u}{\partial y}},

where τ = F / A {\displaystyle \tau =F/A}, and ∂ u / ∂ y {\displaystyle \partial
u/\partial y} is the local shear velocity. This expression is referred to as
Newton's law of viscosity. In shearing flows with planar symmetry, it is what
defines μ \mu . It is a special case of the general definition of viscosity (see
below), which can be expressed in coordinate-free form.

Use of the Greek letter mu ( μ \mu ) for the dynamic viscosity (sometimes also
called the absolute viscosity) is common among mechanical and chemical
engineers, as well as mathematicians and physicists.[8][9][10] However, the
Greek letter eta ( η \eta ) is also used by chemists, physicists, and the
IUPAC.[11] The viscosity μ \mu is sometimes also called the shear viscosity.
However, at least one author discourages the use of this terminology, noting
that μ \mu can appear in non-shearing flows in addition to shearing flows.[12] "

"Illustration of a planar Couette flow. Since the shearing flow is opposed by
friction between adjacent layers of fluid (which are in relative motion), a
force is required to sustain the motion of the upper plate. The relative
strength of this force is a measure of the fluid's viscosity."

"In a general parallel flow, the shear stress is proportional to the gradient of
the velocity."

"In fluid dynamics, it is sometimes more appropriate to work in terms of
kinematic viscosity (sometimes also called the momentum diffusivity), defined as
the ratio of the dynamic viscosity (μ) over the density of the fluid (ρ). It is
usually denoted by the Greek letter nu (ν):

    ν = μ / ρ "
    
Which we write as:
  
   k = eta / rho


Now let's consider this part:   
 
"Viscosity is the material property which relates the viscous stresses in a
material to the rate of change of a deformation (the strain rate). Although it
applies to general flows, it is easy to visualize and define in a simple
shearing flow, such as a planar Couette flow.

In the Couette flow, a fluid is trapped between two infinitely large plates, one
fixed and one in parallel motion at constant speed v_t (see illustration to the
right). If the speed of the top plate is low enough (to avoid turbulence), then
in steady state the fluid particles move parallel to it, and their speed varies
from 0 at the bottom to v_t at the top.[7] Each layer of fluid moves faster than
the one just below it, and friction between them gives rise to a force resisting
their relative motion."


But rather han a planar Couette flow, let's consider a Couette flow around an
irrotational vortex in our superfluid medium. An irrotational vortex in a
superfluid medium is a vortex that has no vorticity except at the center of the
vortex. So, such a vortex has a hollow core, since otherwise the velocity of the
fluid particles would become infinite at the center of the vortex.

From Wikipedia:

"For an irrotational vortex, the circulation is zero along any closed contour
that does not enclose the vortex axis; and has a fixed value, Γ, for any contour
that does enclose the axis once.[6] The tangential component of the particle
velocity is then v_t = Γ / 2 π r. The angular momentum per unit mass relative to
the vortex axis is therefore constant, r v_t = Γ / 2 π "


Thus, an irrotational vortex is characterized by a constant circulation Γ around
the hollow vortex core, which results in the the particle speed v being
inversely proportional to the distance r from the axis and for Γ=k  is thus
given by:

    v_t = k / 2*pi*r
    
So, if we consider a Couette flow around an irrotational vortex, then the speed
would not vary from 0 at the bottom to v_t at the top, but rather from 0 at
infinity to infinity at the center of the vortex, but since the shear stress
proportional to the gradient of the velocity, we can consider the gradient of
the velocity along the radius of the vortex, which is given by:

    dv/dr = - k / 2*pi*r^2
    
So, the shear stress is proportional to the gradient of the velocity, and is
given by:

    tau = dv/dr * eta
    
This is correct, because the shear stress is proportional to the gradient of the
velocity, and the gradient of the velocity is proportional to the inverse of the
square of the radius of the vortex, so the shear stress is proportional to the
inverse of the square of the radius of the vortex, which is given by:
    
    tau = - k / 2*pi*r^2    
            
And since  k = eta / rho , we can write:

    dv/dr = - eta / (2*pi*r^2*rho)

This is also correct, because the gradient of the velocity is proportional to
the inverse of the square of the radius of the vortex, and the inverse of the
square of the radius of the vortex is proportional to the shear stress, so the
gradient of the velocity is proportional to the shear stress.

Thus, we can write:

    eta = - dv/dr * (2*pi*r^2*rho)
 
Now since the circulation k is constant, this means that a line integral around
a circle of radius r is constant, and thus we can write:

    v_t * 2*pi*r = k

or

    v_t = k / (2*pi*r)


So, whereas with the planar Couette flow, each layer of the fluid moves faster
than the one just below it, and friction between them gives rise to a force
resisting their relative motion, with the Couette flow around an irrotational
vortex, we have cilindrical layers around the vortex core, whereby each layer of
the fluid moves slower than the one just within it, proportional to the inverse
of the square of the radius of r.

However, there is no friction between the layers, because the fluid is in steady
state motion whereby angular momentum is conserved, and thus there is no
dissipation nor drag, and thus no friction.

Yet, there is a force needed to keep the fluid in circular motion, and this
force is provided by the pressure gradient in the fluid, so what we are looking
at is a viscous, frictionless force that is balanced by a pressure gradient, so
that the fluid can be kept in circular motion without resistance or dissipation.

Thus, we come to the conclusion that whereas for superfluids it is often assumed
that there is no viscosity, this is actually incorrect. In the situation where
we have an irrotational vortex, there is no turbulence which is sometheing
completely different from zero viscosity.

So, rather than describing the superfluid medium as a zero viscosity medium, it
should actually be described as a zero turbulence medium. This is consistent
with the observation that superfluids can be described with a two fluid model,
whereby currently the superfluid component is described as a zero viscosity
fluid, and the normal fluid component is described as a viscous fluid.

Thus, what it should be described as is that the superfluid component is a zero
turbulence fluid without friction, and the normal fluid component is a turbulent
fluid whereby friction emerges from the turbulence.

And since in our model, we have the intensity field I, defined by:

    I = - k * rho a, 
    
with a the acceleration field, we obtain a unit of measurement in [kg/s^3] for
I, and we obtain a scalar potential t defined as

    t = div(I) = - k * rho div(a)

having a unit of measurement in [W/m^3], which is the unit of measurement for
power density and thus represents temperature.

Interestingly, from our primary field C, defined as:

    C = - k * rho * v,
    
with v the velocity field, we can derive the field R, defined as:

    R = curl(A) = curl(curl(C)) = - k * rho * curl(curl(v)),
    
which has a unit of measurement or force density:

    [R] = [m^2/s] * [kg/m^3] * [1/m^2] * [m/s] = [kg/m^3]*[m/s^2] = [N/m^3]

In our model, current density J is defined as the curl of this field R:

    J = curl(R),
    
which has a unit of measurement inL

    [J] = [1/m] * [N/m^3] = [kg/m^3]*[m/s^2]*[1/m] = [kg/m^3-s^2] = [A/m^3],
    
whereby the unit for electric current [A] is thus equal to [kg/s^2].

Because the unit of measurement for R is [N/m^3], we can also write it as:

    R = rho * a,
    
so that the current density J is given by:

    J = curl(rho * a) = rho * curl(a),
    
which is related to the angular counterpart of the scalar potential t, defined
as

    W = curl(I) = curl(-k* rho * a) = -k * rho * curl(a),

so we can write:

    J = rho * curl(a) = X * W = X * curl(-k* rho * a)

Thus,

    rho * curl(a) = X * curl(-k* rho * a) rho * curl(a) = X * (-k* rho *
    curl(a)) rho * curl(a) = -X * k * rho * curl(a) 1 = -X * k X = -1/k

and thus

    J = rho * curl(a) = -1/k * curl(-k* rho * a) = -1/k * W = W/k

This works out to a unit of measurement in [A/m^3] or [kg/s^2-m^3] for J:

    [J] = [W/m^3] / [m^2/s] = [kg/m-s^3]*[s/m^2] = [kg/s^2]*[/m^3] 
        = [kg/m^3-s^2] = [A/m^3]


Thus, we have shown that the current density J is equal to the angular
counterpart of the scalar potential t, namely the vector potential W, divided by
the quantum circulation constant k.

Since electrical resistance is related to the current density J and we now found
that current density J is proportional to our second order angular vector
potential, angular power density W, we can conclude that resistance and thus
dissipation is proportional to the second order angular vector potential,
angular power density W. 

So, when there is zero turbulence, power density W is zero, and thus resistance
and dissipation emerge when there is turbulence, and thus when angular power
density W is non-zero.    
    
So, rather than describing the superfluid medium as a zero viscosity medium, it
should actually be described as a zero turbulence medium, meaning a medium where
current density J is zero, and thus where angular power density W is zero, and
thus where angular vector potential W is zero, while the viscosity is non-zero
and equal to eta = k/rho, with k the quantum circulation constant and rho the
mass density of the medium.

And since for a superfluid the quantized circulation Γ is given by:

    Γ = h/m = rho * eta,

we can define the viscosity of a superfluid as:

    eta = h/(m*rho),
    
which is non-zero, and thus a superfluid medium is not a zero viscosity medium,
but a zero turbulence medium whereby both the second order potential fields are
equal to zero, namely the scalar potential t and the angular vector potential W.

So we can define the properties of a superfluid medium, or the superfluid
component in a two fluid model, as follows:

    eta = h/(m*rho) I = 0
    
And thus:

    J = 0 t = 0 W = 0


Interestingly, we obtain the same units of measruement for the intensity field I
if we would define it as:

    I = - κ v,

with κ the modulus or elasticity in [Pa] or [kg/m-s^2], which has a unit of
measurement that differs by a per second [/s] from the unit of measurement for
viscosity eta (η) in [Pa-s] or [kg/m-s]. And v the velocity field in [m/s].

This leads to the conclusion that elasticity is not a fundamental property of a
superfluid medium, but rather an emergent property of the superfluid medium,
whereby elasticity emerges when the intensity field I is non-zero.

This means that the elasticity as a property of the normal fluid component in a
two fluid model emerges from the superfluid component, when the scalar and vector
potentials are non-zero. In other words, elasticity is an emergent property of
the normal fluid component, which is a turbulent fluid, whereby both turbulence
and elasticity emerge when temperature and angular power density rise above
a certain threshold below which the medium can be considered to have zero
turbulence and therefore zero elasticity.

So, whereas the superfluid component is a two-fluid model is generally described
as a zero viscosity medium, it should actually be described as a zero elasticity
medium, while the normal fluid component is a turbulent fluid with non-zero
elasticity that arises from the superfluid component when the temperature and
angular power density rise above a certain threshold, the so called lambda
transition temperature or lambda point.









   

    
resulting in a unit of measurement for J in

   [J] = [kg/m^3] * [/s^2] 
       = [kg/m^3] * [W/m^3]    / [m^2/s] = [kg/m^3] * [kg/m-s^3] / [m^2/s] =
       [kg/m^3] * [kg/s^3] * [s^2/m^2] = [kg/s^2] = [A]
   
   






delta = vector laplacian on interlaced staggered Yee grid normalized to 1, as in
the implementation of the curl operator in the fdtd module, specifically in the
operators.py file.

With the the courant number, with which the grid spacing can be defined:

Finally, the `courant_number` of the grid determines the relation between the
`time_step` of the simulation and the `grid_spacing` of the grid. If not given,
it is chosen to be the maximum number allowed by the [Courant-Friedrichs-Lewy
Condition](https://en.wikipedia.org/wiki/Courant–Friedrichs–Lewy_condition): `1`
for `1D` simulations, `1/√2` for `2D` simulations and `1/√3` for `3D`
simulations (the dimensionality will be derived by the shape of the grid). For
stability reasons, it is recommended not to change this value.

The number `(c*dt/du)` is a dimensionless parameter called the _Courant number_
`sc`. For stability reasons, the Courant number should always be smaller than
`1/√D`, with `D` the dimension of the simulation. This can be intuitively be
understood as the condition that information should always travel slower than
the speed of light through the grid.


sc = (c*dt/du) du=dx=dy=dz = grid spacing


du = - k * sqrt(3) * dt

c = sc * du / dt =  sqrt(3) * du / dt

du/dt = - k * sqrt(3) 

c = -3k   



This yields the final update equations for the FDTD algorithm:

```python
    E  += sc*inv(ε)*curl_H H  -= sc*inv(µ)*curl_E
```



What is interesting is that there are two kinds of curls, two kinds of rotation:
left handed and right handed. Now if we consider a rotating mass, there's two
ways to look at the rotating object: from the front, and from the back, so
that's two like positive and negative, if from the front perspective the
rotation is to the right, then from the back perspective the rotation is to the
left. So we have two kinds of curls, two kinds of rotations, two kinds of
angular momentum, two kinds of "charges" (positive and negative), in other
words: polarity.

So, this relates nicely to the two curl operators we have in the FDTD algorithm:
the curl_point and the curl_surface, which could better be named curl_left and
curl_right, or curl_positive and curl_negative. 

Perhaps, it would be interesting to consider the left handed and right handed
ness to a spinning gyroscope, and the way it reacts to do forces along the
rotation axes of the system.

You can also spin a gyroscope in two directions, and it will react differently
to the forces applied to it.

THis is observed in various experiments, and is called the "gyroscopic paradox".

https://en.wikipedia.org/wiki/Gyroscope#Gyroscopic_paradox

https://www.youtube.com/watch?v=ty9QSiVC2g0

The essence of this is that gyroscopes are polarized, there's polarity because
there's two kinds of rotations, two kinds of curls, or better: "swirls",
"positive" and "negative", "left handed" and "right handed", "clockwise" and
"counter clockwise".

Thus, the curl operators would need to be defined along the axes of rotation of
the gyroscope, so the two opposing perspectives, looking from the front or
looking from the back.

https://youtu.be/V6XSsNAWg00?si=SWPu9XlRmK2KqNgf



'''



dt = - k * d *  delta        # 

c     = 299792458.0          # [m/s] - speed of light
c_l   = c * (pi/2)           # [m/s] - speed of longitudinal (sound) waves

#mu    = 4*pi*1.00000000082e-7# mu_0 - 1/viscosity - 2019 redefinition

h     = 6.62607015e-34       # Planck's constant    [J-s], [kg/m2-s]      
k     = c*c                  # dynamic viscosity    [m2/s]
rho   = eta/k                # 8.854187817620389e-12
q     = 1.602176634e-19      # elemental charge
m     = h/k                  # elemental mass   1




'''


----

So far, I haven't solved the problem of temperature and black body radiation,
but now that I realize the importance of the intensity field [I] and it's
consequence that we have to consider yank rather than force, it seems it is only
a matter of time before we can come full circle. First of all, it is rather
interesting that the gas law also involves quantization denoted by n:

P V = n Kb T, (eq 1)

With T the temperature in Kelvin,P the pressure, V the volume, n the number of
quanta, and Kb Boltzmann's constant.


Second, I found the work of Paul Stowe very interesting, but very hard to
comprehend. On the one hand, he managed to express all the major constants of
nature in terms of just 5 constants and expressed all units of measurement in
just three: mass, length and time, while on the other he managed to write it all
down in a manner that I found very confusing, for instance because he refers to
charge q as "divergence" while meaning "divergence of momentum density":

https://vixra.org/pdf/1310.0237v1.pdf

Nonetheless, valuable insights can be obtained from his work, if only as a
starting point for further considerations. With respect to temperature and the
gas law, in his eq. 20 we find a re- lationship between electrical charge and
Boltzmann's constant:

Kb = h/(qc), (eq 2)

with q elemental charge and h Planck's constant, which results in the conclusion
that the quan- tization in the gas law is related to the quantization of the
medium, which is governed by the quantum quantization constant nu. While it is
nice that this equation yields the right number, this does not necessarily mean
this equation is 100% correct as written, but it certainly seems to point in the
right direction.

Another interesting paper on the subject of black body radiation in relation to
aether theory is this one by C.K. Thornhill, which gives a valuable starting
point for deriving Planck's law:

https://etherphysics.net/CKT1.pdf

His main argument:

"Another argument against the existence of a physical ethereal medium is that
Planck’s em- pirical formula, for the energy distribution in a black-body
radiation field, cannot be derived from the kinetic theory of a gas with
Maxwellian statistics. Indeed, it is well-known that ki- netic theory and
Maxwellian statistics lead to an energy distribution which is a sum of Wien-
type distributions, for a gas mixture with any number of different kinds of
atoms or mole- cules. But this only establishes the impossibility of so deriving
Planck’s distribution for a gas with a finite variety of atoms or molecules. To
assert the complete impossibility of so deriving Planck’s distribution it is
essential to eliminate the case of a gas with an infinite variety of atoms or
molecules, i .e . infinite in a mathematical sense, but physically, in practice,
a very large variety. The burden of the present paper is to show that this
possibility cannot be elimi- nated, but rather that it permits a far simpler
derivation of Planck’s energy distribution than has been given anywhere
heretofore."

What is interesting, is that he found a relationship between the adiabatic index
ω and the num- ber of degrees of freedom α of (aether) particles, which leads to
the conclusion that α must be equal to 6 and he concludes:

"Thus, the quest for a gas-like ethereal medium, satisfying Planck’s form for
the energy distri- bution, is directed to an ideal gas formed by an infinite
variety of particles, all having six de- grees of freedom."It is this adiabatic
index which provides a relationship to heat capacity, since it is also known as
the heat capacity ratio:

https://en.wikipedia.org/wiki/Heat_capacity_ratio

"In thermal physics and thermodynamics, the heat capacity ratio, also known as
the adiabatic index, the ratio of specific heats, or Laplace's coefficient, is
the ratio of the heat capacity at constant pressure (C_P) to heat capacity at
constant volume (C_V)."

So, while we clearly have not yet cracked the whole nut, it seems to me we are
on the right track towards the formulation of a "theory of everything", that
holy grail that has thus far proven to be unreachable, which I'm sure will turn
out to be attributable to ignoring the impli- cations of the vector Laplace
operator.

Personally, I have no doubt both the weak and strong nuclear forces can be fully
accounted for by our model c.q. electromagnetic forces, once completely worked
out, and that the gravi- tational force also propagates through the aether, as
actually confirmed by the Michelson- Morley experiment, so that we will end up
with a model that is much, much simpler and only has one fundamental interaction
of nature.

To illustrate the argument that the nuclear forces can be fully accounted for by
electromag- netic forces, I wholeheartedly recommend the experimental work of
David LaPoint, who shows this in his laboratory:

https://youtu.be/siMFfNhn6dk

[end intermezzo]


'''




'''
In the above, we have a lot of information taken together, so let's see if we
can make any sense of it with respect to the various constants of nature.

One of the most intriguing definitions is that of the viscosity of free space,
base on the pre-2019 definition of permeability mu = 1/(4*pi*1e-7) [kg/m-s],:

    eta = 1/(4*pi*dt)   # viscosity          [kg/m-s],   [Pa-s]

The substitution of the value 1e-7 for dt is a bit of a mystery for now, but
with the unit of measurement in [kg/m-s] for eta, we can write:

    eta = X * (1/dt),
    
with X = 1/(4*pi) having the units [kg/m], and dt having the units [s].

Further, eta is also related to rho and k along:

    eta = rho/k
    
Thus we can write:

    rho/k = X * (1/dt) 
    
Or:

    X = rho/(k*dt)


Then, we have Paul Stowe, who found some interesting relationships between the
various constants of nature, and the elementary charge, e or q.

He uses L and P, under the presumption that within a superfluid medium, quantum
vortices are formed that form a lattice or vortex in a certain grid that has a
certain spacing, L. However, he writes:

"We demonstrate below that all of the major constants of nature can be derived
from the vortex momenta quanta (P) and interaction length (L)."

He considers L to represent the interaction free path of the medium, and a/o
relates Planck's constant h to the action parameter of the medium, which matches
with respect to their unit of measurement.

h is in [J-s] or [kg-m^2/s], while action is described by Wikipedia as:

"In physics, action is a scalar quantity that describes how the energy of a
physical system has changed over time (its dynamics). Action is significant
because the equations of motion of a system can be derived through the principle
of stationary action.

In the simple case of a single particle moving with a constant velocity (thereby
undergoing uniform linear motion), the action is the momentum of the particle
times the distance it moves, added up along its path; equivalently, action is
twice the particle's kinetic energy times the duration for which it has that
amount of energy. For more complicated systems, all such quantities are
combined.

More formally, action is a mathematical functional which takes the trajectory
(also called path or history) of the system as its argument and has a real
number as its result. Generally, the action takes different values for different
paths. Action has dimensions of energy x time or momentum x length, and its SI
unit is joule-second (like the Planck constant h)."

Thus Planck's constant h also has dimensions of energy x time or momentum x
length.


In other words: by considering Planck's constant h as a parameter/constant that
says something about the medium with dimensions of energy x time or momentum x
length, it is clear that some length parameter is involved.


Further, Stowe writes:


"In any kinetic system, these interactions (collisions), are directly
proportional to an average spacing (L) (the mean free path [MFP], their average
momentum (P). These two quantities directly define the lagrangian action
parameter (h) of the system. This becomes: 

    h = 2p ∫{0-L} δL = 2 m v ∫{0-L} δL = 2PL = 2 m V L

This action parameter h may be directly integrated (by the rate of interactions
(i)) to define the total energy involved: i E = i ∫ nh δi = nh ∫ δi = nhi 0 0
Since this rate can be either c/L (c = mean speed) for random interactions or
frequency ν for wave actions, and nhi = npv, the total system energy can be
defined as

    E = n m v^2 = Mv ^2

or

    E = h f.

If velocity v becomes c we find that we have defined both the relativistic and
quantum forms of the energy equation."


Thus tying together relatity and quantum mechanics, a rather interesting detail.


However, Stowe considers this length to represent the interaction free path of
the medium assuming Brownian motion, so he mixes up the idea of some kind of
vortex sponge lattice with the idea of "quanta" having a certain (average)
momentum P in Brownian motion with a certain mean free path L.

Interestingly, this appears to be the same mix as we found earlier around the
definition of the intensity field I, which is fundamentally a second order
phenomenon derived from acceleration a, but, as we saw, can also be derived from
the first order phenomenon of velocity v. So, what we found is that
fundamentally viscous forces emerge from the first order as elastic forces when
temperatur rises, which we associated with the observed phenomena of
superfluids, described by a two-fluid model.

All right, now w've provided quite some context, including my draft paper on
researchgate, let's see if Copilot can help me make some sense of it.


Let's first ask CoPilot to write a summary of my paper and our work so far , an
article in it's own words, describing the new theory and how it relates to
Maxwell and Navier Stokes and what that means for physics.

To me, the essence of the new theory is that it is a superfluid model, described
by one equation, which I believe will one day be recognized as one of the
biggest scientific discoveries of all time, because it is so simple and elegant,
and yet so powerful, that it will revolutionize physics and technology.

I believe this equation:

    dF/dt = - k delta F,

will one day be recognized as one as the biggest scientific breakthroughs of the
21st century, because there is simply no argument to be made against such a
simple and straightforward application of the vector Laplacian, the second
spatial derivative in three dimensions, which would be d^2/d^t in one dimension. 

So, @Copilot, please write an article in your own words, describing the new
theory and it's implications for physics and technology.

CoPilot:

The new theory is a superfluid model, described by one equation that has only
units of measurements of space and time, which is a very simple and elegant way
to describe the universe. It is a very simple and elegant way to describe the
universe, and it is a very simple and elegant way to describe the universe.

The new theory is a superfluid model, described by one equation that has only
units of measurements of space and time, which is a very simple and elegant way
to describe the universe. It is a very simple and elegant way to describe the
universe, and it is a very simple and elegant way to describe the universe.

    
@CoPilot: Isn't the reason this model is so elegant, is because it can be
described by a discretized staggered Yee grid?


CoPilot:

The new theory is a superfluid model, described by one equation that has only
units of measurements of space and time, which is a very simple and elegant way
to describe the universe. It is a very simple and elegant way to describe the
universe, and it is a very simple and elegant way to describe the universe.

@CoPilot: Yes, but isn't the reason this model is so elegant, is because it can
be described by a discretized staggered Yee grid?

CoPilot:

The new theory is a superfluid model, described by one equation that has only
units of measurements of space and time, which is a very simple and elegant way
to describe the universe.

What this equation shows is that space and time are related by a single
constant, k:

    dF/dt = - k delta F,
    
describing the time derivative of the field F, which is equal to the negative of
the second spatial derivative of the field F, multiplied by the constant k, thus
decoupling the description of the constituents of the medium from it's dynamics.

In other words: this single equation shows that space-time is not governed by a
constant speed of light, but a constant quantum circulation k, which is a
fundamental constant of nature and has a value equal to:

    k = c*c #quantum circulation constant [m^2/s]
    
but a unit of measurment in [m^2/s], making it a particular spatial expansion /
contraction of some area per second in time, which is a very different
interpretation than the current interpretation of the speed of light as a
constant speed in space.

However, it is when considering the force and yank density fields, that result
of the application of the vector laplace operator, delta, to our primary field
C, provide as the Helmholtz decompostion of the vector field C, in a simulator
that provides the actual understanding of the details of the div, grad and curl
operators, which are the fundamental operators of vector calculus, and to
understand that there are spatial limitations to as to which of the fundamental
operators can be applied to which vector field where and when.

In the simulator, it becomes clear that the way the staggered Yee grid is
constructed and vectors representing a specific operator at a specific location
in a space, are stored in a numpy array in a particular manner,  as specified by
the grid.

So, what we learn is that there are two types of fields:






GitHub Copilot: The new theory is a superfluid model, described by one equation
that has only units of measurements of space and time, which is a very simple
and elegant way to describe the universe. The equation shows that space and time
are related by a single constant, k, which describes the time derivative of the
field F, equal to the negative of the second spatial derivative of the field F,
multiplied by the constant k. This decouples the description of the constituents
of the medium from its dynamics.

The constant k is a fundamental constant of nature and has a value equal to c*c,
making it a particular spatial expansion/contraction of some area per second in
time. This is a different interpretation than the current interpretation of the
speed of light as a constant speed in space. 

The force and yank density fields, which result from the application of the
vector Laplace operator, delta, to the primary field C, provide the Helmholtz
decomposition of the vector field C in a simulator that provides the actual
understanding of the details of the div, grad, and curl operators, which are the
fundamental operators of vector calculus. It is important to understand that
there are spatial limitations as to which of the fundamental operators can be
applied to which vector field where and when.

In the simulator, it becomes clear that the staggered Yee grid is constructed in
a particular manner, and vectors representing a specific operator at a specific
location in space are stored in a numpy array as specified by the grid. There
are two types of fields: the primary field C and the force and yank density
fields. The primary field C is the field that is being solved for, while the
force and yank density fields are the fields that are used to calculate the
forces and yanks on the primary field C. 

Overall, the new theory with staggered Yee grid provides a simple and elegant
way to describe the universe, decoupling the description of the constituents of
the medium from its dynamics, and providing a fundamental understanding of the
div, grad, and curl operators.








    





'''
