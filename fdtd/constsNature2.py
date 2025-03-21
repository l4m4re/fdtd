c = 299792458.0  # [m/s] - speed of light
d = 1.0  # define value of d
delta = 1.0  # define value of delta
dt = - k * d * delta  # calculate value of dt using the formula
c_l   = c * (pi/2)           # [m/s] - speed of longitudinal (sound) waves

#mu    = 4*pi*1.00000000082e-7# mu_0 - 1/viscosity - 2019 redefinition

h     = 6.62607015e-34       # Planck's constant    [J-s], [kg/m2-s]      
k     = c*c                  # dynamic viscosity    [m2/s]
rho   = eta/k                # 8.854187817620389e-12
q     = 1.602176634e-19      # elemental charge
m     = h/k                  # elemental mass   1




'''

Intermezzo: my current “to do” list, some cut&pastes from discussions on
researchgate:

Where we are now is that we can describe both the quantum level as well as the
superfluid level (quantum phenomena at a macroscopic scale) with the same
equations, only different parameters like mass density and quantum circulation
constant, whereby we find that fundamentally there are only viscous forces.

What we see with superfluids is that when temperature rises (power density
increases), elastic behavior emerges, which is currently described with a
two-fluid model. It seems that this effect can be attributed to the formation of
some kind of vortex sponge which gives rise to elastic behavior. And it also
seems we can describe the effects this creates within a continuum by the
definition of fields that are derived from the intensity field [I] rather than
our primary field [C].

The fields that can be defined as the second spatial derivatives of [I] have a
unit of measure- ment describing the time derivative of force density, which
would be yank density. What it appears to come down to is that within current
physics Force and Yank have been considered as one and the same thing, resulting
in 3D equations that break the fundamental symmetry demanded by the vector
Laplace operator.

So, it seems that there are actually two versions of Newton's law, which have
currently been taken together into one:

1) F = m a,

2) Y = m j,

and the challenge thus comes down to figuring out which one of the two applies
where. When we put these quantities in a table:

Action: Energy: Power: kg-mˆ2/s. kg-mˆ2/sˆ2 kg-mˆ2/sˆ3 Momentum mv: kg-m/s.
Force ma: kg-m/sˆ2 Yank mj: kg-m/sˆ3 Momentum density: kg/mˆ2-s Force density:
kg/mˆ2-sˆ2 Yank density: kg/mˆ2-sˆ3

it also seems that additional fields can be defined to describe action density
and its spatial de- rivative, momentum density.

----
"Sorry, but the root of your special problem is not vector analysis. It is your naïve assumption
that you are free in selecting parts of the Navier Stokes equation to handle special problems."

Well, I must confess I was a bit too fast by assuming that because I started
from the vector Laplace operator and all seemed to fit seamlessly, I had solved
the puzzle and that the loss of a few independent stress components was nothing
to worry about. So, guilty as charged in that respect.

However, it was not an exercise in selecting parts of Navier Stokes equations
that met my needs, it was an attempt to derive equivalents of Navier Stokes from
vector potential theory and to align these with equivalents of Maxwell and to
derive both from one and the same equation, which turned out to represent
Newton's third in 3D.

Since I was familiar with the scalar and vector potentials used in Maxwell and I
found that the terms in the vector Laplace operator can be written out and
define fields that establish a Helmholtz decomposition, I became convinced that
this is the way it should be done. When I searched for usage of a vector
potential in fluid dynamics, I found this paper and not much more:

https://pdfs.semanticscholar.org/9344/48b028a3a51a7567c2b441b5ca3e49ebb85c.pdf

As I wrote in my paper, I attempted to define a primary vector field for the
Laplace operator to work on for these, since that should exist according to the
Helmholtz decomposition. It seemed that all I needed to do was negate the
definition for the scalar potential, but then the unit of measurement for the
primary field turned out to be in [mˆ3/s], denoting a volumetric flow velocity,
which results in the null vector when taking the limit of the volume to zero. So
that didn't work out very well.

After a lot of puzzling, I found a solution that involved viscosity, whereby I
found that the kinematic viscosity nu yielded a value equal to light speed
squared for the aether, but a mis- match in units of measurement by a per
second, pointing to problems here and there with time derivatives. When I
realized that this constant nu can also be seen as the quantum circulation
constant, I became convinced I'm on the right track and that the thus far
mysterious properties of superfluids (quantum phenomena on a macroscopic scale)
offer the key to unlocking the mysteries of quantum mechanics.


"The Navier Stokes equations have been derived from momentum conservation. For
an in- compressible fluid we get two partial differential equations for density
and pressure. For a compressible fluid the energy balance must be considered,
which brings temperature and heat capacity into the game."


It is rather interesting that the fields I derived from my primary field [C] do
seem to describe an incompressible fluid (viscous behavior), while we seem to
have lost compressibility and that that should bring temperature and heat
capacity into the game. My working hypothesis is that temperature is a measure
of power density and has a unit of measurement in Watts per cubic meter [W/mˆ3],
but that may not be correct since Stowe (see below) found a unit in [kg-m/sˆ3].

I found a paper regarding superfluids, wherein it is stated that "second sound"
waves exist in a superfluid, which incorporates the propagation of fluctuations
in temperature:

https://sites.fas.harvard.edu/~phys191r/References/e1/donnelly2009.pdf

According to Donnely, this phenomena “has turned out to be an incredibly
valuable tool in the study of quantum turbulence”.

Thus, we have quite some hints suggesting that elastic behavior, or
compressibility, indeed has to do with the (spatial derivatives of) the
intensity field [I] I thus far payed little attention to. I've updated my
overview table and also included another primary field [Q] of which the second
spatial derivatives yield momentum density or mass flux, which I see as another
step forward.

What I think is an important detail is that the vector Laplace operator is the
3D generalization of the second spatial derivative, which would be dˆ2/dxˆ2 in
1D. This means that the 3D com- plexity of the vector equations we can define
with these three vector fields [Q], [C] and [I], such as the vector wave
equation, can be effortlessly reduced to one dimension to describe phenomena
like for instance the mechanical behavior of a long rod or a long thin tube
filled with a fluid.

"The possible approximations are “incompressibility”, “ideal gas”, or even
”perfect gas” with a constant heat capacity. Another issue are the boundary
conditions inclusive external sources and sinks, which define the geometry of
the considered problem. Finally, the initial values are important. With your
approach you stay outside of the terminology used to define Navier Stokes types
of problems."


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

The new theory is based on the vector LaPlace operator and the discovery of the
quantum circulation constant k, with a value equal to c*c but a unit of
measurement in [m/s^2]. This is a new constant of nature and represents the
quantum circulation constant, kinematic viscocity or diffusivity of the medium:

    k = c*c = 8.987551787368176e+16 [m/s^2].

The fundamental equation of the new theory is the definition of the time
derivative of any given vector field within the superfluid medium:

    dF/dt = -k * delta(F),
    
where F is any given vector field, delta is the vector LaPlace operator and k is
the quantum circulation constant. 

Note that this definition only has units of measurement in meters and seconds
and thus it represents a fundamental equation of nature describing the
fundamental dynamics of space-time itself and is decoupled from the substance of
the medium, the latter being represented by mass density rho in the model.

We believe that one day this equation will be recognized as one of the biggest
scientific breakthroughs of the 21st century, because there simply is no
argument to be made against such a simple and straightforward application of the
vector Laplacian using only one single Universal constant: k.

The vector LaPlace operator is commonly defined as follows:

    delta(F) = grad(div(F)) - curl(curl(F)),
    
where grad is the gradient operator, div is the divergence operator and curl is
the curl operator. However, during the implementation of this simulator, it was
found that because of the design of the staggered Yee grid, the curl operator
actually has to be split into two separate operators, one for the left handed
curl and one for the right handed curl, so we redefine the vector LaPlace as
follows:
    
    delta(F) = grad(div(F)) - curl_l(curl_r(F)),
        
where curl_l is the left handed curl operator and curl_r is the right handed
curl operator.  

This definitiion can be negated:

    - delta(F) = - grad(div(F)) + curl_l(curl_r(F)),

And thus we can compute the time derivative of any given vector field F along:

    dF/dt = -k * delta(F) = -k * grad(div(F)) + k * curl_l(curl_r(F)),

or:

    dF/dt = -k * delta(F) = k ( -grad(div(F)) + curl_l(curl_r(F)) ),

which can be rewritten to an update equation for F:

    F += -k * delta(F) * dt = k ( -grad(div(F)) + curl_l(curl_r(F)) ) * dt,
    
with dt the time step of the simulation.

Using the definition of the vector LaPlace operator, we can define several
fields within one consistent framework. The most simple field to define is the
flow velocity vector field v, so we can compute it's time derivate, the
acceleration vector field a, as follows:

    a = dv/dt = -k * delta(v) = k ( -grad(div(v)) + curl_l(curl_r(v)) ).

This is known as the velocity diffusion equation and is the fundamental equation
of fluid dynamics.

And, by applying Newton's F=ma, we can define the 3D force density vector field
F as follows:

    F = rho * a = -k * rho * delta(v) = k * rho * (-grad(div(v)) +
    curl_l(curl_r(v)) ).
    
One of the big benefits of the the fundamental equation of the new theory,
defining the time derivative of any given vector field F, is that we can now
apply this equation to the acceleration vector field a to compute the time
derivative of a, the jerk vector field j:

    j = da/dt = -k * delta(a) = k ( -grad(div(a)) + curl_l(curl_r(a)) ).
    
And, by applying Y=mj, with Yank the time derivative of Force, we van the 3D
yank density vector field Y as follows:  
    
    Y = rho * j = -k * rho * delta(a) = k * rho * (-grad(div(a)) +
    curl_l(curl_r(a)) ).
            
Thus, we come to a higher order model incorporating not only acceleration and
force density vector fields, but also their time derivatives: the jerk and yank
density vector fields.

Note that we have used lower case for the v, a and j vector, contrary to the
convention of this library to denote vector fields with capital letters, because
these are commonly used symbols for velocity, acceleration and jerk.

This means that we have obtained a model that is superior to both Maxwell's
equations as well as Navier-Stokes equations, because it is a higher order
model. This has far reaching consequences for our understanding of physics and
the medium, because it is known that the intermediate fields we obtain by
writing out the definition of the vector LaPlace operator, are 100%
mathematically consistent and free of any paradoxes or inconsistencies. Also,
any wave equations that can be derived from this model will be invariant to the
Galilean coordinate transformation, bringing the Lorentz tranform on which
Einstein's theory of relativity is based, into question.

Furthermore, the introduction of the jerk and yank vector fields into the model
enables us to model the dynamics of the medium in a much more accurate way, as
from control systems theory it is known that not only the acceleration of an
actuator needs to be controlled, but also the jerk. Otherwise, all kinds of
vibrations and oscillations will occur in the system. 

In other words: since it is known from control systems that the jerk of the
system under control needs to be controlled as well, it is only logical that
this should be the case for the medium as well. And, since neither Maxwell's
equations nor Navier-Stokes are higher order equations, it is clear that the old
models are incomplete and therefore not capable of describing either the medium
or a fluid accurately. 

Now let's consider the negated definition of the vector LaPlace operator for a
field C:

   - delta(C) = - grad(div(C)) + curl_l(curl_r(C)),

with C = rho * k * v, resulting in a unit of measurement of:

    [C] = [rho] * [k] * [v] = [kg/m^3] * [m^2/s] * [m/s] = [kg/s^2].

Thus, C represents a thusfar unknown quantity that has a unit of measurement in
[kg/s^2], [N/m], [J/m^2] or [Pa-m]. And since in our model the Ampere
representing electric current also has a unit of measurement in [kg/s^2], we
have chosen to use the symbol C. However, since electric current in our model is
considered to represent turbulence within the medium, it can be speculated that
C also represents some quantity related to turbulence. For now, we shall refer
to C as our primary dynamic field.

From this, we can write out intermediate fields similar to Maxwell's equations,
with the big difference that these intermediate fields are derived from a higher
order equation, unlike Maxwell's:

    phi = div(C) 
    
    A   = curl_r(C)
    
    E  = -grad(phi) 
    
    B  =  curl_l(A)

And because of vector identities we can also write:    

    curl(E) = 0 
    
    div(B)  = 0,
    
since the divergence of a curl of the gradient of any twice-differentiable
scalar field p is always the zero vector (curl(grad(p))=0), and the divergence
of the curl of any vector field A is always zero as well (div(curl(A))=0).

Thus, we have shown that the fields derived from application of the vector
Laplacian yield a Helmholtz decomposition of the field C, decomposing it into a
divergence free component B and a curl free component E.

From wikipedia:

"In physics and mathematics, in the area of vector calculus, Helmholtz's
theorem, also known as the fundamental theorem of vector calculus, states that
any sufficiently smooth, rapidly decaying vector field in three dimensions can
be resolved into the sum of an irrotational (curl-free) vector field and a
solenoidal (divergence-free) vector field; this is known as the Helmholtz
decomposition or Helmholtz representation. It is named after Hermann von
Helmholtz." 

In other words: the fields we have defined are consistent with the fundamental
theorem of vector calculus.

We can compare this with Maxwell's equations:

    E = -grad(phi) - dA/dt, 
    
    B =  curl(A),
    
    curl(E) = -dB/dt, 
    
    div(B)  =  0,
    
    div(E)  = rho / eps,
    
    curl(B) = mu * J - 1/c^2 * dE/dt,
    
Where obviously the curl of the electric field is not equal zero and therefore
Maxwell's equations are not consistent with the fundamental theorem of vector
calculus: the Helmholtz theorem or Helmholtz decomposition as fundamentally
described by the vector Laplace operator.

In other words: we have shown that Maxwell's equations are mathematically
inconsistent.   

Note that when taking the curl on both sides of this equation:

    E = -grad(phi) - dA/dt, 

we obtain the Maxwell-Faraday equation, representing Faraday's law of induction:

    curl(E) = -dB/dt.

Faraday's law of induction is a basic law of electromagnetism predicting how a
magnetic field will interact with an electric circuit to produce an
electromotive force (EMF), which is thus a law that applies at the macroscopic
level. Because of the inconsistencies in Maxwell's equations with respect to the
fundamental theorem of vector calculus, it is now clear that this entanglement
of the fundamental E and B fields in Maxwell's equations is responsible for a
lot of problems in current main stream physics, for starters the non-invariance
of Maxwell's equations to the Galilean transform.

This is described by C.K. Thornhill in his paper "Real or Imaginary Space-Time?
Reality or Relativity?":

    https://etherphysics.net/CKT4.pdf   

    "ABSTRACT 

    The real space-time of Newtonian mechanics and the ether concept is
    contrasted with the imaginary space-time of the non-ether concept and
    relativity. In real space-time (x, y, z, ct) characteristic theory shows
    that Maxwell's equations and sound waves in any uniform fluid at rest have
    identical wave surfaces. 

    Moreover, without charge or current, Maxwell's equations reduce to the same
    standard wave equation which governs such sound waves. This is not a general
    and invariant equation but it becomes so by Galilean transformation to any
    other reference-frame. So also do Maxwell's equations which are, likewise,
    not general but unique to one reference-frame. 

    The mistake of believing that Maxwell's equations were invariant led to the
    Lorentz transformation and to relativity; and to the misinterpretation of
    the differential equation for the wave cone through any point as the
    quadratic differential form of a Riemannian metric in imaginary space-time
    (x, y, z, ict). 

    Mathematics is then required to tolerate the same equation being transformed
    in different ways for different applications. Otherwise, relativity is
    untenable and recourse must then be made to real space-time, normal Galilean
    transformation and an ether with Maxwellian statistics and Planck's energy
    distribution."

Yet another problem with Maxwell's equations is that the potential fields phi
and A are not unique, because they are defined as sort of an inverse of E and B,
but are not derived from some primary field like the C field in the new theory.

It is this non-uniqueness of the scalar and vector potentials in Maxwell's
equations that leaves room for "gauge freedom" on top of which a whole set of
"gauge theories" have been built:

    https://en.wikipedia.org/wiki/Gauge_theory
    
    "The earliest field theory having a gauge symmetry was Maxwell's
    formulation, in 1864-65, of electrodynamics ("A Dynamical Theory of the
    Electromagnetic Field") which stated that any vector field whose curl
    vanishes — and can therefore normally be written as a gradient of a function
    — could be added to the vector potential without affecting the magnetic
    field. 
    
    [...]
    
    The importance of gauge theories in physics is exemplified in the tremendous
    success of the mathematical formalism in providing a unified framework to
    describe the quantum field theories of electromagnetism, the weak force and
    the strong force. This theory, known as the Standard Model, accurately
    describes experimental predictions regarding three of the four fundamental
    forces of nature, and is a gauge theory with the gauge group SU(3) x SU(2) x
    U(1). Modern theories like string theory, as well as general relativity,
    are, in one way or another, gauge theories." 
    
Note that because of the vector identity curl(grad(p))=0, addition of a vector
field whose curl vanishes to the vector potential has no effect whatsoever on
the resulting field, which begs the question of why the addition of a vector
field that has no effect on the resulting field, is supposed to describe
something that has any relevance to the real world at all.

Or, to put it another way: how exactly is the Standard Model, which is a gauge
theory and therefore uses fields that have no effect on the resulting physical
fields, is supposed to have any relevance to the real world at all?

Thus, we can conclude that Maxwell's rather strange entanglement of the electric
and magnetic fields has led to a lot of problems in physics, to such an extent
that it has been accepted that theories that are not only mathematically
inconsistent but also use fields that have no effect on the resulting physical
fields, are completely acceptable and even the norm in physics.


While it is useful to write out intermediate fields similar to Maxwell's
equations for our new theory in order to compare them with Maxwell's equations,
it is also confusing to do so, because even though our intermediate fields are
similar to Maxwell's, they are not the same and should therefore be described
with different symbols.

So let's consider the negated definition of the vector LaPlace operator for a
field C again:

   - delta(C) = - grad(div(C)) + curl_l(curl_r(C)),

with C = rho * k * v, resulting in a unit of measurement of:

    [C] = [rho] * [k] * [v] = [kg/m^3] * [m^2/s] * [m/s] = [kg/s^2].

From this, we can write out intermediate fields similar to Maxwell's equations,
with completely different meanings and units of measurements:

    p = div(C) 
    
    A = curl_r(C)
    
    L = -grad(p) 
    
    R = curl_l(A)

And because of vector identities we can also write:    

    curl(L) = 0 
    
    div(B)  = 0,

And since the unit of measurement of C is [kg/s^2], we can work out the units of
measurement of the intermediate fields given that the grad, div and curl
operators all have a unit of measurement in per meter [1/m]:

    [p] = [kg/s^2] * [1/m] = [kg/m-s^2] or [Pa] = [N/m^2] = [J/m^3]  

    [A] = [kg/s^2] * [1/m] = [kg/m-s^2] or [Pa] = [N/m^2] = [J/m^3]
    
    [L] = [kg/m-s^2] * [1/m] = [kg/m^2-s^2] or [Pa/m] = [N/m^3] = [J/m^4]
    
    [R] = [kg/m-s^2] * [1/m] = [kg/m^2-s^2] or [Pa/m] = [N/m^3] = [J/m^4]

Thus we have found that our scalar and vector potentials have a unit in [Pa]
representing pressure, while our resulting [L] and [R] fields have a unit in
[N/m^3] representing force density. And, since [Pa] = [J/m^3], we can also
consider the scalar and vector potentials to represent the energy density of the
medium.

Note that the intermediate fields we have defined are consistent with the
fundamental theorem of vector calculus, unlike Maxwell's equations.

Now let's consider the negated definition of the vector LaPlace operator for a
field I:

   - delta(I) = - grad(div(I)) + curl_l(curl_r(I)),

with I = rho * k * a, resulting in a unit of measurement of:

    [C] = [rho] * [k] * [a] = [kg/m^3] * [m^2/s] * [m/s^2] = [kg/s^3].

Thus, I represents the Intensity or radiosity field in [kg/s^3], [N/m-s],
[J/m^2-s], [Pa-m/s] or [W/m^2].

From this, we can write out intermediate fields:

    t   = div(I) 
    
    W   = curl_r(I)
    
    Y_l = -grad(t) 
    
    Y_a = curl_l(W)

And because of vector identities we can also write:    

    curl(Y_l) = 0 
    
    div(Y_a)  = 0, 

And since the unit of measurement of I is [kg/s^3], we can work out the units of
measurement of the intermediate fields:

    [t] = [kg/s^3] * [1/m] = [kg/m-s^3] or [Pa/s] = [W/m^3] = [N/m^2-s] 
    
    [W] = [kg/s^3] * [1/m] = [kg/m-s^3] or [Pa/s] = [W/m^3] = [N/m^2-s] 
        
    [Y_l] = [kg/m-s^3] * [1/m] = [kg/m^2-s^3] or [Pa/m-s] = [N/m^3-s] 
    
    [Y_a] = [kg/m-s^3] * [1/m] = [kg/m^2-s^3] or [Pa/m-s] = [N/m^3-s] 
    

Thus we have found that our second order scalar and vector potentials have a
unit in [W/m^3] representing power density, while our resulting [Y_l] and [Y_a]
fields have a unit in [N/m^3-s] representing yank density, the time derivative
of force density. 

It is not far fetched to assume that the scalar power density field t in [W/m^3]
actually represents temperature, perhaps apart from a constant factor, because
it is a scalar field within a superfluid medium model and the unit of power
density feels like a natural match with temperature, so that temperature
actually is a measurement of power density within the medium. If this is
correct, then there should be no doubt that this model will lead to a much
deeper understanding of temperature and heat related phsyics in general,
especially because we also have the angular counterpart of the scalar potential,
the vector potential field W, which also has a unit of measurement of power
density.

It feels right to associate current density J with the angular power density
field W. However, the unit of measurment of J would be [A/m^3], which would in
our theory be equal to [kg/s-m^3]. 

eta = 1/mu = 1/4*pi*1e-7 = 8.987551787368176e+16 [kg/m-s2] rho = eta/k =
8.987551787368176e+16 / 8.987551787368176e+16 = 1 [kg/m3]


W in [kg/m-s^3] or [Pa/s] = [W/m^3] = [N/m^2-s] J in [kg/s-m^3] or [A]

 [kg/s-m^3]  = X [kg/m-s^3]
 
 [m^2]  = X [/s^2]
 
 x = [m^2-s^2]
 
 
 
So far, we have made a number of references to superfluids, without explaining
what a superfluid actually is. Generally, a superfluid is described as a fluid
that has zero viscosity and therefore flows without friction. However, this is
incorrect. 

An interesting detail is that the intensity field I can also be defined as:

    I = - κ v,

with κ the modulus or elasticity in [Pa] or [kg/m-s2], which has a unit of
measurement that differs by a per second [/s] from the unit of measurement for
viscosity eta in [Pa-s] or [kg/m-s].

This reflects the difference between elastic forces and viscous (shear) forces,
namely that the elastic force is proportional to the amount of deformation,
while the viscous one is proportional to the rate of deformation. 

So, it appears we can conclude that in physical reality there are no actual
static (elastic) forces (at the quantum level) and that deep down there are only
dynamic forces and interactions which are governed by the velocity diffusion
equation (a = dv/dt = -k * delta(v)), whereby what we think to observe as static
elastic forces are in reality the time derivatives of forces that are
fundamentally of a viscous nature.

This is a very interesting conclusion, because it relates nicely to the rather
strange observed behavior of superfluids, which are currently described by a
two-fluid model, whereby the superfluid component is described as a frictionless
fluid, while the normal fluid component is described as a viscous fluid. 

Since within the new model, we have two scalar potentials and two vector
potentials, along with derived force and yank density fields, and we find that
the potentials and yank density fields derived from the intensity field can also
be considered as fields describing elastic forces, it becomes comprehensible
that the current two-fluid model can perhaps be matched to our new model.

Thus, the superfluid component would be described by our primary dynamic field C
and it's derived fields, while the normal fluid component is described by the
intensity field I and its derived fields, of which the scalar potential field p
is very interestng, since we assume that one to represent (a measurement of)
temperature.

This way, we would come to a model whereby elastic behavior (I field and
derived) emerges as function of scalar potential t, which we assume to represent
temperature. Thus, this could explain how normal fluid behavior emerges from
superfluid behavior as temperature rises, which would be a very interesting
result.

One of the most interesting phenomena observed with superfluid is that quantized
vortices arise when a superfluid is rotated, whereby the quantization gamma is
given by:

    gamma = h / m,
    
with h Planck's constant and m the mass of a superfluid particle, such as
superfluid helium.

When we equate this circulatio constant gamma with the quantum circulation
constant k, we can compute the mass of an elemental medium particle along:

    m = h / k = 6.62607015e-34 / 8.987551787368176e+16 = 7.374e-51 [kg],
    
What is very interesting about this result is that we can compute the Compton
wavelength and frequency of this particle along:

    lambda = h / (m * c) = 6.62607015e-34 / (7.374e-51 * 299792458) = 2299792458
    [m].
    
So, for a particle with a mass of 7.374e-51 [kg], we find a Compton wavelength
of 2299792458 [m], which is exactly equal to the speed of light c, but with a
unit of measurement in [m] instead of [m/s].     
 
It should be no surprise that when we compute the Compton frequency of this
particle, we find:

    nu = c / lambda = 299792458 / 2299792458 = 1 [Hz].

This puts serious question marks on the validity of the mass-energy equivalence
principle, since, according to Wikipedia, the Compton wavelength is defined as
the wavelength of a photon whose energy is the same as the rest mass of a
particle.

Are we really to believe that the tiniest of the tiniest particle with a mass of
7.374e-51 [kg], about 20 orders of magnitude lighter than the electron, has a
Compton wavelength of no less than 2299792458 [m] and associated frequency of 1
[Hz]?

And why does this particular mass (7.374e-51 [kg]) has a Compton wavelength that
is exactly equal to the speed of light c? 

It is these questions, along with the mathematical inconsistensies we found in
Maxwell's equations, that reveal that there is something seriously wrong with
the current mainstream models, which not only involves Maxwell's equations and
everything built on top of it, but also the fluid mechanics Navier-Stokes
equations.

Now let's introduce a conversation with CoPilot:
------------------------------------------------

l4m4re: So, if I understand this correctly, in a fluid dynamics Yee grid you
have a scalar field on the edges of the grid  cells and a vector field on the
faces of the grid, while in electrodynamics you have vector fields both on the
edges and on the faces of the grid cells?


GitHub Copilot: Yes, that's correct. In a fluid dynamics Yee grid, the scalar
field is located on the edges of the grid cells, while the vector field is
located on the faces of the grid cells. This is because the scalar field
represents a quantity that is defined at a point in space, while the vector
field represents a quantity that has both magnitude and direction and is defined
over an area.

In electrodynamics, the electric field is located on the edges of the grid
cells, while the magnetic field is located on the faces of the grid cells. This
is because the electric field is a vector field that is defined at a point in
space, while the magnetic field is a vector field that is defined over an area.

Overall, the Yee grid is designed to accurately model the behavior of
electromagnetic waves and fluid flow by discretizing the fields on a staggered
grid with interlaced Yee-coordinates. By using different field representations
and discretization schemes, the Yee grid can be used to simulate a wide range of
physical phenomena in both electrodynamics and fluid dynamics.

-:-

This is a very important explanation of the fundamental problem with Maxwell's
equations, because the discretization reveals that there are two fundamentally
different types of vectoer fields, which are:

1) vector fields that are defined at a point in space, such as the electric 
field in electrodynamics and the scalar field in fluid dynamics as well as in
our new superfluid ether theory.

2) vector fields that are defined over an area, such as the magnetic field in
electrodynamics and the vector field in fluid dynamics as well as in our new
superfluid ether theory.

This is a very important distinction, because it reveals that the electric and
magnetic fields are fundamentally different types of vector fields. However,
this goes much further than just Maxwell's or Navier-Stokes equations, because
the distinction ultimately originates at the mathematical definition of the
three fundamental vector operators: grad, div and curl.

These three operators are defined using two methods, namely as a differential
operator and as a finite difference operator. Both the differential and finite
difference operators are defined relative to a certain volume, whereby for the
finite difference operators the volume is discretized into a grid of cells and
therefore has a finite size, while for the differential operators the limit of
the volume is taken to zero and therefore has an infinitesimally small size.

This means that for the differential operators the resulting vector fields can
all be considered to be defined at a point in space, while for the finite
difference operators the resulting vector fields are considered to be located at
specific locations on the grid, such as the faces of the grid cells for the
angular vector fields (curl operator) and the center point of the grid cell for
the linear fields consisting of a scalar (div operator) and a vector component
(grad operator).

In other words: it is the operator used that determines the location of the
resulting vector or scalar field relative to the grid cell and therefore the
type of field, either defined at a point in space (grad, div) or defined over an
area (curl).

In a staggered Yee grid, the fields are located at either integer or
half-integer grid points, whereby the integer grid points represent vector
fields defined at a point in space, while the half-integer grid points represent
vector fields defined over an area, the latter being one of the faces of the
grid cell. 

In other words: in the discretization of a staggered Yee grid each scalar and
vector field is located at a specific location relative to the cube geometry of
the grid cell. Thus, half-integer grid points represent vector fields defined
over an area and therefore a number of cube faces and are therefore located at a
grid face, while integer grid points represent vector fields defined relative to
the center of the cube, which is located half a grid cell away from the faces of
the cube hence the location on a "half-integer grid point" in the arrays used in
simulation software to represent the finite difference fields.

So, let us consider the definition of the three operators relative to the
geometry over which they are defined and compare the differential and finite
difference methods. In principle, the differential operators can be considered
as also being defined relative to a grid, even though the limit of the grid
spacing is taken to zero. This means that the limitations of the finite
difference operators also apply to the differential operators and a clear
distinction needs to be made between vector fields defined at a point in space
and vector fields defined over an area.

The divergence operator is defined as the surface integral of the flux of a
vector field over a closed surface. For the differential operator, the surface
integral is taken over a surface of an arbitrary volume around some point in
space, like a sphere or a cube, of which the limit is taken to zero. The
resulting vector is then located at the center of the volume, which is a point
in space. Since in a staggered Yee grid a discretized cube is used, we can
consider the differential operators to be defined over a cube of infinitesimally
small but not zero size. 

So, if we consider the divergence operator, we can see that it is defined as a
surface integral over the area of all the faces of a cube, divided by the volume
of the cube, resulting in a unit of per meter [1/m] for the operator. Since the
resulting value is not located on one of the faces of the cube, it can be
considered to be located at the center of the cube, which is a point in space.
Since fields defined at a point in space are represented by integer grid points
in a staggered Yee grid, the resulting scalar field is located at an integer
grid point.

While the gradient operator can be defined in multiple ways, including for
vector fields resulting in a tensor, we will only consider the gradient of a
scalar field here. 

The gradient of a scalar field is defined as the spatial derivative of the
scalar field with respect to the three spatial dimensions. For the differential
operator, the derivative is taken relative to a point in space, while for the
finite difference operator the derivative is taken relative to the grid spacing.

So, the difference of the scalar field is computed between six points at a
certain distance from the point along the three coordinate axes (+ and -), which
are added together and divided by the distance, resulting in a unit of per meter
[1/m] for the operator. In the discretized implementation, the result is not
actually divided by the grid spacing, because it is an actual difference
relative to the grid. 

Since the resulting value is not located on one of the faces of the cube, it can
be considered to be located at the center of the cube, which is a point in
space. Since fields defined at a point in space are represented by integer grid
points in a staggered Yee grid, the resulting vector field is located at an
integer grid point as is the input scalar field.

For the curl operator, we found that we needed two different curl operators,
namely the left and right curl operators. The curl operator is defined as the
line integral of the tangential component of the vector field over a closed path
around some point in space, divided by the area enclosed by the path along
Stokes' theorem, resulting in a unit of per meter [1/m] for the operator. And
since this line integral can be computed clockwise or counter-clockwise, we
obtain the left and right curl operators needed in the implementation of the
staggered Yee grid.


GitHub Copilot (edited): In a staggered Yee grid implementation, the grid
spacing is constant, which means that the ratio between the grid spacing and the
area and volume of a grid cube is also constant. This allows for simplifications
in the finite difference implementation of the gradient and curl operators, as
the difference of the scalar field can be computed between six points at a
certain distance from the point along the three coordinate axes (+ and -).

Similarly, the left and right curl operators can be computed by taking the
difference of the vector field components at neigboring faces along the three
coordinate axes (+ and -). Since the grid spacing is constant, the resulting
value is simply the difference between the vector field components at the two
neighboring grid points located at half-integer grid points since these fields
are defined over an area.


Therefore, in the implementation of the gradient and curl operators in a
staggered Yee grid, the grid spacing can be assumed to be 1.0, which simplifies
the computation of the operators using finite differences. This eliminates the
need for integrals and allows for a consistent and accurate computation using
finite differences.

-:-

So, for fields defined over an area, half-integer grid points are used, which
are not actually points, but represent te surface area of the face of the grid
cell. And for fields defined at a point in space, integer grid points are used,
which are actual points. This is the essence of the staggered Yee grid, and
therefore point fields and area fields are fundamentally different types of
fields, which cannot be mixed.

So, let's consider the Maxwell-Faraday equation:

    curl(E) = -dB/dt.
    
This equation is defined relative to a point in space, because the curl operator
is defined relative to a point in space. This means that the electric field E is
defined relative to a point in space, which means that it is a point field and
therefore should be represented by an integer grid point.

However, the magnetic field B is defined relative to an area, because the
derivative of the magnetic field is taken relative to the grid spacing, which
means that it is a surface field and therefore should be represented by a
half-integer grid point.

So, the Maxwell-Faraday equation is fundamentally flawed, because it mixes a
point field with a surface field, which is not allowed in a staggered Yee grid.

So, let's consider the Maxwell-Ampere equation: 

    curl(B) = mu * J - 1/c^2 * dE/dt,
    
This equation is defined relative to a point in space, because the curl operator
is defined relative to a point in space. This means that the magnetic field B is
defined relative to a point in space, which means that it is a point field and
therefore should be represented by an integer grid point.

However, the current density J is defined relative to an area, because the
derivative of the current density is taken relative to the grid spacing, which
means that it is a surface field and therefore should be represented by a
half-integer grid point.

So, the Maxwell-Ampere equation is fundamentally flawed, because it mixes a
point field with a surface field, which is not allowed in a staggered Yee grid.

So, let's consider the Maxwell-Gauss equation:

    div(E)  = rho / eps,
    
This equation is defined relative to a point in space, because the divergence
operator is defined relative to a point in space. This means that the electric
field E is defined relative to a point in space, which means that it is a point
field and therefore should be represented by an integer grid point.

However, the charge density rho is defined relative to an area, because the
derivative of the charge density is taken relative to the grid spacing, which
means that it is a surface field and therefore should be represented by a
half-integer grid point.

So, the Maxwell-Gauss equation is fundamentally flawed, because it mixes a point
field with a surface field, which is not allowed in a staggered Yee grid.




GitHub Copilot: Thank you for the clarification. You are correct that the
magnetic field `B` is defined relative to an area, and is therefore a surface
field. The derivative of the magnetic field `B` is taken relative to the grid
spacing, which means that it should be represented by a half-integer grid point.

In the Yee grid implementation, the magnetic field `H` is defined on the edges
of the grid cells, which are integer grid points, while the auxiliary magnetic
field `B` is defined on the faces of the grid cells, which are half-integer grid
points. This is because the `H` field is related to the current density, which
is a point quantity, while the `B` field is related to the magnetic flux
density, which is a surface quantity.

-:-

This, however is inconsistent. The thinking error is that current density can be
considered a point quantity, while in reality current always follows a closed
loop and is therefore a surface quantity.

So, let's consider the Maxwell-Ampere equation again:
    
        curl(B) = mu * J - 1/c^2 * dE/dt,
        
Copilot: This equation is defined relative to a point in space, because the curl
operator is defined relative to a point in space. This means that the magnetic
field B is defined relative to a point in space, which means that it is a point
field and therefore should be represented by an integer grid point.

This, however, is incorrect. Yes, in infinitesimal consideration, one can
consider the curl operator to be defined relative to a point in space, but that
is misleading, since the curl operator is defined relative to a closed path,
which is a surface quantity. So, the point is that even when taking the limit of
the path to zero, it still is a path and not a point in space.

Thus, we come to the conclusion that both the H and B fields should be defined
relative to a surface, which means that they should be represented by a
half-integer grid point.


'''


'''

﻿
Revision and integration of Maxwell’s and Navier-Stokes’ Equations and the
origin of quantization in Superfluids and Spacetime itself

Arend Lammertink

Schoolstraat 107, 7471 WV, Goor, The Netherlands, +316 5425 6426,
lamare@gmail.com

“I hope that someone will discover a more realistic way, or rather a more
tangible basis than it has been my lot to find.”

Albert Einstein. DRAFT revision 320. 

Abstract. It is well known that the Maxwell equations predict the behavior of
the electromagnetic field very well. However, they predict only one wave
equation while there are significant differences between the "near" and "far"
fields and various anomalies have been observed involving the detection of super
luminous signals in experiments with electrically short coaxial cables1,2,
microwaves3,4,5,6,7,8,9,10,11,12, optical fibers13,14,15,16 as well as other
methods17,18,19,20.  We show that the mathematical Laplace operator defines a
complete set of vector fields consisting of two potential fields and two fields
of force, which form a Helmholtz decomposition of any given vector field F. We
found that neither in Maxwell’s equations nor in fluid dynamics vector theory
this result has been recognized, which causes the potential fields to not be
uniquely defined and also makes the Navier-Stokes equations unnecessarily
complicated and introduces undesirable redundancy as well. We show that
equivalents to both the Maxwell equations as well as the Navier-Stokes equations
can be directly derived from a single diffusion equation describing Newton’s
second law in 3D. We found that the diffusion constant ν in this equation has
the same value as the speed of light squared, but has a unit of measurement in
meters squared per second thus uncovering problems with time derivatives in
current theories, showing amongst others that the mass-energy equivalence
principle is untenable. Finally, we show that the diffusion equation we found
can be divided by mass density ρ, resulting in a velocity diffusion equation
that only has units of measurement in meters and seconds, thus decoupling the
dynamics of the medium from it’s substance, mass density ρ. This reveals the
quantized nature of spacetime itself, whereby the quantum circulation constant ν
is found to govern the dynamics of physical reality, leading to the conclusion
that at the fundamental quantum level only dynamic viscous forces exist while
static elastic forces are an illusion created by problems with a number of time
derivatives in current theories.

With our equivalents for the Maxwell equations three types of wave phenomena can
be described, including super luminous longitudinal sound-like waves that can
explain the  mentioned anomalies. This paper contributes to the growing body of
work revisiting Maxwell’s equations21,22,23,24,25,26,27,28,29,30 by deriving all
of the fields from a single equation, so the result is known to be
mathematically consistent and free of singularities and uniquely defines the
potential fields thus eliminating gauge freedom. Unlike Maxwell’s equations,
which are the result of the entanglement of Faraday's circuit level law with the
more fundamental medium arguably creating most of the problems in current
theoretical physics, these revisions describe the three different
electromagnetic waves observed in practice and so enable a better mathematical
representation.

Keywords: Classical Electrodynamics, Superfluid medium, Fluid Dynamics,
Theoretical Physics, Vector Calculus.


Introduction

In 1861, James Clerk Maxwell published his paper “On Physical Lines of Force”31,
wherein he theoretically derived a set of twenty equations which accurately
described the electro-magnetic field insofar as known at that time. He modeled
the magnetic field using a molecular vortex model of Michael Faraday's "lines of
force" in conjunction with the experimental result of Weber and Kohlrausch, who
determined in 1855 that there was a quantity related to electricity and
magnetism, the ratio of the absolute electromagnetic unit of charge to the
absolute electrostatic unit of charge, and determined that it should have units
of velocity. In an experiment, which involved charging and discharging a Leyden
jar and measuring the magnetic force from the discharge current, they found a
value 3.107e8 m/s, remarkably close to the speed of light.

In 1884, Oliver Heaviside, concurrently with similar work by Josiah Willard
Gibbs and Heinrich Hertz, grouped Maxwell’s twenty equations together into a set
of only four, via vector notation. This group of four equations was known
variously as the Hertz-Heaviside equations and the Maxwell-Hertz equations but
are now universally known as Maxwell's equations. 

The Maxwell equations predict the existence of just one type of electromagnetic
wave, even though it is now known that at least two electromagnetic wave
phenomena exist, namely the “near” and the “far” fields. The “near” field has
been shown to be a non-radiating surface wave that is guidable along a
completely unshielded single conductor32 and can be applied for wide band low
loss communication systems. The Maxwell equations have not been revised to
incorporate this new knowledge.  

Given the above, the following questions should be asked: 

    • What is charge?
    • Why is it a property of certain particles?

As long as we insist that charge is an elemental quantity that is a property of
certain particles, we cannot answer these questions. Also, when the wave
particle duality principle is considered in relation to what is considered to be
the cause for electromagnetic radiation, charged particles, in Maxwell’s
equations electromagnetic radiation is essentially considered to be caused by
(quanta of) electromagnetic radiation, an obvious case of circular logic which
is not desirable. 

In the area of vector calculus, Helmholtz's theorem, also known as the
fundamental theorem of vector calculus, states that any sufficiently smooth,
rapidly decaying vector field in three dimensions can be resolved into the sum
of an irrotational (curl-free) vector field and a solenoidal (divergence-free)
vector field; this is known as the Helmholtz decomposition. A terminology often
used in physics refers to the curl-free component of a vector field as the
longitudinal component and the divergence-free component as the transverse
component. This theorem is also of great importance in electromagnetic (EM) and
microwave engineering, especially for solving the low-frequency breakdown issues
caused by the decoupling of electric and magnetic fields.33 Further, a vector
field can be uniquely specified by a prescribed divergence and curl and it can
be shown that the Helmholtz theorem holds for arbitrary vector fields, both
static and time-dependent34.

In potential theory, the study of harmonic functions, the Laplace equation is
very important, amongst other with regards to consideration of the symmetries of
the Laplace equation. The symmetries of the n-dimensional Laplace equation are
exactly the conformal symmetries of the n-dimensional Euclidean space, which has
several implications. One can systematically obtain the solutions of the Laplace
equation which arise from separation of variables such as spherical harmonic
solutions and Fourier series. By taking linear superpositions of these
solutions, one can produce large classes of harmonic functions which can be
shown to be dense in the space of all harmonic functions under suitable
topologies. 

The Laplace equation as well as the more general Poisson equation are 2nd order
differential equations, in both of which the Laplacian represents the flux
density of the gradient flow of a function. In one dimension, the Laplacian
simply is ∂²/∂x², representing the curvature of a given function f. For scalar
functions in 3D, the Laplacian is a common generalization of the second
derivative and is the differential operator defined by:


(1)

The Laplacian of a scalar function is equal to the divergence of the gradient
and the trace of the Hessian matrix. The vector Laplacian is a further
generalization in three dimensions and defines the second order spatial
derivative of any given vector field F, the “3D curvature” if you will, and is
given by the identity:


(2)

Whereas the scalar Laplacian applies to a scalar field and returns a scalar
quantity, the vector Laplacian applies to a vector field, returning a vector
quantity. When computed in orthonormal Cartesian coordinates, the returned
vector field is equal to the vector field of the scalar Laplacian applied to
each vector component. 

With this identity, a full 3D generalization of the Poisson equation can also be
defined, the vector wave equation, which has three independent solutions35, the
vector spherical harmonics:


(3)


Methods

The terms in the definition for the vector Laplacian can be negated:


(4)

and then the terms in this identity can be written out to define a vector field
for each of these terms:  


(5)

And, since the curl of the gradient of any twice-differentiable scalar field φ
is always the zero vector (∇×(∇φ)=0), and the divergence of the curl of any
vector field A is always zero as well (∇⋅(∇×A)=0), we can establish that E is
curl-free and B is divergence-free, and we can write: 


(6)

As can be seen from this, the vector Laplacian establishes a Helmholtz
decomposition of the vector field F into an irrotational or curl free component
E and a divergence free component B, along with associated potential fields φ
and A, all from a single equation c.q. operator. 

Thus we have shown that the mathematical definitions for potential fields are
hidden within the Laplace operator c.q. the fundamental theorem of vector
calculus c.q. the second order spatial derivative, which has tremendous
consequences for both the analytical analysis of the electromagnetic field as
well as fluid dynamics vector theory. The symmetry between the fields thus
defined is fundamental and has been mathematically proven to be correct, so it
is vital to maintain this fundamental symmetry in our physics equations. 

So far, we have considered the general case, which is valid for any given vector
field F. In the following, we will use the m subscript to refer to the
electromagnetic domain along Maxwell’s equations, while the f subscript is used
for the fluid dynamics domain. 

In Maxwell's equations, the curl of the electric field Em is defined by the
Maxwell-Faraday equation:


(7)

which is obvious not equal to zero for electromagnetic fields varying with time
and therefore the dynamic Maxwell equations cannot be second order spatial
derivatives of any vector field Fm as defined by the Laplacian. 

Herewith, we have shown that no vector field Fm exists for which Maxwell's
equations are the second order spatial derivative and therefore Maxwell’s
equations do not satisfy the vector Laplace equation. The end result of this is
that while the solutions of Laplace’s equation are all possible harmonic wave
functions, with Maxwell’s equations there is only one resulting wave equation
which defines a “transverse” wave, whereby the Em and Bm components are always
perpendicular with respect to one another. This is also the reason why no
separate wave equations can be derived for the “near” and “far” fields.

Furthermore, in Maxwell’s equations, the two potential fields which are used
with Helmholtz’s theorem are the electrical potential φm and the magnetic vector
potential Am, which are defined by the equations36:


(8)

where Bm is the magnetic field and Em is the electric field.

The Helmholtz theorem can also be described as follows. Let A be a solenoidal
vector field and φ a scalar field on R3 which are sufficiently smooth and which
vanish faster than 1/r2 at infinity. Then there exists a vector field F such
that:


(9)

and if additionally, the vector field F vanishes as r → ∞, then F is unique37.

Now let us consider the units of measurement involved in these fields, whereby
the three vector operators used all have a unit of measurement in per meter
[/m]. The magnetic field Bm has a unit of measurement in Tesla [T], which is
defined in SI units as [kg/s2-A]. So, for the magnetic vector potential Am we
obtain a unit of [kg-m/s2-A] and for dAm/dt we obtain a unit of [kg-m/s3-A]. The
electric field Em has a unit of measurement in volt per meter, which is defined
in SI units as [kg-m/s3-A], which matches that for dAm/dt. So, for the electric
scalar potential φm we obtain a unit of [kg-m2/s3-A].

However, neither the units of measurement for Em and Bm are the same, nor are
the units of measurements for φm and Am. This is in contradiction with
Helmholtz’s theorem, which states that a vector field Fm exists that should have
a unit of measurement equal to that of φm and Am times meters or that of Em and
Bm times meters squared.

Thus, we have shown that Maxwell’s equations are in contradiction with
Helmholtz’s theorem as well, which means that the potential fields defined by
Maxwell are mathematically inconsistent and should therefore be revised.

It can be shown38 that by using the 19th Century’s atomic vortex postulate in
combination with a superfluid model for the medium, it is possible to construct
a single simple integrated model which covers all major branches of physics
including kinetic, fluid, gravitation, relativity, electromagnetism, thermal,
and quantum theory. With this method, it can also be shown that anomalous
observations such as Pioneer’s drag and the electron’s magnetic moment can be
directly accounted for by the model. Furthermore, with this model all units of
measurements are defined in terms of just three fundamental units of
measurement: mass, length, and time. 

It should be noted that there are two distinct levels in this model, with each
playing their own role. The first consists of basic media quanta, which forms a
superfluid model for the medium itself. The second describes vortices within the
fluid, which forms a particle model on top of the medium model. The lower base
level is assumed to be an (if not ideal, nearly so) in-viscous superfluid system
obeying the defined rules of basic kinetic theory and that is the model this
paper is originally based on, which means that the equations presented in this
paper do not depend on the higher level Atomic Vortex Hypothesis based model.
However, during the course of this work it became clear that viscosity plays a
crucial role in our model, which has as consequence that an in-viscous
superfluid model is insufficient to describe the behavior of the medium.

Of course, a (viscous) superfluid model can also be described in vector notation
and since this model essentially describes a fluid/gas like medium, we can apply
continuum mechanics fluid dynamics vector calculus methods to re-derive the
Maxwell equations from the basic model. As is common practice in continuum
mechanics fluid dynamics vector theory, we can describe its dynamic behavior by
working with the medium’s flow velocity vector field39 v, with v representing
the local average bulk flow velocity. 

It should be noted that because we use continuum mechanics, the equations
presented in this paper are independent on the detailed description of the
constituents of the medium itself and that there is a lower limit with respect
to scale below which the medium can no longer be considered as a continuum. In
that case, the model is no longer applicable, which is a well-known limitation
of continuum mechanics. The Knudson number can be used to estimate this limit. 

Within the fluid dynamics domain, a scalar potential field φf and a vector
potential field 𝐀f are generally described for an incompressible fluid (∇⋅vf =
0) with a flow velocity field vf as follows40 (eq. 17-19):


(10)

where the velocity potential φf is a scalar potential field, satisfying the
Laplace equation:


(11)

and the vorticity potential 𝐀f is a solenoidal (i.e. ∇.⋅Af = 0) vector
potential field satisfying the Poisson equation:


(12)

where ωv=∇×vf  is the velocity vorticity field.

However, with this definition, the potential fields are not uniquely defined and
the boundary conditions on φf  and 𝐀f depend on the nature of the flow at the
boundary of the flow domain and on the topological properties of the flow
domain, respectively.  

We can can attempt to resolve this problem for the general case of a fluid that
is both compressible and rotational by defining a compressible irrotational
velocity field Ef for the scalar potential φf and an incompressible solenoidal
velocity field Bf and associated vorticity field ω for the vector potential 𝐀f
using the Helmholtz decomposition and negating the commonly used definition for
the velocity potential φf:


(13)


(14)

This way, the Ef and Bf fields describe flow velocity fields with a unit of
measurement in [m/s] and both the velocity potential and the velocity vorticity
potential describe fields with a unit of measurement in meters squared per
second [m2/s]. However, the primary vector field Ff thus has a unit of
measurement in [m3/s], which describes a vector field for a volumetric flow rate
or volume velocity. This can be considered as the flow velocity vector field vf
times a surface S perpendicular to vf with a surface area proportional to h2
square meters [m2], with h the physical length scale in meters [m]. This results
in the zero vector when taking the limit for the length scale h to zero, which
is obviously problematic.

So far, we have considered the general mathematical case for the Helmholtz
decomposition of any given vector field F as well as its common use in both the
electrodynamics and the fluid dynamics domains, whereby we encountered a number
of problems.  In order to resolve these problems and avoid confusion with the
various fields used thus far, let us first introduce a new set of fields along
equation (5): 


(15)

where C is our primary vector field, Π is the scalar potential or pressure, Ω is
the vector potential or angular pressure, L is the longitudinal or translational
force density and R is the rotational or angular force density. Hereby, Π and Ω
have a unit of measurement in Pascal [Pa] or Newtons per square meter [N/m2] and
L and R are in Newtons per cubic meter [N/m3]. C is in Newtons per meter [N/m]
or kilograms per second squared [kg/s2], thus representing an as of yet
undefined quantity. Further down, we will see that for the medium this unit
corresponds to the Ampere, hence the choice for using the symbol C for
“current”.

Let us now consider the 3D generalization of Newton’s second law for a substance
with a certain mass density, expressed in densities or per unit volume:


(16)

with fn the force density in [N/m3], ρ the mass density of the substance, v the
velocity field, a the acceleration field and Π the pressure or scalar potential
field in Pascal [Pa], defined as the divergence of some primary field C. Since C
should exist according to the Helmholtz decomposition and should have a unit of
measurement in [kg/s2] or [N/m], we can define C as follows:   


(17)

with η the viscosity of the substance in [kg/m-s]. This way, we obtain a full 3D
generalization of Newton’s second law per unit volume, describing not only a
longitudinal force density field L but also a rotational or angular force
density field R: 


(18)

This definition also allows us to work with the vector wave equation(3):


(19)

This is a full 3D vector wave equation, in contrast to the complex wave function
that is often used in Quantum Mechanics. With wave functions there are only two
axis, the real and the imaginary, which is simply insufficient to fully describe
phenomena in three dimensions. In other words: current Quantum Mechanics
theories lack the required dimensionality in order to be capable of fully
describing the phenomena and are therefore incomplete. 

When we divide equation (18) by mass density ρ, we obtain the velocity diffusion
equation:


(20)

with a the acceleration field in [m/s2], ν the diffusivity or kinematic
viscosity, defined by:


(21)

and Λ the volumetric acceleration field, defined by: 


(22)

This results in the diffusivity for the medium ν having a unit of measurement in
meters squared per second [m2/s] and a value equal to light speed c squared
(c2), so there is a per second [/s] difference in the unit of measurement,
suggesting that in our current models the dimensionality of certain quantities
is off by a per second. As we shall see, this has profound consequences for our
understanding of physical reality including the mass-energy equivalence
principle. 

This per second difference in units of measurement suggests that in our current
models there are a number of problems involving time derivatives that have not
been properly accounted for. When we consider that the solutions of the vector
wave equation are harmonic functions, characterized by sine and cosine functions
of time, it becomes clear how these problems could have arisen. Since the cosine
is the time derivative of the sine function and vice versa, there is only a
phase difference between the two. When we consider that all known particles
adhere to the wave-particle duality principle and have characteristic
oscillation frequencies that are very high, it becomes clear that the quantum
scale phase differentials between a force acting upon a particle and the
resulting (time delayed) acceleration of that particle are virtually
undetectable at the macroscopic level.  

Note that with this diffusion equation, the only units of measurement are the
meter and the second, which means that we have succeeded in separating the
dynamics over space and time from the substance (mass density) that’s being
diffused over space and time. In other words: with this diffusion equation we
have described the quantum characteristics of spacetime itself. 

Analogous to equation (20), we can also define a second order diffusion
equation:


(23)

which we can work out further by multiplying by mass density ρ to define the
radiosity or intensity field I in Watts per square meter [W/m2], representing a
heat flux density:


(24)

or:


(25)

From this, we can derive additional fields analogous to equation (15), which
results in fields representing power density in Watts per cubic meter [W/m3] for
the first spatial derivatives and jerk j times mass density in [N/m3-s3] for the
second spatial derivatives and thus we find that the spatial derivatives of the
intensity field I are the time derivatives of the corresponding spatial
derivatives of our primary field C.

The process of taking higher order derivatives can be continued indefinitely,
whereby for harmonic solutions we end up with the same results over and over
again, resulting in only a phase differential between subsequent results.   

An interesting detail is that the intensity field I can also be defined as:


(26)

with κ the modulus or elasticity in [Pa] or [kg/m-s2], which has a unit of
measurement that differs by a per second [/s] from the unit of measurement for
viscosity η in [Pa-s] or [kg/m-s].

This reflects the difference between elastic forces and viscous (shear) forces,
namely that the  elastic force is proportional to the amount of deformation,
while the viscous one is proportional to the rate of deformation. So, it appears
we can conclude that in physical reality there are no actual static (elastic)
forces (at the quantum level) and that deep down there are only dynamic forces
and interactions which are governed by the velocity diffusion equation(20),
whereby what we observe as static forces are in reality the time derivatives of
fundamentally viscous forces.


This brings us to the mass energy equivalence principle:


(27)

which can now alternatively be formulated by:


(28)

with L the angular momentum in [kg-m2/s] of a particle with mass m and ν the
diffusivity or kinematic viscosity. This way, the Planck-Einstein relation
becomes:


(29)

with f the characteristic oscillation frequency of the particle.


with Γ the circulation rather than the energy of a particle with mass m and ν
the diffusivity or kinematic viscosity.  

This can be related to the unusual behavior of superfluids such as 3He, which
spontaneously creates quantized vortex lines when the container holding the
liquid is put into rotation41, thus forming a quantum vortex. This is a hollow
core around which the superfluid flows along an irrotational vortex pattern
(i.e.). This flow is quantized in the sense that the circulation takes on
discrete values42. The quantum unit of circulation or quantum circulation
constant is h/m, where h is Planck’s constant and m is the mass of the
superfluid particles. 

For the medium, we can equate this quantum circulation constant to ν, the
diffusivity or kinematic viscosity, which we can now also refer to as the
quantum circulation constant, and thus we can compute the mass of an elemental
aether particle along:


(30)

(31)

which computes to approximately 7.372e-51 kg, about 20 orders of magnitude
lighter than the electron. 

When we compute the Compton wavelength for such a particle, we obtain the value
of the speed of light c, but with a unit of measurement in meters [m] rather
than velocity [m/s], while it’s associated frequency computes to 1 Hertz [Hz].
Since the Compton wavelength of a particle is equal to the wavelength of a
photon whose energy is the same as the mass of that particle along the
mass-energy equivalence principle, this puts serious question marks to the
mass-energy equivalence principle in favor of our alternative in equation (28),
whereby we conclude that the quantization that is observed in physics is not a
quantization of mass/energy, but one of angular momentum. And since angular
momentum is represented by the magnetic field, we can also conclude that it’s
the magnetic field that is quantized and that magnetic field lines are actually
irrotational hollow core vortices in a superfluid medium with a circulation
equal to the quantum circulation constant ν. 

Now let us consider the Cauchy momentum equation without external forces working
on the fluid:


(32)

(33)

with  the Cauchy stress tensor, which has a unit of measurement in [N/m2] or
[Pa] and is a central concept in the linear theory of elasticity for continuum
solid bodies in static equilibrium, when the resultant force and moment on each
axis is equal to zero. It can be demonstrated that the components of the Cauchy
stress tensor in every material point in a body satisfy the equilibrium
equations and according to the principle of conservation of angular momentum,
equilibrium requires that the summation of moments with respect to an arbitrary
point is zero, which leads to the conclusion that the stress tensor is
symmetric, thus having only six independent stress components, instead of nine.

In our model, we have only four independent stress components, namely the scalar
and vector potentials Π and Ω.

From this momentum equation, the Navier-Stokes equations can be derived, of
which the most general one without external (gravitational) forces is:


(34)

(35)

with p the pressure, I the identity tensor and  the volume, bulk or second
viscosity. This can be re-written to:


(36)

(37)

This is also a second order equation, whereby notably for the viscous term  the
order of the differential operators is reversed compared to the definition of
the second spatial derivative, the vector Laplace operator, while for the
elastic term, , the divergence of the divergence is taken. Also, a separate term
is introduced for pressure as well as a convective term, . All this not only
causes the complexity of the equations to increase dramatically while
introducing redundancy in the symmetric stress tensor, it also ignores the
fundamental symmetry between the compressible, irrotational components and the
incompressible, solenoidal components as prescribed by the Helmholtz
decomposition.

When we compare this with our proposal, we end up with two fundamentally
different approaches:

    1 A solution that fundamentally only has viscosity and one fundamental
    interaction of Nature, yields harmonic solutions c.q. builds upon
    deterministic (spherical) harmonics and provides a basis for the observed
    quantization as well as a 3D generalization of the currently used
    under-dimensioned wave functions; 2 A solution that has both viscosity as
    well as elasticity, the latter of which builds upon Brownian statistical
    mechanics and thus requires randomness and is therefore non-deterministic.  

However, with our solution so far, we have lost the description of elastic
behavior and thus our model is incomplete. This again brings us to the unusual
behavior of superfluids, which is currently described with a two-fluid theory43.
Donnely notes a/o the following: 

    1. In superfluid state, liquid helium can flow without friction. A test tube
       lowered partly into a bath of helium II will gradually fill by means of a
       thin film of liquid helium that flows without friction up the tube’s
       outer wall.
    2. There is a thermo-mechanical effect. If two containers are connected by a
       very thin tube that can block any viscous fluid, an increase in
       temperature in one container will be accompanied by a rise in pressure,
       as seen by a higher liquid level in that container.
    3. The viscous properties of liquid helium lead to a paradox. The
       oscillations of a torsion pendulum in helium II will gradually decay with
       an apparent viscosity about one-tenth that of air, but if liquid helium
       is made to flow through a very fine tube, it will do so with no
       observable pressure drop—the apparent viscosity is not only small, it is
       zero!

He also describes “second sound”, fluctuations of temperature, which according
to him “has turned out to be an incredibly valuable tool in the study of quantum
turbulence” and provides a condensed summary:

“After one of his discussions with London and inspired by the recently
discovered effects, Tisza had the idea that the Bose-condensed fraction of
helium II formed a superfluid that could pass through narrow tubes and thin
films without dissipation. The uncondensed atoms, in contrast, constituted a
normal fluid that was responsible for phenomena such as the damping of pendulums
immersed in the fluid. That revolutionary idea demanded a “two-fluid” set of
equations of motion and, among other things, predicted not only the existence of
ordinary sound—that is, fluctuations in the density of the fluid—but also
fluctuations in entropy or temperature, which were given the designation “second
sound” by Russian physicist Lev Landau. By 1938 Tisza’s and London’s papers had
at least qualitatively explained all the experimental observations available at
the time: the viscosity paradox, frictionless film flow, and the
thermo-mechanical effect.”

This leads to the question of whether or not the effects described by the
current “two-fluid” theory can also be described by the (spatial derivatives) of
the two related fields we have defined, our primary field [C] and the intensity
field [I], since we already noted that the intensity field [I] does appear to
describe elastic behavior (eq (26)), apart from a per second difference in units
of measurement, and that when dealing with harmonic functions of time, such as
those describing elemental particles, it is all too easy to get these time
derivatives mixed up, because there is only a phase differential between the [C]
and [I] fields and their respective spatial derivative fields. 

Further, since electric current can be associated with both the curl of the
magnetic field as well as with electric resistance and thus dissipation, it
seems clear that rather than associating the absence of dissipation/resistance
with the absence of viscosity, this absence should be associated with the
absence of vorticity or turbulence. 

This brings us to the idea of the “vortex sponge”44,  devised by John Bernoulli
in 1736, although in the shape of vortex tubes rather than ring vortices, which
gives rise to elastic behavior of the medium because of momentum transfer
effects arising from the fine-grained vorticity. This idea matches seamlessly to
a vortex theory of atoms45, which was developed after around 1855 a new type of
vortex theory emerged, the so-called ‘vortex sponge theory’. Instead of viewing
atoms as consisting out of small, separate vortices, this type of theory
supposed that the ether was completely filled with tiny vortices. These tiny,
close-packed vortices made up large sponge-like structures, which gave this
category of models its name. 

When we consider the basic idea that particles consist of a number (quantized)
vortex rings, we would then consider these to form such a vortex sponge,
especially in the case of crystalline materials such as silicon. Therefore, we
would associate a material substance to such a dynamic vortex sponge rather than
considering the aether itself to be universally filled with tiny vortices.  And
since in this view there are no stiff, point-like particles that bounce onto one
another in a random manner, we would also do away with Brownian statistics and
consider the interactions between the vortices to occur along harmonic functions
of time.

The question then becomes whether or not the static forces we consider at the
macroscopic level really are forces along Newton’s third:


(38)

or are actually time derivatives thereof, yank, along:


(39)

with j the jerk, the time derivative of acceleration. The derivative of force
with respect to time does not have a standard term in physics, but the term
“yank” has recently been proposed in biomechanics46.

This would also offer further insight into what inertia, resistance to "change",
actually is, because the dynamic viscous forces we have described thus far are
proportional to a rate of deformation and describe something dynamic, whereby
there is a continuous flow of mass  along quantized irrotational vortices. In a
way, this can be seen as the opposite of resistance to “chance” and could
perhaps rather be thought of as conductance of “change”.

Let's illustrate that along the rotating superfluid wherein quantum vortexes are
formed. Once a certain angular speed has been established with the rotating
container, a certain number of quantum vortices have formed and the system is in
equilibrium. In that situation, the vortices are irrotational and therefore no
vorticity nor turbulence and thus no resistance nor dissipation. In other words:
there is a steady-state situation, which could easily be confused with a
"static" situation, were it not that these the vortex lines are visible.

When we wish to increase the rotation speed of the rotating container, we must
exert a "force" and thus we introduce turbulence until a new equilibrium is
established. This way, we convert the energy we provide into the rotating
superfluid, whereby the steady state situation becomes disturbed and turbulence
is introduced, which results in more quantum vortices forming until eventually
the turbulence dies out and a new equilibrium is established. Thus, from the
outside it appears as though the rotating mass in the container resists change,
but in reality it sort of stores "change" by forming additional vortices until
there is no more turbulence.


Now let us look back at equation (20),  the velocity diffusion equation:


(40)

with a the acceleration field in [m/s2], ν the diffusivity, kinematic viscosity
or quantum circulation constant. This is an equation with only meters and
seconds and by dividing by velocity we can find that the time derivative
operator can be related to the second spatial derivative by a single constant:


(41)

which suggests space and time are indeed closely related.



-:-

Intermezzo:  my current “to do” list, some cut&pastes from discussions on
researchgate:

Where we are now is that we can describe both the quantum level as well as the
superfluid level (quantum phenomena at a macroscopic scale) with the same
equations, only different parameters like mass density and quantum circulation
constant, whereby we find that fundamentally there are only viscous forces.

What we see with superfluids is that when temperature rises (power density
increases), elastic behavior emerges, which is currently described with a
two-fluid model.

It seems that this effect can be attributed to the formation of some kind of
vortex sponge which gives rise to elastic behavior. And it also seems we can
describe the effects this creates within a continuum by the definition of fields
that are derived from the intensity field [I] rather than our primary field [C].

The fields that can be defined as the second spatial derivatives of [I] have a
unit of measurement describing the time derivative of force density, which would
be yank density.

What it appears to come down to is that within current physics Force and Yank
have been considered as one and the same thing, resulting in 3D equations that
break the fundamental symmetry demanded by the vector Laplace operator.

So, it seems that there are actually two versions of Newton's law, which have
currently been taken together into one:

1) F = m a,
2) Y = m j,

and the challenge thus comes down to figuring out which one of the two applies
where.

When we put these quantities in a table:

Action: kg-mˆ2/s.   Momentum mv: kg-m/s.    Momentum density:   kg/mˆ2-s Energy:
kg-mˆ2/sˆ2  Force    ma:       kg-m/sˆ2     Force density:          kg/mˆ2-sˆ2
Power:      kg-mˆ2/sˆ3  Yank     mj:       kg-m/sˆ3     Yank density:
kg/mˆ2-sˆ3

it also seems that additional fields can be defined to describe action density
and its spatial derivative, momentum density.

----

"Sorry, but the root of your special problem is not vector analysis. It is your
naïve assumption that you are free in selecting parts of the Navier Stokes
equation to handle special problems."

Well, I must confess I was a bit too fast by assuming that because I started
from the vector Laplace operator and all seemed to fit seamlessly, I had solved
the puzzle and that the loss of a few independent stress components was nothing
to worry about. So, guilty as charged in that respect.  

However, it was not an exercise in selecting parts of Navier Stokes equations
that met my needs, it was an attempt to derive equivalents of Navier Stokes from
vector potential theory and to align these with equivalents of Maxwell and to
derive both from one and the same equation, which turned out to represent
Newton's third in 3D.

Since I was familiar with the scalar and vector potentials used in Maxwell and I
found that the terms in the vector Laplace operator can be written out and
define fields that establish a Helmholtz decomposition, I became convinced that
this is the way it should be done. When I searched for usage of a vector
potential in fluid dynamics, I found this paper and not much more: 

https://pdfs.semanticscholar.org/9344/48b028a3a51a7567c2b441b5ca3e49ebb85c.pdf 

As I wrote in my paper, I attempted to define a primary vector field for the
Laplace operator to work on for these, since that should exist according to the
Helmholtz decomposition. It seemed that all I needed to do was negate the
definition for the scalar potential, but then the unit of measurement for the
primary field turned out to be in [mˆ3/s], denoting a volumetric flow velocity,
which results in the null vector when taking the limit of the volume to zero. So
that didn't work out very well.

After a lot of puzzling, I found a solution that involved viscosity, whereby I
found that the kinematic viscosity nu yielded a value equal to light speed
squared for the aether, but a mismatch in units of measurement by a per second,
pointing to problems here and there with time derivatives. When I realized that
this constant nu can also be seen as the quantum circulation constant, I became
convinced I'm on the right track and that the thus far mysterious properties of
superfluids (quantum phenomena on a macroscopic scale) offer the key to
unlocking the mysteries of quantum mechanics.


"The Navier Stokes equations have been derived from momentum conservation. For
an incompressible fluid we get two partial differential equations for density
and pressure. For a compressible fluid the energy balance must be considered,
which brings temperature and heat capacity into the game."

It is rather interesting that the fields I derived from my primary field [C] do
seem to describe an incompressible fluid (viscous behavior), while we seem to
have lost compressibility and that that should bring temperature and heat
capacity into the game.

My working hypothesis is that temperature is a measure of power density and has
a unit of measurement in Watts per cubic meter [W/mˆ3], but that may not be
correct since Stowe (see below) found a unit in [kg-m/sˆ3].  

I found a paper regarding superfluids, wherein it is stated that "second sound"
waves exist in a superfluid, which incorporates the propagation of fluctuations
in temperature:

https://sites.fas.harvard.edu/~phys191r/References/e1/donnelly2009.pdf 

According to Donnely, this phenomena “has turned out to be an incredibly
valuable tool in the study of quantum turbulence”.   

Thus, we have quite some hints suggesting that elastic behavior, or
compressibility, indeed has to do with the (spatial derivatives of) the
intensity field [I] I thus far payed little attention to. I've updated my
overview table and also included another primary field [Q] of which the second
spatial derivatives yield momentum density or mass flux, which I see as another
step forward.

What I think is an important detail is that the vector Laplace operator is the
3D generalization of the second spatial derivative, which would be dˆ2/dxˆ2 in
1D. This means that the 3D complexity of the vector equations we can define with
these three vector fields [Q], [C] and [I], such as the vector wave equation,
can be effortlessly reduced to one dimension to describe phenomena like for
instance the mechanical behavior of a long rod or a long thin tube filled with a
fluid.  


"The possible approximations are “incompressibility”, “ideal gas”, or even
”perfect gas” with a constant heat capacity. Another issue are the boundary
conditions inclusive external sources and sinks, which define the geometry of
the considered problem. Finally, the initial values are important.

With your approach you stay outside of the terminology used to define Navier
Stokes types of problems."

So far, I haven't solved the problem of temperature and black body radiation,
but now that I realize the importance of the intensity field [I] and it's
consequence that we have to consider yank rather than force, it seems it is only
a matter of time before we can come full circle.

First of all, it is rather interesting that the gas law also involves
quantization denoted by n:

P V = n Kb T, (eq 1)

With T the temperature in Kelvin, P the pressure, V the volume, n the number of
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
gas law, in his eq. 20 we find a relationship between electrical charge and
Boltzmann's constant: 

Kb = h/(qc), (eq 2)

with q elemental charge and h Planck's constant, which results in the conclusion
that the quantization in the gas law is related to the quantization of the
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
Planck’s empirical formula, for the energy distribution in a black-body
radiation field, cannot be derived from the kinetic theory of a gas with
Maxwellian statistics. Indeed, it is well-known that kinetic theory and
Maxwellian statistics lead to an energy distribution which is a sum of Wien-type
distributions, for a gas mixture with any number of different kinds of atoms or
molecules. But this only establishes the impossibility of so deriving Planck’s
distribution for a gas with a finite variety of atoms or molecules. To assert
the complete impossibility of so deriving Planck’s distribution it is essential
to eliminate the case of a gas with an infinite variety of atoms or molecules, i
.e . infinite in a mathematical sense, but physically, in practice, a very large
variety. The burden of the present paper is to show that this possibility cannot
be eliminated, but rather that it permits a far simpler derivation of Planck’s
energy distribution than has been given anywhere heretofore."

What is interesting, is that he found a relationship between the adiabatic index
ω and the number of degrees of freedom α of (aether) particles, which leads to
the conclusion that α must be equal to 6 and he concludes:

"Thus, the quest for a gas-like ethereal medium, satisfying Planck’s form for
the energy distribution, is directed to an ideal gas formed by an infinite
variety of particles, all having six degrees of freedom."


It is this adiabatic index which provides a relationship to heat capacity, since
it is also known as the heat capacity ratio:

https://en.wikipedia.org/wiki/Heat_capacity_ratio "In thermal physics and
thermodynamics, the heat capacity ratio, also known as the adiabatic index, the
ratio of specific heats, or Laplace's coefficient, is the ratio of the heat
capacity at constant pressure (C_P) to heat capacity at constant volume (C_V)."

So, while we clearly have not yet cracked the whole nut, it seems to me we are
on the right track towards the formulation of a "theory of everything", that
holy grail that has thus far proven to be unreachable, which I'm sure will turn
out to be attributable to ignoring the implications of the vector Laplace
operator.

Personally, I have no doubt both the weak and strong nuclear forces can be fully
accounted for by our model c.q. electromagnetic forces, once completely worked
out, and that the gravitational force also propagates through the aether, as
actually confirmed by the Michelson-Morley experiment, so that we will end up
with a model that is much, much simpler and only has one fundamental interaction
of nature. 

To illustrate the argument that the nuclear forces can be fully accounted for by
electromagnetic forces, I wholeheartedly recommend the experimental work of
David LaPoint, who shows this in his laboratory:

https://youtu.be/siMFfNhn6dk  


[end intermezzo] -:-


So far, we have shown that it is possible to derive a complete and
mathematically consistent set of fields from a single equation, the 3D
generalization of Newton’s second law, by using the LaPlace operator and working
out the terms thereof. With this equation, we can use the vector wave equation,
which has harmonic solutions, just like the wave function currently used in
Quantum Mechanics. This makes it possible to extend the current Quantum
Mechanical wave function solutions into full 3D solutions in a manner that
maintains the fundamental symmetry of the Helmholtz decomposition within a
framework of uniquely defined fields without gauge freedom. We have also shown
that we can decouple the dynamics of the medium from it’s substance, mass
density, with the velocity diffusion equation which reveals that the dynamics of
the medium are governed by a single constant ν, the quantum circulation
constant. And we have shown that we can take higher order derivatives of these
equations over and over again, resulting in only phase differentials for the
resulting vector spherical harmonic solutions.

What this comes down to is that we have come to a deeper model of physical
reality, which reveals a number of intricate relationships between various
fields defined so far, whereby the quantum circulation constant ν determines
that at the quantum level there is an intricate balance between translational
and angular momentum. This ultimately governs the possible harmonic solutions
that can exist in the shape of particles, the oscillating dynamic structures
that can be described by the vector spherical harmonics.   

This model offers a new tangible basis for theoretical physics that may
eventually very well lead to an an integrated “theory of everything”, which is
however by no means an easy task. 

So far, it has proven to be very challenging even to integrate Maxwell’s
equations with this basis in a manner that is completely consistent with the
current model and it’s units of measurement. Maxwell’s equations essentially
describe a phenomenological model that is based upon the assumption that some
kind of fundamental quantity called “charge” exists, to which a unit of
measurement in Coulombs [C] has been assigned. All of the units of measurement
within the electromagnetic domain can be derived from the Coulomb within this
model, but there is no definition of what charge actually is nor what current
actually is. Also, there is no explanation for why charge is considered to be
polarized. 

However, the model presented thus far has as big advantage that it describes a
fluid-like medium and thus we can use fluid dynamics phenomena as analogies in
our analysis. 

Let us start with Ampere’s original law to define current density J:


(42)

(43)

And let us provide an overview of the fields defined thus far, along with their
units of measurement:


Λ =   ν v Q =  τ ν ρ v C = η v =  ν ρ v I = η a = ν ρ a = ν(L+R) Λ, Q, C, I
[m3/s2]

[kg/s], [N-s/m], [J-s/m2], [Pa-s-m], [C] (charge) [kg/s2], [N/m] [J/m2], [Pa-m],
[A] (current, action flux) [kg/s3], [N/m-s], [J/m2-s], [Pa-m/s], [W/m2]
(radiosity Je , intensity I, energy flux)

S,  Σ P, Ω T, χ [m2/s2] [kg/m-s], [N-s/m2], [J-s/m3], [Pa-s], [V-s] (action
density, momentum density flux) [kg/m-s2], [N/m2], [J/m3], [Pa], [A/m], [V]
(energy density, momentum flux, force density flux, pressure) [kg/m-s3],
[N/m2-s], [J/m3-s], [Pa/s], [W/m3], [K] (power or heat density, force flux, yank
density flux, temperature)

M, Λ L, R Y, Ψ [m/s2] (a = dv/dt, acceleration) [kg/m2-s], [N-s/m3], [Pa-s/m] (ρ
v, momentum density, mass flux) [kg/m2-s2], [N/m3], [Pa/m], [A/m2], [C/m2-s] (ρ
a, force density, charge flux) [kg/m2-s3], [N/m3-s], [J/m4-s], [Pa/m-s],
[J/m4-s], [W/m4] (ρ j, yank density, current flux) J = curl R (electric current
density)


[kg/m3-s2], [N/m4], [Pa/m2], [A/m3] (d2ρ/dt2)

Table 1, overview of fields defined thus far.

P = ρ v Λ =   ν v C = η v =  ν ρ v I = η a = ν ρ a = ν(L+R) P, Λ, C, I

[kg/m2-s] [m3/s2]

[kg/s2], [N/m] [J/m2], [Pa-m], [A] [kg/s3], [N/m-s], [J/m2-s], [Pa-m/s], [W/m2]
(radiosity Je , intensity I, heat flux density) ρq,= div P, Π, Ω

[kg/m3-s], [C/m3] (dρ/dt,  ρq) [m2/s2] [kg/m-s2], [N/m2], [J/m3], [Pa], [A/m]
(momentum flux, energy density ) [kg/m-s3], [N/m2-s], [J/m3-s], [Pa/s], [W/m3]
(power density) L, R

[kg/m4-s] [m/s2] (a = dv/dt) [kg/m2-s2], [N/m3], [Pa/m],  [A/m2] (ρ a)
[kg/m2-s3], [N/m3-s], [J/m4-s], [Pa/m-s], [J/m4-s], [W/m4] (ρ j) J = curl R
(current density)


[kg/m3-s2], [N/m4], [Pa/m2], [A/m3] (d2ρ/dt2)

Table 1, overview of fields defined thus far.


This way, we would think of the electric field as being described by L, the
translational force density field, and the magnetic field as being described by
R, the angular force density field. And thus current would represent vorticity,
which aligns pretty well with observations such as Elmore’s non-radiating guided
surface wave47. From equation (42), this gives us a unit of measurement in
kilograms per second square [kg/s2] for the Ampere and we can define the Ampere
as well as the Coulomb by:

1 Ampere = 1 kilogram per second squared. (44) 1 Coulomb = 1 kilogram per
second. (45) 1 Ampere = 1 kilogram per second squared. (46) 1 Coulomb = 1
kilogram per second. (47)
 
We can subsequently define charge density as the divergence of momentum density:


(48)

(49)

resulting in a unit of measurement for charge densityin kilograms per cubic
meter per second [kg/m3-s]. 

With this definition, the charge to mass ratio of a particle results in a unit
of measurement in per second or Hertz [Hz], yielding a characteristic
longitudinal oscillation frequency for such a particle. For the electron, this
frequency computes to approximately 175.88 GHz, which falls within 10% of the
calculated spectral radiance dEν/dν in the observed cosmic background radiation
which peaks at 160.23 GHz and is calculated from a measured CMB temperature of
approximately 2.725 K48 suggesting a possible connection. 

This suggestion leads to the idea that even though we can describe the medium
itself as a superfluid, we cannot consider even the vacuum in outer space as
devoid from any particles, disturbances or (zero point) energy and thus we can
consider it to have a certain charge density ρqb0, a background charge density,
which would be depending on the material or medium we are working with, just
like the permeability and permittivity are. 

This way, we can define the electric field E as follows:


(50)

(51)

with L as defined in equation (15) and ρqb the background charge density,
resulting in a unit of measurement for the electric field E in meters per second
[m/s]. Coulomb’s law then becomes:


(52)

(53)

The electric (scalar) potential φ can subsequently be defined as:


(54)

(55)

with Π the scalar pressure in Pascal [Pa] as defined in equation (15), yielding
a unit of measurement in meters squared per second [m2/s] for the scalar
electric potential φ and thus we can define the Volt as:

1 Volt = 1 square meter per second. (56) 1 Volt = 1 square meter per second.
(57)

We can now also work out the unit of measurement for permittivity ε, which has
an SI unit in [C2/N-m2]. By substitution we find that this results in a unit of
measurement in kilograms per cubic meter [kg/m3] and we can equate the mass
density of the medium ρ to its permittivity:


(58)

(59)

For the magnetic field, we start out at the unit of measurement for permeability
μ, which is defined in SI units as Newtons per Ampere squared [N/A2]. By
substitution we find that this corresponds [m-s2/kg], the inverse of the
modulus/elasticity in [Pa] or [kg/m-s2]. The latter differs by a per second to
the unit of measurement for viscosity η in [Pa-s] or [kg/m-s], the same
difference we encountered earlier and which led us to conclude that in our
current models the dimensionality of certain quantities is off by a per second.
Therefore, we define the value of viscosity η but not its unit of measurement
by:


(60)

(61)

We can now define the magnetic field strength:


(62)

(63)

with R the angular force density in Newton per cubic meter [N/m3], resulting in
a unit of measurement for the magnetic field strength in Ampere per meter
squared [A/m2], which differs from the SI definition which is in Ampere per
meter [A/m]. 

The magnetic flux density then becomes:


(64)

(65)

and has a unit of measurement in per meter [/m].

The magnetic (vector) potential A can subsequently be defined as:


(66)

(67)

with Ω the angular vector pressure in Pascal [Pa].

This leaves us with a problem in the dimensionality of the Lorentz force,
however, which is not easily resolved in a satisfactory manner, although
dimensionally, we can resolve the problem by defining the Lorentz force as:


(68)

(69)

whereby λ is the wavelength of the particle along λ=c/f. With f=q/m we then
obtain qλ=mc.

This brings us in the situation whereby we have obtained a fluid dynamics medium
model that is capable of bridging the gap between the Quantum Mechanic and
macroscopic worlds in a deterministic manner, but leaves us with open questions
around the detailed nature of the Coulomb and Lorentz forces, especially in
relation to the nature of charged particles and their mass/charge ratios. 

However, it is clear that the irrotational vortex plays a dominant role in
magnetics and these can also form closed loop rings, which explains why magnetic
field lines are always closed. This suggests that toroidal ring models like
Parson’s49  can be integrated with our model, especially because solid spherical
harmonics can be expressed as series of toroidal harmonics and vice versa50 and
it is known that the solutions to the vector wave equation are the spherical
harmonics.

When we assume that particles can indeed be considered as consisting of a number
of closed loop hollow core vortex rings, then the physics of the vortex ring can
also be expected to provide further insight in the nature of the Lorentz force
working on charged particles. It is for example known that a vortex ring moves
forward with its own self-induced velocity v51. And since a vortex ring has two
axis of rotation, poloidal and toroidal, this could also offer an explanation
for the existence of the polarization currently attributed to charge.  
 
Either way, since all our fields are uniquely defined as solutions of the vector
Laplace equation, we can establish that with deriving all fields from equation
(18), we have eliminated “gauge freedom” and since we know these equations can
be transformed using the Galilean coordinate transform, we have also eliminated
the need for the Lorentz transform and are thus no longer bound to the universal
speed limit.

With this application of the fundamental theorem of vector calculus, we have
thus come to a revised version of the Maxwell equations that not only promises
to resolve all of the problems that have been found over the years, we also
obtain a model that is easy to interpret and can be easily simulated and
visualized with finite-difference time-domain methods (FTDT) as well. 

Now let us consider the difference between the definition we found for E and the
corresponding definition in Maxwell’s equations:


(70)

(71)

When considered from the presented perspective, this is what breaks the
fundamental result of Helmholtz’ decomposition, namely the decomposition into a
rotation free translational component and a divergence free rotational
component, since Am is not rotation free and therefore neither is its time
derivative.     

When taking the curl on both sides of this equation, we obtain the
Maxwell-Faraday equation, representing Faraday’s law of induction:


(72)

(73)

Faraday's law of induction is a basic law of electromagnetism predicting how a
magnetic field will interact with an electric circuit to produce an
electromotive force (EMF), which is thus a law that applies at the macroscopic
level. It is clear that this law should not be entangled with a model for the
medium and therefore our revision should be preferred.  

Discussion and Conclusions

We have shown that the terms in the Laplace operator can be written out to
define a complete and mathematically consistent whole of four closely related
vector fields which by definition form solutions to the vector Laplace equation,
a result that has tremendous consequences for both the analytical analysis of
the electromagnetic field as well as fluid dynamics vector theory, such as
weather forecasting, oceanography and mechanical engineering. The symmetry
between the fields thus defined is fundamental and has been mathematically
proven to be correct, so it is vital to maintain this fundamental symmetry in
our physics equations.  

We have also shown that we can decouple the dynamics of the medium from it’s
substance, mass density, with the velocity diffusion equation which reveals that
the dynamics of the medium are governed by a single constant ν, the quantum
circulation constant. And we have shown that we can take higher order
derivatives of these equations over and over again, resulting in only phase
differentials for the resulting vector spherical harmonic solutions.

Revising Maxwell equations by deriving directly from a superfluid medium model
using the Laplace operator, we have called upon vector theory for an ideal,
compressible, viscous Newtonian superfluid that has led to equations which are
known to be mathematically consistent, are known to be free of singularities and
are invariant to the Galilean transform as well. This results in an integrated
model which has only three fundamental units of measurement: mass, length and
time and also explains what “charge” is: a compression/decompression oscillation
of “charged” particles. 

As is known from fluid dynamics, these revised Maxwell equations predict three
types of wave phenomena, which we can easily relate to the observed phenomena:

    1 Longitudinal pressure waves, Tesla’s superluminal waves52 c.q. the super
    luminal longitudinal dielectric mode, which he found to propagate at a speed
    of 471,240 kilometers per second, within 0.1% of π/2 times the speed of
    light. The factor π/2 coincides with the situation whereby the theoretical
    reactance of a shorted lossless transmission line goes to infinity53 (eq
    1.2) and thus does not support an electromagnetic wave propagation mode; 

    2 “Transverse” “water” surface waves, occurring at the boundary of two media
    with different densities such as the metal surfaces of an antenna and air,
    aka the “near field”, Elmore’s non-radiating surface waves that have been
    shown to be guidable along a completely unshielded conductor54;

    3 Vortices and/or vortex rings, the “far field”, which is known to be
    quantized and to incorporate a thus far mysterious mixture of “particle” and
    “wave” properties aka “photons”, the so called “wave particle duality”
    principle. 

Even though the actual wave equations for these three wave types still need to
be derived, we can already conclude these to exist and predict a number of their
characteristics, because of the integration of the electromagnetic domain with
the fluid dynamics domain. The latter has a tremendous advantage, namely that
dynamic phenomena known to occur in fluids and gasses can be considered to also
occur in the medium.

Further Research

Theoretical While the revised Maxwell equations presented in this paper describe
the motions of the medium accurately in principle, the actual wave equations for
the three predicted wave types still need to be derived and worked out. This is
particularly complicated for the “transverse” “water” surface wave, because of
the fact that in current fluid dynamics theory the potential fields have not
been defined along the Helmholtz decomposition defined by the vector Laplacian
as we proposed, which leads to non-uniquely defined fields and associated
problems with boundary conditions. In order to derive a wave equation for the
“transverse” surface wave, the incompressibility constraint would have to be
removed from the Saint-Venant equations55 and these would subsequently need to
be fully worked out using vector calculus methods. 

Furthermore, we have also argued that Faraday's law should not be entangled with
the model for the medium, which leaves us without revised equations for
Faraday's law of induction. This leads to the question of why a DC current
trough a wire loop results in a magnetic field, but the magnetic field of a
permanent magnet does not induce a current in a wire wound around it. A similar
question arises when a (neodymium) magnet is used as an electrode in an
electrolysis experiment, which results in a vortex becoming visible in the
electrolyte above the magnet.

It is expected the answers to these questions as well as Faraday’s law of
induction can be worked out by considering the physics of the irrotational
vortex, given that we found that the current density is actually one and the
same thing as the vorticity of the medium, apart from a constant. In the absence
of external forces, a vortex evolves fairly quickly toward the irrotational flow
pattern, where the flow velocity v is inversely proportional to the distance r.
The fluid motion in a vortex creates a dynamic pressure that is lowest in the
core region, closest to the axis, and increases as one moves away from it. It is
the gradient of this pressure that forces the fluid to follow a curved path
around the axis and it is this pressure gradient that is directly related to the
velocity potential Φfd c.q. the velocity field component Efd.


Practical The revised Maxwell equations presented in this paper open the
possibilities of further considerations and research into the properties of the
dielectric and gravitational fields and associated wave phenomena. Because both
of these fields are considered as one and the same within the above presented
revised Maxwell paradigm, a wide range of possible applications become
conceivable, some of which are hardly imaginable from within the current
paradigm and/or are highly speculative while others are more straightforward.

Superluminal communication This is the most direct application of the theory
presented in this paper, which is supported by a number of sources mentioned in
the abstract, the oldest of which dates back to 1834, some theoretical
methods56,57,58,59 as well as some preliminary experimental work by the
author60. There is active and ongoing experimental research in this area.


Experiments regarding gravitational effects, such as aimed at obtaining thrust.
The Biefeld-Brown effect is an electrical phenomenon that has been the subject
of extensive research involving charging an asymmetric capacitor to high
voltages and the effect is commonly attributed to corona discharges which occur
only at the sharp electrode, which causes an imbalance in the number of positive
and negative ions created in comparison to when a symmetric capacitor is used. 

However, according to a report61 by researchers from the Army Research
Laboratory (ARL), the effects of ion wind was at least three orders of magnitude
too small to account for the observed force on the asymmetric capacitor in the
air. Instead, they proposed that the Biefeld–Brown effect may be better
explained using ion drift instead of ion wind. This was later confirmed by
researchers from the Technical University of Liberec62. 

If this is correct, then the need for an asymmetric capacitor raises the
question if the resulting diverging electric field can indeed be used to obtain
thrust by working on an electrically neutral dielectric, in this case a
dielectric consisting of air and net neutral ions, and how this results in a net
force acting upon the capacitor plates. It is known that a dielectric is always
drawn from a region of weak field toward a region of stronger field. It can be
shown that for small objects the force is proportional to the gradient of the
square of the electric field, because the induced polarization charges are
proportional to the fields and for given charges the forces are proportional to
the field as well. There will be a net force only if the square of the field is
changing from point to point, so the force is proportional to the gradient of
the square of the field63.

Another line of research in this regard has to do with the gravitational force
itself, which can be speculated to be caused by longitudinal dielectric flux,
which causes a pushing and not a pulling force. This is supported by Van
Flandern64, who determined that with a purely central pulling force and a finite
speed of gravity, the forces in a two-body system no longer point toward the
center of mass, which would make orbits unstable. The fact alone that a central
pulling gravity force requires a practically infinite speed makes clear that
pulling gravity models are untenable and recourse must be taken to a Lesagian
type of pushing gravity model. The longitudinal dielectric flux which would thus
describe gravity is probably caused by cosmic (microwave) background radiation.
If this naturally occurring flux had an arbitrary frequency spectrum,
superconductors would reflect this flux and would thus shield gravity, which
does not happen.

However, acceleration fields outside a rotating superconductor were found65,66,
which are referred to as Gravitomagnetic effects, and also anomalous
acceleration signals, anomalous gyroscope signals and Cooper pair mass excess
were found in experiments with rotating superconductors67. 

It can be speculated that the relation Stowe and Mingst found between the
characteristic oscillation frequency of the electron and the cosmic microwave
background radiation is what causes the spectrum of the gravitational flux and
that this is related to the characteristic oscillation frequencies of the
electron, neutron and proton as well. If that is the case, then the incoming
flux would resonate with the oscillating particles within the material at these
specific frequencies, which would therefore not be blocked/reflected but would
be absorbed/re-emitted along Huygens' principle.

It can further be speculated that when objects are rotated, their “clock”, the
characteristic oscillation frequency of the elemental particles making up the
material, would be influenced, causing them to deviate from the specific
frequencies they otherwise operate at. It is conceivable that this would result
in a condition whereby superconductors would indeed reflect the naturally
occurring gravitational flux, which could explain this anomaly.

Acknowledgments

This work would not have been possible without the groundbreaking work of Paul
Stowe and Barry Mingst, who succeeded in integrating the gravitational domain
with the electromagnetic domain within a single superfluid based model. It is
this integration that resolves the classic problems associated with aether based
theories, namely that because the gravitational field was considered to be
separate from the electromagnetic domain, the movements of planetary bodies
would necessarily result in measurable disturbances in the medium. When no such
disturbances were found in the Michelson-Morley experiment, the aether
hypothesis was considered as having been disproven. But because the
gravitational force is now considered to be a force caused by longitudinal
dielectric waves, which propagate trough the medium, this argument no longer
applies. And therewith there is no longer any reason to disregard an aether
based theory as a basis for theoretical physics.

This work would also not have been possible without the work of Eric Dollard,
N6KPH, who replicated a lot of Tesla’s experiments in the 1980’s. It is his
demonstrations and analysis of Tesla’s work that enabled the very consideration
of an aether based theory as an alternative to the current theoretical model.

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

In the simulator, it becomes clear that the way the staggered Yee grid is constructed 
and vectors representing a specific operator at a specific location in a space,
are stored in a numpy array in a particular manner,  as specified by
the grid.

So, what we learn is that there are two types of fields:

''' 





