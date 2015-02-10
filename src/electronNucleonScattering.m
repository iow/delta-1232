(* ::Package:: *)

BeginPackage["electronNucleonScattering`"]


chiFirst::usage =
"\t chiFirst[helicity] - for spin in positive z-axis direction
\t chiFirst[theta, phi, helicity] - for spin in positive (theta, phi) direction
theta:    polar angle of scattering \
(for in-particle use zero). \n\
phi:      azimuthal angle of scattering \
(for in-particle use zero). \n\
helicity: spin projection over (theta, phi) axis. \n
Two-component spinor for first \
particle of 2->2 process. \n\
Defined in c.m. reference frame.
This definition follows Jacob-Wick convention \
for spinors. \
(second state will be with opposite up/down \
conventions).";

chiSecond::usage =
"\t chiSecond[helicity] - for spin in negative z-axis direction
\t chiSecond[theta, phi, helicity] - for spin in negative (theta, phi) direction
theta:    polar angle of scattering \
(for in-particle use zero). \n\
phi:      azimuthal angle of scattering \
(for in-particle use zero). \n\
helicity: spin projection over (theta, phi) axis.

Two-component spinor for second \
particle of 2->2 process. \n\
Defined in c.m. reference frame.
This definition follows Jacob-Wick convention \
for spinors -- rather specific \
(first state will be with opposite up/down \
conventions, but it is not all things!).\n
For further references consult M. Jacob, G.C. Wick, \
Annals of Physics: 7, 404-428 (1959) (this paper has reprint, \
but it has minor typos). ";



polarVector::usage =
"\t polarVector[abs,theta,phi]
Converts vector from polar to Cartesian coordinates \n
abs: absolute value of the vector.
theta: polar angle.
phi: azimuthal angle.";


leptonOut::usage =
"\t leptonOut[mn, s, theta, phi, helicity]
mn:       mass of the compliment scattering particle \
(mass of the lepton is taken zero)\n\
s:        Mandelstam invariant in c.m. frame\
theta:    polar angle of scattering \
(for in-particle use zero). \n\
phi:      azimuthal angle of scattering \
(for in-particle use zero). \n\
helicity: spin projection over (theta, phi) axis.
Outcoming 1/2-spin lepton state spinor \
in certain helicity state.
Taken in the mL<<mN limit.
For using in X->eN process, where N is \
nucleon with mass \
mN (either proton or neutron)";

leptonIn::usage =
"Incoming 1/2-spin lepton state spinor.
Taken in the mL<<mN limit.
For using in eN-> process, where N is nucleon with mass \
mN (either proton or neutron)";


nucleonOut::usage =
"\t nucleonOut[mn, s, theta, phi, helicity]
mn:       mass of this nucleon scattering particle \
(mass of the compliment lepton is taken zero)\n\
s:        Mandelstam invariant in c.m. frame\
theta:    polar angle of scattering \
(for in-particle use zero). \n\
phi:      azimuthal angle of scattering \
(for in-particle use zero). \n\
helicity: spin projection over (theta, phi) axis.

Outcoming 1/2-spin nucleon state spinor in certain helicity state.
For using in X->eN process, where N is nucleon with mass \
mN (either proton or neutron). Taken in the mL<<mN limit.";

nucleonIn::usage =
	"Incoming 1/2-spin nucleon state spinor in certain helicity state.
For using in eN->X process, where N is nucleon with mass \
mN (either proton or neutron). Taken in the mL<<mN limit.";


g::usage =
"Minkowski metric tensor";

gamma::usage =
"Dirac gamma matrices (Dirac basis)";

e::usage =
"Subscript[e,0] - 4x4 unit matrix";

sliced::usage="Feynman slash notation";

dotMinkowski::usage =
"dotMinkovski[p, q] returns a convolution\
 of two four-vectors _p_ and _q_\
 with Minkowski metric";


jAngular::usage=
"  jAngular[i], i = 1, 2, 3
Spin-1 angular momentum operators";
jAngularVector::usage=
	"   jAngularVector[[i]]";


cepsilonFirst::usage=
"\t cepsilonFirst[helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along positive z-axis direction.
\t cepsilonFirst[theta, phi, helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along positive (theta, phi) axis direction.";

cepsilonSecond::usage=
"\t cepsilonSecond[helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along negative z-axis direction.
\t cepsilonSecond[theta, phi, helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along negative (theta, phi) axis direction.";


epsilonFirst::usage =
"\t Subscript[epsilonFirst[mass, s, theta, phi, helicity], mu]";

epsilonSecond::usage =
"\t Subscript[epsilonSecond[mass, s, theta, phi, helicity],mu]";

epsilonSecondConjugate::usage =
"\t Subscript[epsilonSecondConjugate[mass, s, theta, phi, helicity],mu]";


raritaOut::usage =
"  Subscript[raritaOut[mDelta, s, theta, phi, helicity],mu]";

raritaIn::usage =
"Subscript[raritaIn[mDelta, s, theta, phi, helicity],mu]";

raritaOutConjugate::usage=
	"  Subscript[raritaOutConjugate[mDelta, s, theta, phi, helicity],mu]";


conj::usage =
"User-defined complex conjugation. Assumes all symbolic variables are real, \
resolves only \"a+Ib->a-Ib\" substitution in all expressions.";


leptonCurrent::usage=
	"  Subscript[leptonCurrent[mOut_, mIn_, s_, theta_, phi_, h_], mu_]";


Begin["`Private`"];


chiFirst[+1/2] = {1, 0};
chiFirst[-1/2] = 1/2 (PauliMatrix[1] - I PauliMatrix[2]) . chiFirst[+1/2];

chiFirst[theta_, phi_, #1] = 
MatrixExp[-I phi PauliMatrix[3] / 2] .
 MatrixExp[-I theta PauliMatrix[2] / 2] .
 MatrixExp[+I phi PauliMatrix[3] / 2] .
 chiFirst[#1] // Simplify; & /@
{1/2, -1/2};


chiSecond[#1] =
MatrixExp[-I Pi PauliMatrix[2] / 2] .
 chiFirst[#1] //
 Simplify; & /@ {1/2, -1/2};

chiSecond[theta_, phi_, #1] = 
MatrixExp[-I phi PauliMatrix[3] / 2] .
  MatrixExp[-I theta PauliMatrix[2] / 2] .
  MatrixExp[+I phi PauliMatrix[3] / 2] .
  chiSecond[#1] // Simplify; & /@
{1/2, -1/2};


polarVector[abs_, theta_, phi_] =
{abs Sin[theta] Cos[phi],
 abs Sin[theta] Sin[phi],
 abs Cos[theta]};


pauliVector::usage = "Row of three Pauli matrices";
pauliVector = Table[PauliMatrix[k],{k,3}];


leptonOut[mn_, s_, theta_, phi_, #1] =
Module[{ml, k0, kAbs},
	ml = 0;
	k0 = kAbs = (s - mn^2)/(2 Sqrt[s]);
	Sqrt[k0 + ml] *
	  Join[chiFirst[theta, phi, #1],
		 2 * #1 * kAbs * chiFirst[theta, phi, #1]/(k0 + ml)
	  ]//Simplify[#, s>0 ]&
]; & /@
{1/2, -1/2};

leptonIn[mn_, s_, h_] :=
leptonOut[mn, s, 0, 0, h]

leptonIn[mn_, s_, theta_, phi_, h_] :=
leptonOut[mn, s, theta, phi, h]


nucleonOut[mn_, s_, theta_, phi_, #1] =
Module[{ml, p0, pAbs},
	ml = 0;
	p0 = (s + mn^2)/(2 Sqrt[s]);
	pAbs = (s - mn^2)/(2 Sqrt[s]);
	(-1)^(1/2 - #1) Sqrt[p0 + mn] *
	 Join[chiSecond[theta, phi, #1],
	       2 #1 pAbs chiSecond[theta, phi, #1]/(p0 + mn)
	 ] // Simplify[#, s>0 && mn>0 && s > mn^2]&
]; & /@
{1/2, -1/2};

nucleonIn[mn_, s_, l_] :=
nucleonOut[mn, s, 0, 0, l]

nucleonIn[mn_, s_, theta_, phi_, l_] :=
nucleonOut[mn, s, theta, phi, l]


g /: Subscript[g, 0, 0] = 1;
g /: Subscript[g, i_Integer, i_Integer] = -1;
g /: Subscript[g, i_Integer, j_Integer] = 0;

gamma /: Subscript[gamma,0] =
{{1, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, -1, 0},
 {0, 0, 0, -1}};
gamma /: Subscript[gamma,1] =
{{0, 0, 0, 1},
 {0, 0, 1, 0},
 {0, -1, 0, 0},
 {-1, 0, 0, 0}};
gamma /: Subscript[gamma,2] =
{{0, 0, 0, -I},
 {0, 0, I, 0},
 {0, I, 0, 0},
 {-I, 0, 0, 0}};
gamma /: Subscript[gamma,3] =
{{0, 0, 1, 0},
 {0, 0, 0, -1},
 {-1, 0, 0, 0},
 {0, 1, 0, 0}};
gamma /: Subscript[gamma,5] =
{{0, 0, 1, 0},
 {0, 0, 0, 1},
 {1, 0, 0, 0},
 {0, 1, 0, 0}};

e /: Subscript[e, 0] =
{{1, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, 1, 0},
 {0, 0, 0, 1}};

sliced[p_] =
Sum[Subscript[g,mu,nu] Subscript[p,mu]*Subscript[gamma,nu], {mu, 0, 3}, {nu, 0, 3}];

dotMinkowski[p_, q_] =
Sum[Subscript[g, mu, nu] Subscript[p, mu] Subscript[q, nu], {mu, 0, 3}, {nu, 0, 3}];


jAngular[1] =
{{0, 0, 0},
 {0, 0, -I},
 {0, I, 0}};
jAngular[2] =
{{0, 0, I},
 {0, 0, 0},
 {-I, 0, 0}};
jAngular[3] =
{{0, -I, 0},
 {I, 0, 0},
 {0, 0, 0}};

jAngularVector = Table[jAngular[k],{k,3}];





cepsilonFirst[0] =
	{0, 0, 1};
cepsilonFirst[+1] =
	1/Sqrt[2] * (jAngular[1] + I jAngular[2]) .
	cepsilonFirst[0];
cepsilonFirst[-1] =
	1/Sqrt[2] * (jAngular[1] - I jAngular[2]) .
	cepsilonFirst[0];


cepsilonFirst[theta_, phi_, #1] =
MatrixExp[-I phi jAngular[3]] .
MatrixExp[-I theta jAngular[2]] .
MatrixExp[+I phi jAngular[3]] .
cepsilonFirst[#1] //.
{Cos[phi] -> TrigToExp[Cos[phi]],
 Sin[phi] -> TrigToExp[Sin[phi]]} //
Simplify; & /@
{0, 1, -1};


cepsilonSecond[#1] =
MatrixExp[-I Pi jAngular[2]] .
cepsilonFirst[#1] //
Simplify; & /@
{0, 1, -1};

cepsilonSecond[theta_, phi_, #1] =
MatrixExp[-I phi jAngular[3]] .
MatrixExp[-I theta jAngular[2]] .
MatrixExp[+I phi jAngular[3]] .
cepsilonSecond[#1] //.
{Cos[phi] -> TrigToExp[Cos[phi]],
 Sin[phi] -> TrigToExp[Sin[phi]]} //
Simplify; & /@
{0, 1, -1};



Apply[
epsilonFirst /: Subscript[epsilonFirst[m_, s_, theta_, phi_, #1], #2] =
Prepend[cepsilonFirst[theta, phi, #1], 0][[#2 + 1]]; &,
Table[{t, mu}, {t, -1, 1, 2}, {mu, 0, 3}],
{2}];

epsilonFirst /: Subscript[epsilonFirst[m_, s_, theta_, phi_, 0], #1] =
Module[{p0, pAbs},
	p0 = (s + m^2)/(2Sqrt[s]);
	pAbs = (s - m^2)/(2Sqrt[s]);
	Prepend[p0/m cepsilonFirst[theta, phi, 0], pAbs/m][[#1 + 1]]
]; & /@
Table[mu, {mu, 0, 3}];







Apply[
epsilonSecond /: Subscript[epsilonSecond[m_, s_, theta_, phi_, #1], #2] =
Prepend[cepsilonSecond[theta, phi, #1], 0][[#2 + 1]]; &,
Table[{t, mu}, {t, -1, 1, 2}, {mu, 0, 3}],
{2}];

epsilonSecond /: Subscript[epsilonSecond[m_, s_, theta_, phi_, 0], #1] =
Module[{p0, pAbs},
	p0 = (s + m^2)/(2Sqrt[s]);
	pAbs = (s - m^2)/(2Sqrt[s]);
	Prepend[p0/m cepsilonSecond[theta, phi, 0], pAbs/m][[#1 + 1]]
]; & /@
Table[mu, {mu, 0, 3}];



conj[exp__] := exp /. {Complex[re_,im_] :> Complex[re, -im]};


Apply[
	epsilonSecondConjugate /: 
	Subscript[epsilonSecondConjugate[m_, s_, theta_, phi_, #1], #2] =
	conj[Subscript[epsilonSecond[m, s, theta, phi, #1], #2]]; &,
Table[{t, mu}, {t, -1, 1, 1}, {mu, 0, 3}],
{2}];



Function[mu,
{raritaOut /: Subscript[raritaOut[md_, s_, theta_, phi_,  3/2], mu] =
	Subscript[epsilonSecond[md, s, theta, phi, +1],mu] *
	(-1)^(1/2 - 1/2) nucleonOut[md, s, theta, phi, +1/2];,
raritaOut /: Subscript[raritaOut[md_, s_, theta_, phi_,  1/2], mu] =
	Sqrt[2/3] Subscript[epsilonSecond[md, s, theta, phi, 0],mu] *
	(-1)^(1/2 - 1/2) nucleonOut[md, s, theta, phi, 1/2] +
	Sqrt[1/3] Subscript[epsilonSecond[md,s, theta, phi, +1],mu] *
	(-1)^(1/2 + 1/2) nucleonOut[md, s, theta, phi, -1/2];,
raritaOut /: Subscript[raritaOut[md_, s_, theta_, phi_, -1/2], mu_] =
	Sqrt[2/3] Subscript[epsilonSecond[md,s, theta, phi, 0],mu] *
	(-1)^(1/2 + 1/2) nucleonOut[md, s, theta, phi, -1/2] +
	Sqrt[1/3] Subscript[epsilonSecond[md, s, theta, phi, -1],mu] *
	(-1)^(1/2 - 1/2) nucleonOut[md, s, theta, phi, +1/2];,
raritaOut /: Subscript[raritaOut[md_, s_, theta_, phi_, -3/2], mu] =
	Subscript[epsilonSecond[md, s, theta, phi, -1],mu] *
	(-1)^(1/2 + 1/2) nucleonOut[md, s, theta, phi, -1/2];}] /@
Table[mu, {mu, 0, 3, 1}];






raritaIn /:
Subscript[raritaIn[md_, s_, theta_, phi_, t_], mu_]:=
Subscript[raritaOut[md, s, theta, phi, t],mu];


Function[mu,
{raritaOutConjugate /: Subscript[raritaOutConjugate[md_, s_, theta_, phi_,  3/2], mu] =
	Subscript[epsilonSecondConjugate[md, s, theta, phi, +1],mu] *
	(-1)^(1/2 - 1/2) conj[nucleonOut[md, s, theta, phi, +1/2]];,
raritaOutConjugate /: Subscript[raritaOutConjugate[md_, s_, theta_, phi_,  1/2], mu] =
	Sqrt[2/3] Subscript[epsilonSecondConjugate[md, s, theta, phi, 0],mu] *
	(-1)^(1/2 - 1/2) conj[nucleonOut[md, s, theta, phi, 1/2]] +
	Sqrt[1/3] Subscript[epsilonSecondConjugate[md,s, theta, phi, +1],mu] *
	(-1)^(1/2 + 1/2) conj[nucleonOut[md, s, theta, phi, -1/2]];,
raritaOutConjugate /: Subscript[raritaOutConjugate[md_, s_, theta_, phi_, -1/2], mu_] =
	Sqrt[2/3] Subscript[epsilonSecondConjugate[md,s, theta, phi, 0],mu] *
	(-1)^(1/2 + 1/2) conj[nucleonOut[md, s, theta, phi, -1/2]] +
	Sqrt[1/3] Subscript[epsilonSecondConjugate[md, s, theta, phi, -1],mu] *
	(-1)^(1/2 - 1/2) conj[nucleonOut[md, s, theta, phi, +1/2]];,
raritaOutConjugate /: Subscript[raritaOutConjugate[md_, s_, theta_, phi_, -3/2], mu] =
	Subscript[epsilonSecondConjugate[md, s, theta, phi, -1],mu] *
	(-1)^(1/2 + 1/2) conj[nucleonOut[md, s, theta, phi, -1/2]];}] /@
Table[mu, {mu, 0, 3, 1}];



Apply[
leptonCurrent /: Subscript[leptonCurrent[mOut_, mIn_, s_, theta_, phi_, #1], #2] =
	Simplify[(conj[leptonOut[mOut, s, theta, phi, #1]].Subscript[gamma, 0].Subscript[gamma, #2].leptonIn[mIn,s,#1])];
leptonCurrent /: Subscript[leptonCurrent[mOut_, mIn_, s_, xi_, theta_, phi_, #1], #2] =
	Simplify[(conj[leptonOut[mOut, s, theta, phi, #1]].Subscript[gamma, 0].Subscript[gamma, #2].leptonIn[mIn, s, xi, 0, #1])]; &,
Table[{h, mu}, {h, -1/2, 1/2, 1}, {mu, 0, 3}],
{2}];


End[];


(*epsilonPhoton::usage =
"\t Subscript[epsilonPhoton[q0, qz, helicity], mu]
\t q0 - zero component;
\t qz - z-component;
Polarization 4-vector for (virtual) photon along z-axis.
Includes possible longitudinal polarization (helicity == 0).";

Subscript[epsilonPhoton[q0_,qz_,l_], mu_]:=
	If[l != 0,
		If[mu == 0,
			Return[0],
			Return[cepsilonFirst[l][[mu]]]
		],
		If[mu == 0,
			Return[1],
			Return[q0/qz cepsilonFirst[0][[mu]]]
		]
	];
epsilonPhotonConjugate::usage =
"\t Subscript[epsilonPhotonConjugate[q0, qz, helicity], mu]
\t q0 - zero component;
\t qz - z-component;
Polarization 4-vector for (virtual) photon along z-axis.
Conjugated version.
Includes possible longitudinal polarization (helicity == 0).";
Subscript[epsilonPhotonConjugate[q0_,qz_,l_], mu_]:=
	conj[Subscript[epsilonPhoton[q0, qz, l], mu]];

*)



EndPackage[]
