(* ::Package:: *)

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

chiFirst[h_]:=
	Switch[h,
		+1/2,
			{1, 0},
		-1/2,
			1/2 (PauliMatrix[1] - I PauliMatrix[2]).chiFirst[+1/2]
	];

chiFirst[theta_, phi_, h_]:=
	Return[MatrixExp[-I phi PauliMatrix[3] / 2]
		 . MatrixExp[-I theta PauliMatrix[2] / 2]
		 . MatrixExp[+I phi PauliMatrix[3] / 2] . chiFirst[h]
	//Simplify
	];


								
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

chiSecond[h_]:=
	Return[MatrixExp[-I Pi PauliMatrix[2] / 2]
		 . chiFirst[h]
	//Simplify
	];

chiSecond[theta_, phi_, h_]:=
	Return[MatrixExp[-I phi PauliMatrix[3] / 2]
		 . MatrixExp[-I theta PauliMatrix[2] / 2]
		 . MatrixExp[+I phi PauliMatrix[3] / 2] . chiSecond[h]
	//Simplify
	];


polarVector::usage =
"\t polarVector[abs,theta,phi]
Converts vector from polar to Cartesian coordinates \n
abs: absolute value of the vector.
theta: polar angle.
phi: azimuthal angle.";
polarVector[abs_,theta_,phi_] :=
	{abs Sin[theta] Cos[phi],abs Sin[theta] Sin[phi],abs Cos[theta]};


pauliVector::usage = "Row of three Pauli matrices.";
pauliVector:= Table[PauliMatrix[k],{k,3}];


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
Outcoming 1/2-spin lepton state spinor in certain helicity state.
Taken in the mL<<mN limit.
For using in X->eN process, where N is nucleon with mass \
mN (either proton or neutron)";
leptonOut[mn_, s_, theta_, phi_, h_]:=
	Block[
		{  
			mL = 0,
			k0 = (s - mn^2)/(2 Sqrt[s]),
			kAbs = (s - mn^2)/(2 Sqrt[s])
		},
			Return[Sqrt[k0 + mL]
					*Join[chiFirst[theta, phi, h],
					      2 * h * kAbs *
					       chiFirst[theta, phi, h]/(k0 + mL)]
					//Simplify[#, s>0 ]&
			]
	];

leptonIn::usage =
"Incoming 1/2-spin lepton state spinor.
Taken in the mL<<mN limit.
For using in eN-> process, where N is nucleon with mass \
mN (either proton or neutron)";

leptonIn[mn_, s_, h_]:=
	leptonOut[mn, s, 0, 0, h];

leptonIn[mn_, s_, theta_, phi_, h_]:=
	leptonOut[mn, s, theta, phi, h];


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

nucleonOut[mn_, s_, theta_, phi_, l_]:=
	Block[
		{  
			ml = 0,
			p0 = (s + mn^2)/(2 Sqrt[s]),
			pAbs = (s - mn^2)/(2 Sqrt[s])
		},
			Return[(-1)^(1/2 - l) Sqrt[p0 + mn]
					* Join[chiSecond[theta, phi, l],
					       2 * l * pAbs *
							chiSecond[theta, phi, l]/(p0 + mn)]
			//Simplify[#, s>0 && mn>0 && s > mn^2]&
			]
	];
nucleonIn::usage =
	"Incoming 1/2-spin nucleon state spinor in certain helicity state.
For using in eN->X process, where N is nucleon with mass \
mN (either proton or neutron). Taken in the mL<<mN limit.";
nucleonIn[mn_, s_, l_]:=
	nucleonOut[mn, s, 0, 0, l];
nucleonIn[mn_, s_, theta_, phi_, l_]:=
	nucleonOut[mn, s, theta, phi, l];


g::usage="Minkowski metric tensor";
Subscript[g, mu_,nu_] :=
	If[ mu != nu,
		(*true*)  0,
		(*false*) If[mu == 0, 1, -1]
	]; 

gamma::usage="Dirac gamma matrices (Dirac basis)";
Subscript[gamma,0]:={{1,0,0,0},{0,1,0,0},{0,0,-1,0},{0,0,0,-1}};
Subscript[gamma,1]:={{0,0,0,1},{0,0,1,0},{0,-1,0,0},{-1,0,0,0}};
Subscript[gamma,2]:={{0,0,0,-I},{0,0,I,0},{0,I,0,0},{-I,0,0,0}};
Subscript[gamma,3]:={{0,0,1,0},{0,0,0,-1},{-1,0,0,0},{0,1,0,0}};
Subscript[gamma,5]:={{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}};

e::usage="Subscript[e,0] - 4x4 unit matrix"
Subscript[e,0]:={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

sliced::usage="Feynman slash notation";
sliced[p_]:= Sum[Subscript[g,mu,nu] Subscript[p,mu] Subscript[gamma,nu],{mu, 0, 3},{nu,0,3}];
dotMinkowski[p_, q_] := Sum[Subscript[g, mu, nu] Subscript[p, mu] Subscript[q, nu], {mu, 0, 3}, {nu,0,3}];


jAngular::usage=
"  jAngular[i], i = 1, 2, 3
Spin-1 angular momentum operators";
jAngular[1] =  {{0,0,0},{0,0,-I},{0,I,0}};
jAngular[2] =  {{0,0,I},{0,0,0},{-I,0,0}};
jAngular[3] = {{0,-I,0},{I,0,0},{0,0,0}};
jAngularVector::usage=
	"   jAngularVector[[i]]";
jAngularVector:= Table[jAngular[k],{k,3}];


cepsilonFirst::usage=
"\t cepsilonFirst[helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along positive z-axis direction.
\t cepsilonFirst[theta, phi, helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along positive (theta, phi) axis direction.";

cepsilonFirst[t_]:=
	Switch[t,
		0,
			{0, 0, 1},
		+1,
			1/Sqrt[2] * (jAngular[1] + I jAngular[2]) . cepsilonFirst[0],
		-1,
			1/Sqrt[2] * (jAngular[1] - I jAngular[2]) . cepsilonFirst[0]];

cepsilonFirst[theta_, phi_, t_]:=
	Return[MatrixExp[-I phi jAngular[3]]
		 . MatrixExp[-I theta jAngular[2]]
		 . MatrixExp[+I phi jAngular[3]] . cepsilonFirst[t]
	//.{Cos[phi] -> TrigToExp[Cos[phi]], Sin[phi] -> TrigToExp[Sin[phi]]}
	//Simplify 
	];


cepsilonSecond::usage=
"\t cepsilonSecond[helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along negative z-axis direction.
\t cepsilonSecond[theta, phi, helicity] -- Returns the polarization 3-vector of particle with spin 1 and helictity _helicity_ \
along negative (theta, phi) axis direction.";

cepsilonSecond[t_]:=
	Return[MatrixExp[-I Pi jAngular[2]]
		 . cepsilonFirst[t]
	//Simplify
	];
cepsilonSecond[theta_, phi_, t_]:=
	Return[MatrixExp[-I phi jAngular[3]]
		 . MatrixExp[-I theta jAngular[2]]
		 . MatrixExp[+I phi jAngular[3]] . cepsilonSecond[t]
	//.{Cos[phi] -> TrigToExp[Cos[phi]], Sin[phi] -> TrigToExp[Sin[phi]]}
	//Simplify 
	];


epsilonFirst::usage=
"\t Subscript[epsilonFirst[mass, s, theta, phi, helicity], mu]";

Subscript[epsilonFirst[m_, s_, theta_, phi_, t_], mu_]:=
	If[t != 0,
		If[mu == 0,
			Return[0],
			Return[cepsilonFirst[theta, phi, t][[mu]]]
		],
		Block[
			{
			 p0 = (s + m^2)/(2Sqrt[s]),
			 pAbs = (s - m^2)/(2Sqrt[s])
			},
			 If[mu == 0,
				Return[pAbs/m],
				Return[p0/m cepsilonFirst[theta, phi, 0][[mu]]]
			 ]
		]
	];
			


epsilonSecond::usage=
"\t Subscript[epsilonSecond[mass, s, theta, phi, helicity],mu]";
Subscript[epsilonSecond[m_, s_, theta_, phi_, t_], mu_]:=
	If[t != 0,
		If[mu == 0,
			Return[0],
			Return[cepsilonSecond[theta, phi, t][[mu]]]
		],
		Block[
			{
			 p0 = (s + m^2)/(2Sqrt[s]),
			 pAbs = (s - m^2)/(2Sqrt[s])
			},
			 If[mu == 0,
				Return[pAbs/m],
				Return[p0/m cepsilonSecond[theta, phi, 0][[mu]]]
			 ]
		]
	];
		


conj[x_]:= Refine[Conjugate[x]];


epsilonSecondConjugate::usage=
"\t Subscript[epsilonSecondConjugate[mass, s, theta, phi, helicity],mu]";
Subscript[epsilonSecondConjugate[m_, s_, theta_, phi_, t_], mu_]:=
	conj[Subscript[epsilonSecond[m, s, theta, phi, t], mu]];


raritaOut::usage=
	"  Subscript[raritaOut[mDelta, s, theta, phi, helicity],mu]";

Subscript[raritaOut[md_, s_, theta_, phi_, t_], mu_]:=
	Return[(-1)^(3/2 - t) *
	 Switch[t,
	+3/2,
		(Subscript[epsilonSecond[md, s, theta, phi, +1],mu]
			*(-1)^(1/2 - 1/2) nucleonOut[md, s, theta, phi, +1/2]),
	+1/2,
		(Sqrt[2/3] Subscript[epsilonSecond[md, s, theta, phi, 0],mu]
			*(-1)^(1/2 - 1/2) nucleonOut[md, s, theta, phi, 1/2]
	   + Sqrt[1/3] Subscript[epsilonSecond[md,s, theta, phi, +1],mu]
			*(-1)^(1/2 + 1/2) nucleonOut[md, s, theta, phi, -1/2]),	
	-1/2,
		(Sqrt[2/3] Subscript[epsilonSecond[md,s, theta, phi, 0],mu]
			*(-1)^(1/2 + 1/2) nucleonOut[md, s, theta, phi, -1/2]
	   + Sqrt[1/3] Subscript[epsilonSecond[md, s, theta, phi, -1],mu]
			*(-1)^(1/2 - 1/2) nucleonOut[md, s, theta, phi, +1/2]),
	-3/2,
		(Subscript[epsilonSecond[md, s, theta, phi, -1],mu]
			* (-1)^(1/2 + 1/2) nucleonOut[md, s, theta, phi, -1/2])
		]
	//Simplify
	];


(*raritaIn::usage=
	"  Subscript[raritaIn[mDelta, s, theta, phi, helicity],mu]";
Subscript[raritaIn[mDelta_, s_, theta_, phi_, helicity_],mu_]:=
	Subscript[raritaOut[mDelta, s, theta, phi, helicity],mu];*)


raritaOutConjugate::usage=
	"  Subscript[raritaOutConjugate[mDelta, s, theta, phi, helicity],mu]";

Subscript[raritaOutConjugate[md_, s_, theta_, phi_, t_], mu_]:=
	Return[(-1)^(3/2 - t) *
	 Switch[t,
	+3/2,
		(Subscript[epsilonSecondConjugate[md, s, theta, phi, +1],mu]
			*(-1)^(1/2 - 1/2) conj[nucleonOut[md, s, theta, phi, +1/2]]),
	+1/2,
		(Sqrt[2/3] Subscript[epsilonSecondConjugate[md, s, theta, phi, 0],mu]
			*(-1)^(1/2 - 1/2) conj[nucleonOut[md, s, theta, phi, 1/2]]
	   + Sqrt[1/3] Subscript[epsilonSecondConjugate[md,s, theta, phi, +1],mu]
			*(-1)^(1/2 + 1/2) conj[nucleonOut[md, s, theta, phi, -1/2]]),	
	-1/2,
		(Sqrt[2/3] Subscript[epsilonSecondConjugate[md,s, theta, phi, 0],mu]
			*(-1)^(1/2 + 1/2) conj[nucleonOut[md, s, theta, phi, -1/2]]
	   + Sqrt[1/3] Subscript[epsilonSecondConjugate[md, s, theta, phi, -1],mu]
			*(-1)^(1/2 - 1/2) conj[nucleonOut[md, s, theta, phi, +1/2]]),
	-3/2,
		(Subscript[epsilonSecondConjugate[md, s, theta, phi, -1],mu]
			* (-1)^(1/2 + 1/2) conj[nucleonOut[md, s, theta, phi, -1/2]])
		]
	//Simplify
	];


leptonCurrent::usage=
	"  Subscript[leptonCurrent[mOut_, mIn_, s_, theta_, phi_, h_], mu_]";
Subscript[leptonCurrent[mOut_, mIn_, s_, theta_, phi_, h_], mu_]:=
	Simplify[(Conjugate[leptonOut[mOut, s, theta, phi, h]].Subscript[gamma, 0].Subscript[gamma, mu].leptonIn[mIn,s,h])];
