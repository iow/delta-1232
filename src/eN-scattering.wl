(* ::Package:: *)

chiFirst::usage =
	"Two-component spinor for first \
particle of 2->2 process. \n\
Defined in c.m. reference frame.
This definition follows Jacob-Wick convention \
for spinors. \
(second state will be with opposite up/down \
conventions).\n\n
theta:    polar angle of scattering \
(for in-particle use zero). \n\
phi:      azimuthal angle of scattering \
(for in-particle use zero). \n\
helicity: spin projection over (theta, phi) axis.";
chiFirst[theta_, phi_, helicity_] :=
	Switch[helicity,
		+1/2,
			{Cos[theta/2], Exp[I phi] Sin[theta/2]},
		-1/2,
			{-Exp[-I phi] Sin[theta/2], Cos[theta/2]},
		_,
			Message[chiFirst::nnarg, helicity]
		];
chiFirst[helicity_] := chiFirst[0, 0, helicity];
chiFirst::nnarg="Helicity `1` is not equal to +1/2 or -1/2";

								
chiSecond::usage =
	"Two-component spinor for second \
particle of 2->2 process. \n\
Defined in c.m. reference frame.
This definition follows Jacob-Wick convention \
for spinors. \
(first state will be with opposite up/down \
conventions).\n\n
theta:    polar angle of scattering \
(for in-particle use zero). \n\
phi:      azimuthal angle of scattering \
(for in-particle use zero). \n\
helicity: spin projection over (theta, phi) axis.";
chiSecond[theta_, phi_, helicity_] :=
	chiFirst[theta, phi, -helicity];
chiSecond[helicity_] := chiSecond[0, 0, helicity];


polarVector::usage =
	"Converts vector from polar to Cartesian coordinates \n
abs: absolute value of the vector.
theta: polar angle.
phi: azimuthal angle.";
polarVector[abs_,theta_,phi_] :=
	{abs Sin[theta] Cos[phi],abs Sin[theta] Sin[phi],abs Cos[theta]};


pauliVector::usage = "Row of three Pauli matrices.";
pauliVector:= Table[PauliMatrix[k],{k,3}];


leptonOut::usage =
	"Outcoming 1/2-spin lepton state spinor in certain helicity state.
Taken in the mL<<mN limit.
For using in X->eN process, where N is nucleon with mass \
mN (either proton or neutron)";
leptonOut[mN_, s_, theta_, phi_, helicity_]:=
	Block[
		{  
			mL = 0,
			k0 = (s - mN^2)/(2 Sqrt[s]),
			kAbs = (s - mN^2)/(2 Sqrt[s]),
			kVector = polarVector[(s - mN^2)/(2 Sqrt[s]),theta, phi]
		},
			Return[Sqrt[k0 + mL]
					*Join[chiFirst[theta,phi,helicity],
					      (kVector.sigmaVector/(k0 + mL)).
					       chiFirst[theta,phi,helicity]]
			]
	];
leptonIn::usage =
	"Incoming 1/2-spin lepton state spinor.
Taken in the mL<<mN limit.
For using in eN-> process, where N is nucleon with mass \
mN (either proton or neutron)";
leptonIn[mN_, s_, helicity_]:=
	leptonOut[mN, s, 0, 0, helicity];



nucleonOut::usage =
	"Outcoming 1/2-spin nucleon state spinor in certain helicity state.
For using in X->eN process, where N is nucleon with mass \
mN (either proton or neutron). Taken in the mL<<mN limit.";
nucleonOut[mN_, s_, theta_, phi_, helicity_]:=
	Block[
		{  
			mL = 0,
			p0 = (s + mN^2)/(2 Sqrt[s]),
			pAbs = (s - mN^2)/(2 Sqrt[s]),
			pVector = polarVector[(s - mN^2)/(2 Sqrt[s]),theta, phi]
		},
			Return[Sqrt[p0 + mN]
					* Join[chiSecond[theta,phi,helicity],
					      (pVector.sigmaVector/(p0 + mN)).
					       chiSecond[theta,phi,helicity]]
			]
	];
nucleonIn::usage =
	"Incoming 1/2-spin nucleon state spinor in certain helicity state.
For using in eN->X process, where N is nucleon with mass \
mN (either proton or neutron). Taken in the mL<<mN limit.";
nucleonIn[mN_, s_, helicity_]:=
	nucleonOut[mN, s, 0, 0, helicity];



Subscript[g, mu_,nu_] :=
	If[ mu != nu,
		(*true*)  0,
		(*false*) If[mu == 0, 1, -1]
	]; 
Subscript[gamma,0]={{1,0,0,0},{0,1,0,0},{0,0,-1,0},{0,0,0,-1}};
Subscript[gamma,1]={{0,0,0,1},{0,0,1,0},{0,-1,0,0},{-1,0,0,0}};
Subscript[gamma,2]={{0,0,0,-I},{0,0,I,0},{0,I,0,0},{-I,0,0,0}};
Subscript[gamma,3]={{0,0,1,0},{0,0,0,-1},{-1,0,0,0},{0,1,0,0}};
Subscript[e,0]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
Subscript[gamma,5]={{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}};

