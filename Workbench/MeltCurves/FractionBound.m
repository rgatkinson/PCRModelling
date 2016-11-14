(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: FractionBound *)
(* :Context: FractionBound` *)
(* :Author: Robert Atkinson *)
(* :Date: 2016-09-06 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 Robert Atkinson *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["FractionBound`"]
<<"Utilities.m"
Clear @ Evaluate[Context[] <> "*"]

unitlessParameters::usage="to come"

Begin["`Private`"]
Clear @ Evaluate[Context[] <> "*"]
publishSymbol[name_] := Symbol["FractionBound`" <> symbolName[name]]
publishSymbol /@ { Atot, Btot, t, tau, deltaG, deltaH, deltaS, R, atot, btot, deltas, deltah }

(* Equilibrium Constant from Energy ----------------------------------------------------------------------------------*)
(*
We start from basic energy relations and derive an expression for DNA percent
bound as a function of temperature. We can look at this from a deltaH and deltaS
perspective. To a first approximation at least, deltaH and deltaS
are constants.
*)

eqns = {
    deltaG[t] == (deltaH - t deltaS) / 1 ,
    deltaG[t] == -R t Log[K] }
eqns2 = Eliminate[eqns, {deltaG[t]}]
kRule0 = FullSimplify[Solve[eqns2, K, Reals], R > 0][[1, 1]]

(*
We add units to that expression. That will introduce constants into the expression in
order to compensate for conversion from the units in which our parameters are typically
specified. We convert to SI base units in order to simplify as much as possible.
*)

system = "SIBase";
units = {
    t      -> UnitConvert[Quantity[t, "Kelvin"], system],
    deltaS -> UnitConvert[Quantity[deltaS, "cal / (mol Kelvin)"], system],
    deltaH -> UnitConvert[Quantity[deltaH, "kcal / mol"], system],
    R      -> UnitConvert[Quantity[1, "gas constant"], system]
}
eunits = { deltaG[t] -> UnitConvert[Quantity[deltag, "kcal/mol"], system] }
kRule1 = kRule0 /. units

(* Equilibrium Constant from Concentrations --------------------------------------------------------------------------*)
(*
    Our system is A + B <-> AB, where A & B are the individual strands, and AB is the hybridized pair. K is thus:
*)
chemicalSystem = { K == (A B) / AB, Atot == A + AB, Btot == B + AB  }
kEqn0 = Solve[Eliminate[chemicalSystem, {A, B}], K][[1, 1]] /. Rule -> Equal

(*
    We assume that the initial concentrations of A & B are the same.
*)

kEqn1 = kEqn0 /. {Btot -> Atot}
kRoot1 = kEqn1 /. Equal -> Subtract

(* Parameters --------------------------------------------------------------------------------------------------------*)
(*
We try some plausible parameters (see http://biotools.nubic.northwestern.edu/OligoCalc.html)
*)

unitizeUsingUsualUnits = { Btot -> Atot,
    Atot   -> Quantity[atot, "Molar"],
    t      -> Quantity[tau, "Celsius"],
    deltaH -> Quantity[deltah, "kcal / mol"],
    deltaS -> Quantity[deltas, "cal / (mol Kelvin)"],
    R      -> Quantity[1, "gas constant"]
}

convertUsualToSIUnits = UnitConvert[#, system] & /@ unitizeUsingUsualUnits

publishSymbol[numericalUsualParameters]
numericalUsualParameters = {
    deltah ->  778.5 (*kcal/mol*),
    deltas -> 2045.5 (*cal / (mol K)*),
    atot   -> 10^-5  (*molar*)
}
parameters = convertUsualToSIUnits /. numericalUsualParameters

(* Fraction Bound ----------------------------------------------------------------------------------------------------*)

(* We rewrite AB as a fraction rho of the common total amount of reagents, and solve for rho. *)
kEqn2 = kEqn1 /. { AB -> rho Atot }

(* We choose the appropriate root *)
rhoRule = Solve[kEqn2, rho][[1, 1]]

(* We substitute in what we know about K. *)
rhoRule1 = rhoRule /. kRule0

(*
Interesting! There remains a dependence on Atot as well as temperature and energy.

We substitute in our common parameters. rhoRule1 is the abstract fraction bound relation, which
assumes that all the unknowns are in appropriate units.

However, we don't usually like to work in those units. unitizeUsingUsualUnits is a set of rules adds
units in our usual units. convertUsualToSIUnits then converts those to SI. parameters then
plugs in some specific useful values.

We capture this in a functional dependence: fractionBound.
*)

unitlessParameters = convertUsualToSIUnits /. numericalUsualParameters // removeUnits

rhoRule3 = rhoRule1 /. {R -> (R /. unitlessParameters), t -> (t /. unitlessParameters)}
rhoRule3 = rhoRule3/. convertUsualToSIUnits
rhoRule3 = rhoRule3 // removeUnits

publishSymbol[fractionBound]
fractionBound[tautau_, deltahh_, deltass_, atottot_] := (rho /. rhoRule3) //. {tau -> tautau, deltah -> deltahh, deltas -> deltass, atot -> atottot}

End[] (* `Private` *)

EndPackage[]















































