(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: Utilities *)
(* :Context: Utilities` *)
(* :Author: Robert Atkinson *)
(* :Date: 2016-09-04 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 Robert Atkinson *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["Utilities`"]

(* Clear existing definitions so that this package can easily be reloaded*)
Clear["Utilities`*"]

generateColors::usage="generates a discrete set of colors"
pairTemps::usage="usage to come"
selectIndices::usage="usage to come"
pairUp::usage="usage to come"
pairUpAll::usage="usage to come"
undataset::usage="usage to come"
domain::usage="usage to come"
prune::usage="usage to come"
removeUnits::usage="removes units from an expression"
preempt::usage="usage to come"
isUnknown::usage="to come"
findUnknowns::usage="to come"
boundAboveClosed::usage="to come"
boundAboveOpen::usage="to come"
boundBelowClosed::usage="to come"
boundBelowOpen::usage="to come"
nonNeg::usage="to come"
nonPos::usage="to come"
mustNeg::usage="to come"
mustPos::usgage="to come"
symbolRules::usage="to come"
toSymbol::usage="to come"
symbolName::usage="to come"
toList::usage="to come"

Begin["`Private`"]
Clear["Utilities`Private`*"]

(*--------------------------------------------------------------------------------------------------------------------*)
(* Utilities *)
(*--------------------------------------------------------------------------------------------------------------------*)

generateColors[n_] := ColorData["Rainbow"][#/n] & /@ Range[1,n]
generateColors[n_, all_] := ConstantArray[#, all/n] & /@ generateColors[n] // Flatten

pairTemps[temps_, sample_] := Transpose[{temps, sample}]

selectIndices[list_, predicate_] := Select[Range[1, Length[list]], predicate[list[[#]]] &]

pairUp[x_, vector_] := Transpose[{x, Normal[vector]}]
pairUpAll[x_, cols_] := pairUp[x, #] & /@ Transpose[cols]

undataset[ds_] := Keys[#] /. # & /@ Normal[ds]

domain[ds_Dataset] := Module[{domains},
  domains = domain /@ ds[All, "ts"];
  {Max[domains[[All, 1]]], Min[domains[[All, 2]]]}]
domain[assoc_Association]         := domain[assoc["ts"]]
domain[series_TemporalData]       := {series["FirstTime"], series["LastTime"]}
domain[series : {__List}]         := { Min[series[[All, 1]]], Max[series[[All, 1]]]}
domain[fn_InterpolatingFunction]  := Flatten[{fn["Domain"]}]
                                      Quiet[domain[fn : Function[_,Piecewise[{{_, _ < min_}, {_, _ > max_}} , _]]] := {min, max}]
                                      Quiet[domain[fn : Function[_, Piecewise[{{_, _ <= min_}, {_, _ >= max_}} , _]]] := {min, max}]
domain[var_Symbol, thing : _] := {var} ~Join~ domain[thing]

(* prune adjusts the time domain of data *)
prune[ds_Dataset, min_, max_] := prune[#, min, max] & /@ ds
prune[assoc_Association, min_, max_] := assoc[[Complement[Keys[assoc], {"ts"}]]] ~ Append~ {"ts" -> prune[assoc["ts"], min, max]}
prune[series_TemporalData, min_, max_] := TimeSeriesWindow[series, {min, max}]
prune[{var_Symbol, lower_, upper_}, min_, max_] := {var, Max[lower, min], Min[upper, max]}
prune[series : {__List}, min_, max_] := Select[series, #[[1]] >= min && #[[1]] <= max &]
prune[fn_InterpolatingFunction, min_, max_] := Module[{lower, upper, xx},
      lower = Max[min, domain[fn][[1]]];
      upper = Min[max, domain[fn][[2]]];
      Function[{\[FormalX]}, Piecewise[{
        {Null, \[FormalX] <= lower},
        {Null, \[FormalX] >= upper}} ,
        fn[\[FormalX]]]]
    ]

removeUnits[expr_] := expr /. {Quantity[x_, y_] :> x}

preempt[newRules_List, defaultRules_List] := Module[{keys, rules},
  rules = newRules ~Join~ defaultRules;
  keys = rules[[All, 1]] // Union;
  # -> (# /. rules) & /@ keys
]
preempt[newRules___, defaultRules_List] := preempt[{newRules}, defaultRules]

(*--------------------------------------------------------------------------------------------------------------------*)
(* Unknowns *)
(*--------------------------------------------------------------------------------------------------------------------*)

(* Support for finding / quering the unknowns / free variables of an expression *)
isUnknown[ \[ExponentialE] ] := False;
isUnknown[sym_Symbol] := True;
isUnknown[sym_Symbol[t_?NumericQ]] := True
isUnknown[other_] := False;
findUnknowns[expr_] := Reap[Scan[If[isUnknown[#], Sow[#]] &, expr, Infinity]][[2, 1]] // Union

(*--------------------------------------------------------------------------------------------------------------------*)
(* Constraints *)
(*--------------------------------------------------------------------------------------------------------------------*)

boundBelowClosed[bound_, args_List] := Sequence @@ (# >= bound & /@ args)
boundAboveClosed[bound_, args_List] := Sequence @@ (# <= bound & /@ args)
boundBelowOpen[bound_, args_List] := Sequence @@ (# > bound & /@ args)
boundAboveOpen[bound_, args_List] := Sequence @@ (# < bound & /@ args)
boundBelowClosed[bound_, args : ___] := boundBelowClosed[bound, {args}]
boundAboveClosed[bound_, args : ___] := boundAboveClosed[bound, {args}]
boundBelowOpen[bound_, args : ___] := boundBelowOpen[bound, {args}]
boundAboveOpen[bound_, args : ___] := boundAboveOpen[bound, {args}]

nonNeg[args_List] := boundBelowClosed[0, args];
nonPos[args_List] := boundAboveClosed[0, args];
mustNeg[args_List] := boundAboveOpen[0, args];
mustPos[args_List] := boundBelowOpen[0, args];
nonNeg[args : ___] := nonNeg[{args}]
nonPos[args : ___] := nonPos[{args}]
mustNeg[args : ___] := mustNeg[{args}]
mustPos[args : ___] := mustPos[{args}]

(*--------------------------------------------------------------------------------------------------------------------*)
(* Symbols *)
(*--------------------------------------------------------------------------------------------------------------------*)

symbolRules[sym_Symbol -> value_] := sym -> value
symbolRules[string_String -> value_] := Symbol[string] -> value
symbolRules[rules_List] := symbolRules /@ rules
symbolRules[rules_Sequence] := symbolRules /@ rules

toSymbol[sym_Symbol] := sym
toSymbol[string_String] := Symbol[string]

symbolName[sym_Symbol] := SymbolName[sym]
symbolName[string_String] := string

toList[list_List] := list
toList[seq_Sequence] := { seq }

End[] (* `Private` *)

EndPackage[]