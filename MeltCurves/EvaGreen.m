(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: EvaGreen *)
(* :Context: EvaGreen` *)
(* :Author: bob *)
(* :Date: 2016-09-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2016 bob *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["EvaGreen`"]
<<"Utilities.m"
<<"FractionBound.m"
Clear @ Evaluate[Context[] <> "*"]

Begin["`Private`"]
Clear @ Evaluate[Context[] <> "*"]
publishSymbol[name_] := Symbol["EvaGreen`" <> symbolName[name]]
publishSymbol /@ { evaSignal, FractionBound, SecondaryStructure, LinearParams, ConstParams,
    secondaryStructure, evaDS, evaSS, thermoParameters, hybridize,
    alpha, beta, sigma }

Options[evaSignal] = {
    FractionBound -> fractionBound,
    SecondaryStructure -> secondaryStructure,
    LinearParams -> {},
    ConstParams -> {},
    Volume -> removeUnits[UnitConvert[Quantity[25, "\[Mu]L"], "Liters"]]
};

(* concentration is in molar *)
evaSignal[tau_, hybridize[merLength_, truncationLength_], concentration_ (*molar*), optionsSeq : OptionsPattern[]] :=
    Module[{linearized, constants, deltah, deltas, fb, volume,
        molesEachStrand, molesDuplex, molesSingleStrandedMer, molesSingleStrandedTruncation,
        basePairsPerDuplex, basesUnpairedPerDuplex,
        basePairsPerSingleStrandedMer, basesUnpairedPerSingleStrandedMer,
        basePairsPerSingleStrandedTruncation, basesUnpairdPerSingleStrandedTrunction,
        molesBasePairs, molesUnpaired
    },

        linearized = linearize[tau,#]&/@ OptionValue[LinearParams];
        constants = constize[tau, #]& /@ OptionValue[ConstParams];

        (* We retrieve thermodynamic parameters appropriate to this situation *)
        {deltah, deltas} = thermoParameters[hybridize[merLength, truncationLength]];

        (* From that we can compute the fraction bound *)
        fb = OptionValue[FractionBound][tau, deltah, deltas, concentration];

        (*  OK, now comes the hardish part. We need to think about this carefully.
            We have three populations of strands: hybridized duplexes, single-stranded mer's, and single-stranded
            truncates. Within the duplexes, we have some double stranded region and some single stranded regions.
        *)
        volume = OptionValue[Volume];  (* in liters *)
        molesEachStrand = volume * concentration;

        molesDuplex = fb * molesEachStrand;
        molesSingleStrandedMer = molesEachStrand - molesDuplex;
        molesSingleStrandedTruncation = molesEachStrand - molesDuplex;

        (* Of the tail that's not part of the hybridization complex, some amount of it will hybridize *)
        basePairsPerDuplex = truncationLength + OptionValue[SecondaryStructure][tau, merLength - truncationLength];
        basesUnpairedPerDuplex = (merLength + truncationLength) - basePairsPerDuplex * 2;

        basePairsPerSingleStrandedMer = OptionValue[SecondaryStructure][tau, merLength];
        basesUnpairedPerSingleStrandedMer = merLength - basePairsPerSingleStrandedMer * 2;

        basePairsPerSingleStrandedTruncation = OptionValue[SecondaryStructure][tau, truncationLength];
        basesUnpairdPerSingleStrandedTrunction = truncationLength - basePairsPerSingleStrandedTruncation * 2;

        molesBasePairs  = molesDuplex * basePairsPerDuplex
            + molesSingleStrandedMer * basePairsPerSingleStrandedMer
            + molesSingleStrandedTruncation * basePairsPerSingleStrandedTruncation;

        molesUnpaired = molesDuplex * basesUnpairedPerDuplex
            + molesSingleStrandedMer * basesUnpairedPerSingleStrandedMer
            + molesSingleStrandedTruncation * basesUnpairdPerSingleStrandedTrunction;

        (* The ultimate signal is a function of the amount of ds and ss nucleotides, together with background *)
        (evaDS[tau, molesBasePairs] + evaSS[tau, molesUnpaired] + evaBackground[tau]) /. linearized /. constants
    ]

(* thermoParameters returns the deltaH and deltaS symbols appropriate for the indicated hybridization *)
thermoParameters[hybridize[merLengthIn_?NumericQ, truncationLengthIn_?NumericQ]] := Module[
    {merLength, truncationLength, dh, ds},
    merLength = Round[merLengthIn];
    truncationLength = Round[truncationLengthIn];
    If [merLength==0,
        {0,0},

        dh = publishSymbol["deltaH" <> ToString[merLength] <> ToString[truncationLength]];
        ds = publishSymbol["deltaS" <> ToString[merLength] <> ToString[truncationLength]];
        { dh, ds }
        ]
    ]

(* We know some things about behavior when DNA is absent *)
secondaryStructure[tau_,0] := 0
evaDS[tau_, 0]  := 0
evaSS[tau_, 0]  := 0

(* We model the strength of EvaGreen on SS vs DS substrates as a constant ratio sigma. The fluorescence signal
   on DS substrates is proportional to the length thereof; this proportionality also in effect incorporates the
   signal gain in our instruments. The background is independent of any strand lengths. All of these relationships
   are also temperature dependent. *)

evaSS[tau_, n_] := evaDS[tau, n] / sigma[tau]
evaDS[tau_, n_] := n * alpha[tau]
evaBackground[tau_] := beta[tau]

End[] (* `Private` *)
EndPackage[]




























