(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: LoadingData *)
(* :Context: LoadingData` *)
(* :Author: Robert Atkinson *)
(* :Date: 2016-09-05 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) Robert Atkinson *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["LoadingData`", {"Utilities`"}]
Clear @ Evaluate[Context[] <> "*"]

loadMeltCurves::usage="to come"
headersFromMeltCurves::usage="to come"
replicateTagsFromMeltCurves::usage="to come"
datasetFromMeltCurves::usage="to come"
plotMeltCurves::usage="to come"
groupMeltCurvesByReplicates::usage="to come"

Begin["`Private`"]
Clear @ Evaluate[Context[] <> "*"]

(* Loads the indicated spreadsheet containing melt curve data from the CFX *)
(* Returns loaded data as an association *)
loadMeltCurves[path_, nHeaderRows_] := Module[
    { allData, headerNames, headerData, headerRows, woHeaders, temps, data },
    allData = Import[path, {"Sheets", "SYBR"}];
    headerNames = allData[[1;;nHeaderRows,1]];
    headerData = allData[[1;;nHeaderRows, 2;;]];
    headerRows = #[[1]] -> #[[2]] &/@ Transpose[{headerNames, headerData}] // Association;
    woHeaders = allData[[nHeaderRows+1;;]];
    temps = Transpose[woHeaders][[1]];
    data = Transpose[Transpose[woHeaders][[2 ;;]]];

    <|  "headerNames" -> headerNames,
        "headerRows" -> headerRows,   (* association mapping header name to header data with that name *)
        "temps" -> temps,
        "data" -> data
    |>
]

(* Returns a list of associations, one for each header in the data
 *)
headersFromMeltCurves[meltCurves_] := Module[
    { headerData, headerDataByCol, headers },
    headerData = #[[2]] & /@ Normal[meltCurves["headerRows"]];
    headerDataByCol = Transpose[headerData];
    Association[(# /. List->Rule &/@ Transpose[{meltCurves["headerNames"], #}])]& /@ headerDataByCol
]

replicateTagsFromMeltCurves[meltCurves_, makeTag_] := Module[ {},
    makeTag /@ headersFromMeltCurves[meltCurves]
]

(* Constructs a dataset from the load. Each row in the dataset has
 *      well
 *      tag
 *      color
 *      time (temp) series
 *      columns added by headerFunc
 *)
datasetFromMeltCurves[meltCurves_, makeTag_, headerFunc_] := Module[
    { tags, uniqueTags, colors, colorRules, makeDatum, tempAndValues, tsList, dataset },
    tags = replicateTagsFromMeltCurves[meltCurves, makeTag];
    uniqueTags = Union[tags] // Sort;
    colors = generateColors[uniqueTags // Length];
    colorRules = (# /. List->Rule)&/@ Transpose[{uniqueTags, colors}];
    makeDatum[header_, tag_, ts_] :=
        Association[
                {
                "well" -> header["Well"],
                "tag" -> tag,
                "color" -> (tag /. colorRules),
                "ts" -> ts
                }
            ~Join~ headerFunc[header]];
    tempAndValues = Transpose[{meltCurves["temps"], #}] & /@ Transpose[meltCurves["data"]];

    (* tsList is list of temperature series, one for each column / well *)
    tsList = TimeSeries[#, ResamplingMethod->{"Interpolation",InterpolationOrder->3}]&/@ tempAndValues;
    dataset = makeDatum @@ # & /@ Transpose[{headersFromMeltCurves[meltCurves], tags, tsList}] // Dataset;
    dataset
]

(* Plots the dataset using  the internally associated colors *)
plotMeltCurves[ds_Dataset, providedOptions : OptionsPattern[]] := Module[
    {tss, colors, localOptions, options},
    tss = ds[All,"ts"] // Normal;
    colors = ds[All, "color"] // Normal;
    localOptions = { ImageSize->Large };
    options = preempt[{providedOptions}, localOptions];
    ListLinePlot[tss, PlotStyle->colors, FilterRules[options, Options[ListLinePlot]]]
]

(* Returns a list of rules: LHS==tag, RHS==samples having that tag*)
groupMeltCurvesByReplicates[ds_Dataset] := Module[
    { tags },
    tags = ds[[All, "tag"]] // Normal // Union;
    Function[{tag}, tag -> ds[Select[#tag==tag&]]] /@ Union[tags]
]

End[] (* `Private` *)

EndPackage[]
























