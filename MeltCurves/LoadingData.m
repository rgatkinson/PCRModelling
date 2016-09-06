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

loadMelt::usage="to come"
headersFromLoadedMelt::usage="to come"
replicateTagsFromLoadedMelt::usage="to come"
datasetFromLoadedMelt::usage="to come"

Begin["`Private`"]
Clear @ Evaluate[Context[] <> "*"]

(* Loads the indicated spreadsheet containing melt curve data from the CFX *)
(* Returns loaded data as an association *)
loadMelt[path_, nHeaderRows_] := Module[
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
headersFromLoadedMelt[loadedMelt_] := Module[
    { headerData, headerDataByCol, headers },
    headerData = #[[2]] & /@ Normal[loadedMelt["headerRows"]];
    headerDataByCol = Transpose[headerData];
    Association[(# /. List->Rule &/@ Transpose[{loadedMelt["headerNames"], #}])]& /@ headerDataByCol
]

replicateTagsFromLoadedMelt[loadedMelt_, makeTag_] := Module[ {},
    makeTag /@ headersFromLoadedMelt[loadedMelt]
]

(* Constructs a dataset from the load. Each row in the dataset has
 *      well
 *      tag
 *      color
 *      time (temp) series
 *)
datasetFromLoadedMelt[loadedMelt_, makeTag_, headerFunc_] := Module[
    { tags, uniqueTags, colors, colorRules, makeDatum, tempAndValues, tsList, dataset },
    tags = replicateTagsFromLoadedMelt[loadedMelt, makeTag];
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
    tempAndValues = Transpose[{loadedMelt["temps"], #}] & /@ Transpose[loadedMelt["data"]];

    (* tsList is list of temperature series, one for each column / well *)
    tsList = TimeSeries[#, ResamplingMethod->{"Interpolation",InterpolationOrder->3}]&/@ tempAndValues;
    dataset = makeDatum @@ # & /@ Transpose[{headersFromLoadedMelt[loadedMelt], tags, tsList}] // Dataset;
    dataset
]

End[] (* `Private` *)

EndPackage[]