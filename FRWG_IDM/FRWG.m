BeginPackage["FRWG`"];
(* includes *)
Get["HiggsSector.m"];
$ContextPath = DeleteCases[$ContextPath, "HiggsSector`"];
Get["BosonSector.m"];
$ContextPath = DeleteCases[$ContextPath, "BosonSector`"];
Get["Unitarity.m"];
$ContextPath = DeleteCases[$ContextPath, "Unitarity`"];

(* script control *)
RemoveTempFiles::usage = "removes all temporary files, that prevent recalculation";
RemoveContextNames::usage = "removes context names from a given expression";
RCN::usage = "removes default context names from a given expression";
RCM::usage = "removes default context names from a given expression and displays it in matrix form";
RGN::usage = "removes default context names from a given expression without Goldstones";
RGM::usage = "removes default context names from a given expression without Goldstones and displays it in matrix form";

(* global definitions *)
(* Feynman rules (usage: FR3Vertex // RCM) *)
hFR3Vertex = HiggsSector`FRthreevertex;
hFR4Vertex = HiggsSector`FRfourvertex;
bFR3Vertex = BosonSector`FRthreevertex;
bFR4Vertex = BosonSector`FRfourvertex;
brFR3Vertex = BosonSector`rFRthreevertex;
brFR4Vertex = BosonSector`rFRfourvertex;

rFR3Vertex = Join[hFR3Vertex, brFR3Vertex];
rFR4Vertex = Join[hFR4Vertex, brFR4Vertex];
FR3Vertex = Join[hFR3Vertex, bFR3Vertex];
FR4Vertex = Join[hFR4Vertex, bFR4Vertex];


mmatrix = HiggsSector`massmatrix;                   (* higgs mass matrix *)
hpars = HiggsSector`modelpars;                      (* higgs model parameters *)
hvars = HiggsSector`modelvars;                      (* higgs model variables *)
hVTM = HiggsSector`VTM;                             (* higgs variables transformation matrix *)
HSAssumptions = HiggsSector`HiggsSectorAssumptions; (* assumptions for all parameters and variables in the Higgs sector *)

bpars = BosonSector`modelpars;
bvars = BosonSector`modelvars;
BSAssumptions = BosonSector`BosonSectorAssumptions;

mvars = Join[hvars, bvars];                         (* all field variables *)
mpars = Join[hpars, bpars];                         (* all model parameter *)

smatrix = Unitarity`allvert;
CheckUnitarity = Unitarity`CheckUnitarityFunc;
SetUnitarityCond::usage = "sets the unitarity threshold";


Begin["`Private`"];
(* script control *)
RemoveTempFiles[] := Module[{},
    HiggsSector`RemoveTempFile[];
    BosonSector`RemoveTempFile[];
    Unitarity`RemoveTempFile[];
];

RemoveContextNames[expr_, context_List] := Module[{removerule},
    removerule = Table[context[[i]] -> "", {i,1,Length[context]}];
    ToExpression[StringReplace[ToString[expr, InputForm], removerule]]
];

RCN[expr_, context_List:{"HiggsSector`Private`", "BosonSector`Private`"}] := RemoveContextNames[expr,context];
RCM[expr_] := MatrixForm[RCN[expr]];
RGN[expr_] := RCN[HiggsSector`RemoveGoldstones[expr]];
RGM[expr_] := RCM[HiggsSector`RemoveGoldstones[expr]];

(* private definitions *)
SetUnitarityCond[x_] := Module[{}, Unitarity`unitaritycond = x];


End[];
EndPackage[];

