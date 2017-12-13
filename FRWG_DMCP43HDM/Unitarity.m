BeginPackage["Unitarity`"];
(* includes *)
Needs["HiggsSector`"];
$ContextPath = DeleteCases[$ContextPath, "HiggsSector`"];

(* script control *)
RemoveTempFile::usage = "removes the temporary file, that prevents recalculation";

(* global definitions *)
unitaritycond::usage = "set unitarity condition with this variable";
allvert::usage = "all possible 2to2 processes with their corresponding Feynman vertex factor";
EvaluateAllVert::usage = "numarically evaluate the 2to2 process matrix allvert corresponding to HiggsSector`modelpars";
EigenvaluesAllVert::usage = "gives the eigenvalues of the 2to2 process matrix allvert";
CheckUnitarityFunc::usage = "checks unitarity conditions";

unitaritycond = N[8*Pi];


Begin["`Private`"];
(* script control *)
(* True = always recalculate everything, False = use previous calculation results if available *)
recalc = False;
tempfile = FileNameJoin[{($InputFileName // DirectoryName)<>"/temp", StringReplace[$Context,"`"->""]<>"TempFile"}];
RemoveTempFile[] := Quiet[DeleteFile[tempfile]];
If[FileExistsQ[tempfile] && !recalc, Print["using "<>tempfile]; Get[tempfile]; ,

(* private definitions *)

vertex = HiggsSector`FRfourvertex;

vars = HiggsSector`modelvars;
cvars = HiggsSector`ConjugateHV[vars];

realvars = HiggsSector`modelrealvars;
compvars = HiggsSector`modelcomplexvars;

pars = HiggsSector`modelpars;

(* make all 2to2 vertices *)
allinitial = Flatten[Table[{i,j}, {i, 1, Length[vars]}, {j, 1, Length[vars]}], 1];
allfinal = Flatten[Table[{i,j}, {i, 1, Length[vars]}, {j, 1, Length[vars]}], 1];

(* TODO: check double occurence *)

allvert = Outer[List, allinitial, allfinal, 1];

Do[
    curinitial = allvert[[i,j,1]];
    curinitialsymb = curinitial /. x_/;NumberQ[x] :> vars[[x]];
    
    curfinal = allvert[[i,j,2]];
    curfinalsymb = curfinal /. x_/;NumberQ[x] :> cvars[[x]];
    curfinal = Join[Flatten[Position[vars, curfinalsymb[[1]]]], Flatten[Position[vars, curfinalsymb[[2]]]]];
    
    curvert = Sort[Join[curinitial, curfinal]];
    curvertsymb = curvert /. x_/;NumberQ[x] :> vars[[x]];
    
    matchvert = Select[vertex, MatchQ[#[[2]], curvertsymb] &];
    If[Length[matchvert] == 1,
        allvert[[i,j]] = matchvert[[1,1]];
        ,
        allvert[[i,j]] = 0;
    ];
    ,
    {i, Length[allinitial]}, {j, Length[allfinal]}
];

allvert = 1/2 * allvert;


EvaluateAllVert[x1_, x2_, x3_, x4_, x5_, x6_, x7_, x8_, x9_, x10_, x11_, x12_, x13_, x14_] := \
    allvert /. {pars[[1]]  :> Sqrt[x1],    \
                pars[[2]]  :> Sqrt[x2],    \
                pars[[3]]  :> x3,          \
                pars[[4]]  :> x4,          \
                pars[[5]]  :> x5,          \
                pars[[6]]  :> x6,          \
                pars[[7]]  :> x7,          \
                pars[[8]]  :> x8,          \
                pars[[9]]  :> x9,          \
                pars[[10]] :> x10,         \
                pars[[11]] :> x11 + I*x13, \
                pars[[12]] :> x12 + I*x14} ;

EvaluateAllVert[arglist_List] := EvaluateAllVert[arglist[[1]], arglist[[2]], arglist[[3]], arglist[[4]], arglist[[5]], arglist[[6]], arglist[[7]], arglist[[8]], arglist[[9]], arglist[[10]], arglist[[11]], arglist[[12]], arglist[[13]], arglist[[14]]];

EigenvaluesAllVert[x__] := Eigenvalues[N[EvaluateAllVert[x]]];

CheckUnitarityFunc[x__] := Apply[And, Map[# <= unitaritycond &, Abs[EigenvaluesAllVert[x]]]];
                


(* save variables to tempfile *)
Print["saving "<>tempfile];
RemoveTempFile[];
Save[tempfile, allvert];
Save[tempfile, EvaluateAllVert];
Save[tempfile, EigenvaluesAllVert];
Save[tempfile, CheckUnitarityFunc];


]; (* end of If[FileExistsQ[...] *)
End[];
EndPackage[];

