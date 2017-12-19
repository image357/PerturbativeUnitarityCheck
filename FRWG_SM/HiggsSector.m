BeginPackage["HiggsSector`"];
(* includes *)

(* script control *)
RemoveTempFile::usage = "removes the temporary file, that prevents recalculation";

(* global definitions *)
massmatrix::usage = "mass matrix of the Higgs sector";
higgsvars::usage = "generic higgs variables";
modelpars::usage = "model parameters";
modelvars::usage = "model variables";
modelrealvars::usage = "model real variables";
modelcomplexvars::usage = "model complex variables";
modelgoldstones::usage = "model goldstone variables";
FRthreevertex::usage = "three-vertex Feynman rules";
FRfourvertex::usage = "four-vertex Feynman rules";
VTM::usage = "transformation matrix from higgsvars to modelvars";
HiggsSectorAssumptions::usage = "assumptions about the parameters and variables in the Higgs sector";

HD1m::usage = "first Higgs doublet expressed in real variables";
HD2m::usage = "second Higgs doublet expressed in real variables";
HD3m::usage = "third Higgs doublet expressed in real variables";

ConjugateHV::usage = "gives the conjugated list of mass and cp eigenstates";
subtrigrules::usage = "some simplification and replacement rules which mathematica can't handle";
RemoveGoldstones::usage = "removes goldstone interaction terms";

Begin["`Private`"];
(* script control *)
(* True = always recalculate everything, False = use previous calculation results if available *)
recalc = False;
tempfile = FileNameJoin[{($InputFileName // DirectoryName)<>"/temp", StringReplace[$Context,"`"->""]<>"TempFile"}];
RemoveTempFile[] := Quiet[DeleteFile[tempfile]];
If[FileExistsQ[tempfile] && !recalc, Print["using "<>tempfile]; Get[tempfile]; ,

(* private definitions *)

HiggsDoublet[a_,c_] := {a,c};
HD[a_,c_] := HiggsDoublet[a,c];
HiggsProduct[a_,b_] := Conjugate[a].b;

subtrigrules = { \
                   Cos[1/2*ArcTan[x_]] -> Sqrt[1/Sqrt[1+x^2]/2 + 1/2],   \
                   Sin[1/2*ArcTan[x_]] -> x/Sqrt[2*Sqrt[1+x^2]+2+2*x^2], \
                   Sin[2*ArcTan[x_]] -> 2*x/(1 + x^2),                   \
                   Cos[2*ArcTan[x_]] -> (1 - x^2)/(1 + x^2)              \
};

(* readable parameters with subscript *)
m11 = Subscript[m,1,1];
m22 = Subscript[m,2,2];
l1 = Subscript[\[Lambda],1];

modelpars = {m11,l1};
modelrealpars = {m11,l1};
modelcomplexpars = {};


(* general higgs doublet variables *)
a1 = Subscript[a,1];
b1 = Subscript[b,1];
c1 = Subscript[c,1];
d1 = Subscript[d,1];

higgsvars = {a1,b1,c1,d1};
zerosubhiggs = Table[higgsvars[[i]]->0, {i,1,Length[higgsvars]}];


(* model variables *)
G0 = Superscript[G,0];
GP = Superscript[G,"+"];
GM = Superscript[G,"-"];
hSM = Subscript[h,"SM"];

modelvars = {hSM,G0,GP,GM};
modelrealvars = {hSM,G0};
modelcomplexvars = {GP,GM};
modelgoldstones = {G0,GP,GM};

ConjugateHiggsVars[vars_] := vars[[{1,2,4,3}]];
ConjugateHV[mcv_] := ConjugateHiggsVars[mcv];


(* Assumptions *)
HiggsSectorAssumptions = {
    Element[modelrealpars, Reals],
    l1 > 0,
    Element[modelcomplexpars, Complexes],
    Element[higgsvars, Reals],
    Element[modelrealvars, Reals],
    Element[modelcomplexvars, Complexes]
};
Assuming[HiggsSectorAssumptions,


(* Higgs-Potential *)
HiggsPotential[HD1_] := Simplify[ \
    - m11^2 * HiggsProduct[HD1,HD1] \
    + l1 * HiggsProduct[HD1,HD1]^2 \
];


(* derivation of mass-matrices *)
VEV = Sqrt[m11^2/l1/2];

HD1m = HD[1/Sqrt[2]*(c1 + I*d1), VEV+1/Sqrt[2]*(a1 + I*b1)];

V3Dm = HiggsPotential[HD1m];

Print["HiggsSector: calculating mass matrix and vertices"];
massmatrix = Simplify[D[V3Dm, {higgsvars, 2}] /. zerosubhiggs];


(* calculate vertex factors *)
zerovertex = Simplify[V3Dm /. zerosubhiggs];
onevertex = Simplify[D[V3Dm, {higgsvars, 1}] /. zerosubhiggs]; (* is zero anyway *)
twovertex = massmatrix;
threevertex = Simplify[D[V3Dm, {higgsvars, 3}] /. zerosubhiggs];
fourvertex = Simplify[D[V3Dm, {higgsvars, 4}] /. zerosubhiggs];


(* transformation of variables *)
Print["HiggsSector: calculating variable transformations"];
ChargedFieldTransformationMatrix[dim_,indreal_,indimag_] := Table[
    If[MemberQ[{indreal,indimag},i] && MemberQ[{indreal,indimag},j],
        (* True goes here *)
        1/Sqrt[2] * If[i==indreal, If[j==indreal, 1, I], If[j==indreal, 1, -I]]
        ,
        (* False goes here *)
        If[i==j, 1, 0]
    ],
    {i,1,dim},{j,1,dim}
];


(* from c1,d1 to GP,GM *)
CFTMGPI = ChargedFieldTransformationMatrix[Length[higgsvars],3,4];
CFTMI = CFTMGPI;

(* all transformations combined *)
VTMI = CFTMI; (* VariableTransformationMatrixInverse *)
VTM = ConjugateTranspose[VTMI];


(* apply transformations to tensors *)
Print["HiggsSector: applying transformations to mass matrix"];
massmatrix = TensorContract[TensorProduct[massmatrix,Conjugate[VTM]], {1,3}];
massmatrix = TensorContract[TensorProduct[massmatrix,VTM], {1,3}];
massmatrix = Simplify[massmatrix];
twovertex = massmatrix;

Print["HiggsSector: applying transformations to three vertex"];
threevertex = Simplify[TensorContract[TensorProduct[threevertex,VTM], {1,4}], Trig -> False];
Print["HiggsSector: Step 1"];
threevertex = Simplify[TensorContract[TensorProduct[threevertex,VTM], {1,4}], Trig -> False];
Print["HiggsSector: Step 2"];
threevertex = Simplify[TensorContract[TensorProduct[threevertex,VTM], {1,4}], Trig -> False];
Print["HiggsSector: Step 3"];
threevertex = Simplify[threevertex /. subtrigrules];

Print["HiggsSector: applying transformations to four vertex"];
fourvertex = Simplify[TensorContract[TensorProduct[fourvertex,VTM], {1,5}], Trig -> False];
Print["HiggsSector: Step 1"];
fourvertex = Simplify[TensorContract[TensorProduct[fourvertex,VTM], {1,5}], Trig -> False];
Print["HiggsSector: Step 2"];
fourvertex = Simplify[TensorContract[TensorProduct[fourvertex,VTM], {1,5}], Trig -> False];
Print["HiggsSector: Step 3"];
fourvertex = Simplify[TensorContract[TensorProduct[fourvertex,VTM], {1,5}], Trig -> False];
Print["HiggsSector: Step 4"];
fourvertex = Simplify[fourvertex /. subtrigrules];


(* Feynman rules *)
Print["HiggsSector: extracting three vertex Feynman rules"];
nonzerothreevertexind = \
    DeleteDuplicates[
        Select[
            Select[
                Position[threevertex, _?(! (# === 0) &), 3],
            Length[#] == 3 &],
        ! MemberQ[#, 0] &],
    Sort[#1] == Sort[#2] &];

FRthreevertex = Table[{
        -I*Extract[threevertex, nonzerothreevertexind[[i]]],
        modelvars[[nonzerothreevertexind[[i]]]]
    }, {i,1,Length[nonzerothreevertexind]}];


Print["HiggsSector: extracting four vertex Feynman rules"];
nonzerofourvertexind = \
    DeleteDuplicates[
        Select[
            Select[
                Position[fourvertex, _?(! (# === 0) &), 4],
            Length[#] == 4 &],
        ! MemberQ[#, 0] &],
    Sort[#1] == Sort[#2] &];

FRfourvertex = Table[{
        -I*Extract[fourvertex, nonzerofourvertexind[[i]]],
        modelvars[[nonzerofourvertexind[[i]]]]
    }, {i,1,Length[nonzerofourvertexind]}];


(* remove Goldstones from Feynman rules *)
RemoveGoldstones[rules_] := Module[{expr=rules, sgs},
    sgs = Table[ToString[modelgoldstones[[i]], InputForm], {i,1,Length[modelgoldstones]}];
    Select[expr, ! Apply[Or, Table[StringContainsQ[ToString[#, InputForm], sgs[[i]]], {i, 1, Length[sgs]}]] &]
];


(* save variables to tempfile *)
Print["saving "<>tempfile];
RemoveTempFile[];
Save[tempfile, massmatrix];
Save[tempfile, modelpars];
Save[tempfile, modelvars];
Save[tempfile, modelrealvars];
Save[tempfile, modelcomplexvars];
Save[tempfile, modelgoldstones];
Save[tempfile, FRthreevertex];
Save[tempfile, FRfourvertex];
Save[tempfile, VTM];
Save[tempfile, HiggsSectorAssumptions];

Save[tempfile, HD1m];

Save[tempfile, ConjugateHiggsVars];
Save[tempfile, ConjugateHV];
Save[tempfile, subtrigrules];
Save[tempfile, RemoveGoldstones];


]; (* end of Assuming[...] *)
]; (* end of If[FileExistsQ[...] *)
End[];
EndPackage[];

