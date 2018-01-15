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
l2 = Subscript[\[Lambda],2];
l3 = Subscript[\[Lambda],3];
l3h = OverHat[Subscript[\[Lambda],3]];
l4 = Subscript[\[Lambda],4];
l4h = OverHat[Subscript[\[Lambda],4]];
l5 = Subscript[\[Lambda],5];
l6 = Subscript[\[Lambda],6];
l8 = Subscript[\[Lambda],8];
l9 = Subscript[\[Lambda],9];

modelpars = {m11,m22,l1,l2,l3,l3h,l4,l4h,l5,l6,l8,l9};
modelrealpars = {m11,m22,l1,l2,l3,l3h,l4,l4h,l5,l6};
modelcomplexpars = {l8,l9};


(* general higgs doublet variables *)
a1 = Subscript[a,1];
b1 = Subscript[b,1];
c1 = Subscript[c,1];
d1 = Subscript[d,1];

a2 = Subscript[a,2];
b2 = Subscript[b,2];
c2 = Subscript[c,2];
d2 = Subscript[d,2];

a3 = Subscript[a,3];
b3 = Subscript[b,3];
c3 = Subscript[c,3];
d3 = Subscript[d,3];

higgsvars = {a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3};
zerosubhiggs = Table[higgsvars[[i]]->0, {i,1,Length[higgsvars]}];


(* model variables *)
G0 = Superscript[G,0];
GP = Superscript[G,"+"];
GM = Superscript[G,"-"];
hSM = Subscript[h,"SM"];
H2P = Subsuperscript[H,2,"+"];
H3P = Subsuperscript[H,3,"+"];
H2M = Subsuperscript[H,2,"-"];
H3M = Subsuperscript[H,3,"-"];
h2 = Subscript[h,2];
h3 = Subscript[h,3];
a2 = Subscript[a,2];
a3 = Subscript[a,3];

modelvars = {hSM,h2,a2,h3,a3,H2P,H2M,H3P,H3M,G0,GP,GM};
modelrealvars = {hSM,h2,a2,h3,a3,G0};
modelcomplexvars = {H2P,H2M,H3P,H3M,GP,GM};
modelgoldstones = {G0,GP,GM};

ConjugateHiggsVars[vars_] := vars[[{1,2,3,4,5,7,6,9,8,10,12,11}]];
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
l5 = 0;
HiggsPotential[HD1_,HD2_,HD3_] := Simplify[ \
    - m11^2 * HiggsProduct[HD1,HD1] \
    - m22^2 * (HiggsProduct[HD2,HD2] + HiggsProduct[HD3,HD3]) \
    + l1 * HiggsProduct[HD1,HD1]^2 \
    + l2 * (HiggsProduct[HD2,HD2]^2 + HiggsProduct[HD3,HD3]^2) \
    + l3 * HiggsProduct[HD1,HD1] * (HiggsProduct[HD2,HD2] + HiggsProduct[HD3,HD3]) \
    + l3h * HiggsProduct[HD2,HD2]*HiggsProduct[HD3,HD3] \
    + l4 * (HiggsProduct[HD1,HD2]*HiggsProduct[HD2,HD1] + HiggsProduct[HD1,HD3]*HiggsProduct[HD3,HD1]) \
    + l4h * HiggsProduct[HD2,HD3]*HiggsProduct[HD3,HD2] \
    \
    + l5 * HiggsProduct[HD3,HD1]*HiggsProduct[HD2,HD1] \
    + l6/2 * (HiggsProduct[HD2,HD1]^2 - HiggsProduct[HD1,HD3]^2) \
    + l8 * HiggsProduct[HD2,HD3]^2 \
    + l9 * HiggsProduct[HD2,HD3]*(HiggsProduct[HD2,HD2] - HiggsProduct[HD3,HD3]) \
    \
    + l5 * HiggsProduct[HD1,HD2]*HiggsProduct[HD1,HD3] \
    + l6/2 * (HiggsProduct[HD1,HD2]^2 - HiggsProduct[HD3,HD1]^2) \
    + Conjugate[l8] * HiggsProduct[HD3,HD2]^2 \
    + Conjugate[l9] * (HiggsProduct[HD2,HD2] - HiggsProduct[HD3,HD3])*HiggsProduct[HD3,HD2] \
];


(* derivation of mass-matrices *)
VEV = Sqrt[m11^2/l1/2];

HD1m = HD[1/Sqrt[2]*(c1 + I*d1), VEV+1/Sqrt[2]*(a1 + I*b1)];
HD2m = HD[1/Sqrt[2]*(c2 + I*d2), 1/Sqrt[2]*(a2 + I*b2)];
HD3m = HD[1/Sqrt[2]*(c3 + I*d3), 1/Sqrt[2]*(a3 + I*b3)];

V3Dm = HiggsPotential[HD1m,HD2m,HD3m];

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


(* from c1,d1,c2,d2,c3,d3 to GP,GM,H2P,H2M,H3P,H3M *)
CFTMGPI = ChargedFieldTransformationMatrix[Length[higgsvars],3,4];
CFTMH2PI = ChargedFieldTransformationMatrix[Length[higgsvars],7,8];
CFTMH3PI = ChargedFieldTransformationMatrix[Length[higgsvars],11,12];
CFTMI = CFTMH3PI.CFTMH2PI.CFTMGPI;

(* from a1,b1,GP,GM,a2,b2,H2P,H2M,a3,b3,H3P,H3M to hSM,a2,a3,b2,b3,H2P,H2M,H3P,H3M,G0,GP,GM *)
PermuteVariablesI1 = IdentityMatrix[Length[higgsvars]];
PermuteVariablesI1 = PermuteVariablesI1[[{1,5,9,6,10,7,8,11,12,2,3,4}]];

(* from a2,a3 to H,h *)
alpha = ArcTan[l5/l6]/2;
MDMHhI = IdentityMatrix[Length[higgsvars]];
MDMHhI[[{2,3},{2,3}]] = {{Cos[alpha],Sin[alpha]},{-Sin[alpha],Cos[alpha]}};

(* from b2,b3 to a,A *)
MDMaAI = IdentityMatrix[Length[higgsvars]];
MDMaAI[[{4,5},{4,5}]] = {{Cos[alpha],Sin[alpha]},{-Sin[alpha],Cos[alpha]}};

(* from H,h,a,A to h,a,H,A *)
PermuteVariablesI2 = IdentityMatrix[Length[higgsvars]];
PermuteVariablesI2 = PermuteVariablesI2[[{1,3,4,2,5,6,7,8,9,10,11,12}]];

(* all transformations combined *)
VTMI = PermuteVariablesI2.MDMaAI.MDMHhI.PermuteVariablesI1.CFTMI;
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
Save[tempfile, HD2m];
Save[tempfile, HD3m];

Save[tempfile, ConjugateHiggsVars];
Save[tempfile, ConjugateHV];
Save[tempfile, subtrigrules];
Save[tempfile, RemoveGoldstones];


]; (* end of Assuming[...] *)
]; (* end of If[FileExistsQ[...] *)
End[];
EndPackage[];

