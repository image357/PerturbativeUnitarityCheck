BeginPackage["BosonSector`"];
(* includes *)
Needs["HiggsSector`"];
$ContextPath = DeleteCases[$ContextPath, "HiggsSector`"];

(* script control *)
RemoveTempFile::usage = "removes the temporary file, that prevents recalculation";

(* global definitions *)
modelvars::usage = "model field variables";
modelpars::usage = "model parameters";
FRthreevertex::usage = "three-vertex Feynman rules";
FRfourvertex::usage = "four-vertex Feynman rules";
rFRthreevertex::usage = "more readable three-vertex Feynman rules";
rFRfourvertex::usage = "more readable four-vertex Feynman rules";
BosonSectorAssumptions::usage = "assumptions about the parameters and variables in the Boson sector";

Begin["`Private`"];
(* script control *)
(* True = always recalculate everything, False = use previous calculation results if available *)
recalc = False;
tempfile = FileNameJoin[{($InputFileName // DirectoryName)<>"/temp", StringReplace[$Context,"`"->""]<>"TempFile"}];
RemoveTempFile[] := Quiet[DeleteFile[tempfile]];
If[FileExistsQ[tempfile] && !recalc, Print["using "<>tempfile]; Get[tempfile]; , 

(* private definitions *)
(* readable variable names *)
gw = Subscript[g,"W"];
gh = Subscript[g,"H"];
thetaw = Subscript[\[Theta],"W"];
eem = Subscript[e,"EM"];

modelpars = {gw,gh,thetaw,eem};
modelrealpars = {gw,gh,thetaw,eem};

(* pauli matrices *)
pauli1 = {{0,1},{1,0}};
pauli2 = {{0,-I},{I,0}};
pauli3 = {{1,0},{0,-1}};
pauli = {pauli1, pauli2, pauli3};

W11 = Subscript["W",1,1];
W12 = Subscript["W",1,2];
W13 = Subscript["W",1,3];
W14 = Subscript["W",1,4];
W1 = {W11,W12,W13,W14};

W21 = Subscript["W",2,1];
W22 = Subscript["W",2,2];
W23 = Subscript["W",2,3];
W24 = Subscript["W",2,4];
W2 = {W21,W22,W23,W24};

W31 = Subscript["W",3,1];
W32 = Subscript["W",3,2];
W33 = Subscript["W",3,3];
W34 = Subscript["W",3,4];
W3 = {W31,W32,W33,W34};

W = {W1,W2,W3};

B1 = Subscript["B",1];
B2 = Subscript["B",2];
B3 = Subscript["B",3];
B4 = Subscript["B",4];
B = {B1,B2,B3,B4};

WP1 = Subsuperscript["W",1,"+"];
WP2 = Subsuperscript["W",2,"+"];
WP3 = Subsuperscript["W",3,"+"];
WP4 = Subsuperscript["W",4,"+"];

WP = {WP1, WP2, WP3, WP4};

WM1 = Subsuperscript["W",1,"-"];
WM2 = Subsuperscript["W",2,"-"];
WM3 = Subsuperscript["W",3,"-"];
WM4 = Subsuperscript["W",4,"-"];

WM = {WM1, WM2, WM3, WM4};

Z1 = Subscript["Z",1];
Z2 = Subscript["Z",2];
Z3 = Subscript["Z",3];
Z4 = Subscript["Z",4];
Z = {Z1,Z2,Z3,Z4};

A1 = Subscript["A",1];
A2 = Subscript["A",2];
A3 = Subscript["A",3];
A4 = Subscript["A",4];
A = {A1,A2,A3,A4};

vectorvars = Join[W1,W2,W3,B,Z,A,WP,WM];
vectorrealvars = Join[W1,W2,W3,B,Z,A];
vectorcomplexvars = Join[WP,WM];
modelvars = Join[WP, WM, Z, A];

(* metric *)
gyv = DiagonalMatrix[{1,-1,-1,-1}];
metric = gyv;


(* Assumptions *)
BosonSectorAssumptions = {
    Element[vectorrealvars, Reals],
    Element[vectorcomplexvars, Complexes],
    Element[modelrealpars, Reals]
};
Assuming[Join[BosonSectorAssumptions,HiggsSector`HiggsSectorAssumptions],

CovariantDerivative[hd_,w_,b_,vars_List:HiggsSector`higgsvars] := \
    + { \
        Table[D[hd[[1]],{{t,x,y,z}}[[1,i]],NonConstants->vars], {i,1,4}], \
        Table[D[hd[[2]],{{t,x,y,z}}[[1,i]],NonConstants->vars], {i,1,4}]  \
      } \
    + I*gw/2*Sum[TensorProduct[pauli[[i]].hd, metric.(w[[i]])], {i,1,3}] \
    + I*gh/2*TensorProduct[hd,metric.b]
;

KineticTerm[hd_,w_,b_,vars_List:HiggsSector`higgsvars,complexvars_List:{}] := Simplify[ComplexExpand[ \
    Sum[ \
            Conjugate[CovariantDerivative[hd,w,b,vars]][[i,j]] * \
            TensorContract[TensorProduct[CovariantDerivative[hd,w,b,vars],metric],{{2,3}}][[i,j]], \
            \
            {i,1,2}, {j,1,4}
       ]
, complexvars]];

KineticLagrangian[HD1_,w_,b_,vars_List:HiggsSector`higgsvars,complexvars_List:{}] := \
    KineticTerm[HD1,w,b,vars,complexvars] \
;

subtrigrules = HiggsSector`subtrigrules;
transformedrealvars = Simplify[(HiggsSector`VTM.HiggsSector`modelvars) /. subtrigrules, Assumptions->HiggsSector`HiggsSectorAssumptions];
subtransformedvars = Table[HiggsSector`higgsvars[[i]] -> transformedrealvars[[i]], {i,1,Length[HiggsSector`modelvars]}];

(* kinetic Lagrangian *)
K3Dm = KineticLagrangian[HiggsSector`HD1m, W, B, HiggsSector`higgsvars];

(* electroweak mixing *)
subchargedvb1 = Table[W1[[i]] ->   (WP[[i]]+WM[[i]])/Sqrt[2], {i,1,4}];
subchargedvb2 = Table[W2[[i]] -> I*(WP[[i]]-WM[[i]])/Sqrt[2], {i,1,4}];
subneutralvb3 = Table[W3[[i]] ->  Cos[thetaw]*Z[[i]] + Sin[thetaw]*A[[i]], {i,1,4}];
subneutralvb4 = Table[B[[i]]  -> -Sin[thetaw]*Z[[i]] + Cos[thetaw]*A[[i]], {i,1,4}];
subewmixing = Join[subchargedvb1,subchargedvb2,subneutralvb3,subneutralvb4];
subgtoeem = {gw -> eem/Sin[thetaw], gh -> eem/Cos[thetaw]};

K3Dm = (K3Dm /. subewmixing) /. subgtoeem;

(* mass and cp eigenstates *)
HK3Dm = K3Dm /. subtransformedvars;
HK3Dm = HK3Dm /. HoldPattern[D[s1_,s2_,___]] :> D[s1,s2,NonConstants->HiggsSector`modelvars];
K3Dm = HK3Dm;

(* write derivatives in short form *)
GetDSeq[var_] := Module[{pos,temp},
    pos = Position[{t,x,y,z},var];
    temp = {0, 0, 0, 0};
    temp[[pos[[1, 1]]]] = 1;
    Apply[Sequence, temp]
];
K3Dm = K3Dm /. HoldPattern[D[s1_,s2_,___]] :> Derivative[GetDSeq[s2]][s1];

(* make variables to functions in order to have distinct behaviour of D[...] on variables and their derivatives *)
K3Dm = K3Dm /. var_/;MemberQ[HiggsSector`modelvars, var] :> var[t,x,y,z];
K3Dm = K3Dm /. HoldPattern[Derivative[n__][var_[arg__]]] :> Derivative[n][var][arg];
higgsfuncvars = Join[
    Table[HiggsSector`modelvars[[i]][t,x,y,z], {i, 1, Length[HiggsSector`modelvars]}],
    Flatten[Table[D[HiggsSector`modelvars[[i]][t,x,y,z], {{t,x,y,z}}], {i, 1, Length[HiggsSector`modelvars]}]]
];
modelfuncvars = Join[higgsfuncvars, WP, WM, Z, A];

(* Feynman rules *)
zerosubmodel = Table[modelfuncvars[[i]] -> 0, {i,1,Length[modelfuncvars]}];

zerovertex = Simplify[K3Dm /. zerosubmodel];
onevertex = Simplify[D[K3Dm, {modelfuncvars, 1}] /. zerosubmodel]; (* is zero anyway *)
twovertex = Simplify[D[K3Dm, {modelfuncvars, 2}] /. zerosubmodel];
threevertex = Simplify[D[K3Dm, {modelfuncvars, 3}] /. zerosubmodel];
fourvertex = Simplify[D[K3Dm, {modelfuncvars, 4}] /. zerosubmodel];

nonzerothreevertexind = \
    DeleteDuplicates[
        Select[
            Select[
                Position[threevertex, _?(! (# === 0) &), 3],
            Length[#] == 3 &],
        ! MemberQ[#, 0] &],
    Sort[#1] == Sort[#2] &];

nonzerofourvertexind = \
    DeleteDuplicates[
        Select[
            Select[
                Position[fourvertex, _?(! (# === 0) &), 4],
            Length[#] == 4 &],
        ! MemberQ[#, 0] &],
    Sort[#1] == Sort[#2] &];

FRthreevertex = Table[{
        I*Extract[threevertex, nonzerothreevertexind[[i]]],
        modelfuncvars[[nonzerothreevertexind[[i]]]]
    }, {i,1,Length[nonzerothreevertexind]}];

FRfourvertex = Table[{
        I*Extract[fourvertex, nonzerofourvertexind[[i]]],
        modelfuncvars[[nonzerofourvertexind[[i]]]]
    }, {i,1,Length[nonzerofourvertexind]}];

(* remove function arguments *)
FRthreevertex = FRthreevertex /. var_[t,x,y,z]/;MemberQ[HiggsSector`modelvars, var] :> var;
FRthreevertex = FRthreevertex /. HoldPattern[Derivative[n__][var_][t,x,y,z]/;MemberQ[HiggsSector`modelvars, var]] :> Derivative[n][var];
FRfourvertex = FRfourvertex /. var_[t,x,y,z]/;MemberQ[HiggsSector`modelvars, var] :> var;
FRfourvertex = FRfourvertex /. HoldPattern[Derivative[n__][var_][t,x,y,z]/;MemberQ[HiggsSector`modelvars, var]] :> Derivative[n][var];

(* trigonometirc optimisation *)
FRthreevertex = Simplify[ExpToTrig[Simplify[TrigToExp[FRthreevertex]]]];
FRfourvertex = Simplify[ExpToTrig[Simplify[TrigToExp[FRfourvertex]]]];

(* make Feynman rules more readable *)
hvectorvars = Flatten[Table[D[HiggsSector`modelvars[[i]][t,x,y,z], {{t,x,y,z}}], {i, 1, Length[HiggsSector`modelvars]}]] /. HoldPattern[Derivative[n__][var_][t,x,y,z]/;MemberQ[HiggsSector`modelvars, var]] :> Derivative[n][var];
bvectorvars = Join[WP, WM, Z, A];
hbvvars = Join[hvectorvars, bvectorvars];
zerovvars = Table[hbvvars[[4*i+1]], {i,0,Length[hbvvars]/4-1}];

subhrzvvars = Table[zerovvars[[i]] -> "\[PartialD]"*HiggsSector`modelvars[[i]], {i, 1, Length[HiggsSector`modelvars]}];
submu = {
            "\[PartialD]" -> Subscript["\[PartialD]", \[Mu]],
            bvectorvars[[0*4+1]] -> Subsuperscript["W", \[Mu], "+"],
            bvectorvars[[1*4+1]] -> Subsuperscript["W", \[Mu], "-"],
            bvectorvars[[2*4+1]] -> Subscript["Z", \[Mu]],
            bvectorvars[[3*4+1]] -> Subscript["A", \[Mu]]
        };
subnu = {
            "\[PartialD]" -> Subscript["\[PartialD]", \[Nu]],
            bvectorvars[[0*4+1]] -> Subsuperscript["W", \[Nu], "+"],
            bvectorvars[[1*4+1]] -> Subsuperscript["W", \[Nu], "-"],
            bvectorvars[[2*4+1]] -> Subscript["Z", \[Nu]],
            bvectorvars[[3*4+1]] -> Subscript["A", \[Nu]]
        }; 

rFRthreevertex = FRthreevertex;
For[ind=1, ind<=Length[rFRthreevertex], ind++,
    curcol = rFRthreevertex[[ind]];
    curvars = curcol[[2]];
    isvectorlist = Map[MemberQ[hbvvars, #] &, curvars];
    isvector = Apply[Or, isvectorlist];
    If[isvector,
        (* True goes here *)
        iszerovector = Apply[Or, Map[MemberQ[zerovvars, #] &, curvars]];
        If[iszerovector,
            (* True goes here *)
            (* insert metric in vertex factor *)
            curcol[[1]] = StringForm["`` ``", curcol[[1]], Superscript[g,\[Mu] \[Nu]]];
            (* replace derivative by partial derivative *)
            curcol[[2]] = Map[HoldForm, curcol[[2]]] /. subhrzvvars;
            
            (* add mu and nu contraction indices *)
            curcol[[2,1]] = curcol[[2,1]] /. submu;
            For[ind2=2, ind2<=Length[curcol[[2]]], ind2++,
                isaftervector = Apply[Or, isvectorlist[[1;;(ind2-1)]]];
                If[isaftervector,
                    curcol[[2,ind2]] = curcol[[2,ind2]] /. subnu
                    ,
                    curcol[[2,ind2]] = curcol[[2,ind2]] /. submu
                ];
            ];
             
            rFRthreevertex[[ind]] = curcol;
            ,
            (* False goes here *)
            rFRthreevertex = Drop[rFRthreevertex, {ind}];
            ind = ind-1;
        ];
        ,
        (* False goes here *)
        False;
    ];
];

rFRfourvertex = FRfourvertex;
For[ind=1, ind<=Length[rFRfourvertex], ind++,
    curcol = rFRfourvertex[[ind]];
    curvars = curcol[[2]];
    isvectorlist = Map[MemberQ[hbvvars, #] &, curvars];
    isvector = Apply[Or, isvectorlist];
    If[isvector,
        (* True goes here *)
        iszerovector = Apply[Or, Map[MemberQ[zerovvars, #] &, curvars]];
        If[iszerovector,
            (* True goes here *)
            (* insert metric in vertex factor *)
            curcol[[1]] = StringForm["`` ``", curcol[[1]], Superscript[g,\[Mu] \[Nu]]];
            (* replace derivative by partial derivative *)
            curcol[[2]] = Map[HoldForm, curcol[[2]]] /. subhrzvvars;
            
            (* add mu and nu contraction indices *)
            curcol[[2,1]] = curcol[[2,1]] /. submu;
            For[ind2=2, ind2<=Length[curcol[[2]]], ind2++,
                isaftervector = Apply[Or, isvectorlist[[1;;(ind2-1)]]];
                If[isaftervector,
                    curcol[[2,ind2]] = curcol[[2,ind2]] /. subnu
                    ,
                    curcol[[2,ind2]] = curcol[[2,ind2]] /. submu
                ];
            ];
            
            rFRfourvertex[[ind]] = curcol;
            ,
            (* False goes here *)
            rFRfourvertex = Drop[rFRfourvertex, {ind}];
            ind = ind-1;
        ];
        ,
        (* False goes here *)
        False;
    ];
];


(* save variables to tempfile *)
Print["saving "<>tempfile];
RemoveTempFile[];
Save[tempfile, modelvars];
Save[tempfile, modelpars];
Save[tempfile, FRthreevertex];
Save[tempfile, FRfourvertex];
Save[tempfile, rFRthreevertex];
Save[tempfile, rFRfourvertex];
Save[tempfile, BosonSectorAssumptions];


]; (* end of Assuming[...] *)
]; (* end of If[FileExistsQ[...] *)
End[];
EndPackage[];

