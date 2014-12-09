(* ::Package:: *)

BeginPackage["GenChaosGame`"]


sparseFCGR::usage="sparseCGR[seq_,alphabet_List,k_] returns the FCGR plot of a genome as a sparse matrix. assumes linear, non-circular dna string"
FCGR::usage="FCGR[seq_,alphabet_List,k_] returns the grayscale Chaos Game Representation plot of a DNA sequence in a unit square in 2D using 2^k resolution.assumes linear, non-circular dna string"

Begin["`Private`"]

sparseFCGR[seq_,alphabet_List,k_]:=
Module[{chars,shifts,helpfunc1,helpfunc2,helpfunc3,pts,arraypts},
chars=StringCases[seq,alphabet[[1]]|alphabet[[2]]|alphabet[[3]]|alphabet[[4]]];
shifts=chars/.{alphabet[[1]]->{0.,0.},alphabet[[2]]->{0.,1.0*2^k},alphabet[[3]]->{1.0*2^k,1.0*2^k},alphabet[[4]]->{1.0*2^k,0.}};
helpfunc1=Compile[{{a,_Real,1}},FoldList[IntegerPart[(#+#2)/2]&,2^(k-1),a]];
helpfunc2=Compile[{{a,_Integer,1}},Map[2^k-#&,a]];
helpfunc3=Compile[{{a,_Integer,1}},Map[1+#&,a]];
pts=Transpose[helpfunc1/@Transpose[shifts]];
pts=Drop[pts,k];
arraypts=Transpose[{helpfunc2[Transpose[pts][[2]]],helpfunc3[Transpose[pts][[1]]]}];
Return[{Length[pts],Tally[arraypts]}];
];



FCGR[seq_,alphabet_List,k_]:=
Module[{resimg,chars,shifts,helpfunc1,helpfunc2,helpfunc3,pts,arraypts,setgrayscale},
setgrayscale[{a_,b_}]:=Module[{oldval,newval},oldval=resimg[[a,b]];newval=oldval+1;resimg[[a,b]]=newval;];
resimg=Table[0,{i,1,2^k},{j,1,2^k}];
chars=StringCases[seq,alphabet[[1]]|alphabet[[2]]|alphabet[[3]]|alphabet[[4]]];
shifts=chars/.{alphabet[[1]]->{0.,0.},alphabet[[2]]->{0.,1.0*2^k},alphabet[[3]]->{1.0*2^k,1.0*2^k},alphabet[[4]]->{1.0*2^k,0.}};
helpfunc1=Compile[{{a,_Real,1}},FoldList[IntegerPart[(#+#2)/2]&,2^(k-1),a]];
helpfunc2=Compile[{{a,_Integer,1}},Map[2^k-#&,a]];
helpfunc3=Compile[{{a,_Integer,1}},Map[1+#&,a]];
pts=Transpose[helpfunc1/@Transpose[shifts]];
pts=Drop[pts,k];
arraypts=Transpose[{helpfunc2[Transpose[pts][[2]]],helpfunc3[Transpose[pts][[1]]]}];
Map[setgrayscale[{#[[1]],#[[2]]}]&,arraypts];
Return[{Length[pts],resimg}];
];


End[]

EndPackage[]




