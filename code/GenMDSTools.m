(* ::Package:: *)

BeginPackage["GenMDSTools`"]

getUpperHalf::usage="getUpperHalf[matr_List] accepts as input an nxn matrix and returns the upper half of it (excluding the diagonal), flattened"
mds::usage="mds[delta_List,dimensions_,dbg_] gets as input the distances delta_{ij}, and performs MDS"

Begin["`Private`"]

getUpperHalf[matr_List]:=
Module[{i,res={}},
For[i=1,i<=Length[matr]-1,i++,AppendTo[res,Take[matr[[i]],{i+1,Length[matr]}]]];Return[Flatten[res]]
];

mds[delta_List,dimensions_,dbg_,dim1_:1,dim2_:2,dim3_:3]:= (* DEBUG options= -1: PRINT NOTHING 0: all - 1: eigenvals+sol+MDSError *)
Module[{n,deltasq,deltatotals,sumOfDelta,bMatr,eigensys,diagMatr,solpts,dMatr,mdserrors,pointstress,nofmdserrors,pos,step,finalpos},
Monitor[
n=Length[delta];
deltasq=delta^2;
step=1;
deltatotals=Plus@@deltasq;
sumOfDelta=Plus@@deltatotals;
step=2;
bMatr=Table[-(1/2)*(deltasq[[i, j]]-(1/n)*deltatotals[[i]]-(1/n)*deltatotals[[j]]+(1/(n^2))*sumOfDelta), {i, 1, n}, {j, 1, n}];
(*jMatr=IdentityMatrix[Length[delta]]-(1/Length[delta])*Table[1,{k,1,Length[delta]},{l,1,Length[delta]}];*)
(*bMatr=-1/2*jMatr.deltasq.jMatr;*)
step=3;
eigensys=N[Eigensystem[N[bMatr,20],Min[10,n]],20];  (*eigensys=N[Eigensystem[N[bMatr,15]],15];*)
Print["10_First_Eigenvalues= ",eigensys[[1]]];
step=4;
pos=Select[Transpose[eigensys],#[[1]]>0&];
If[dimensions==2,finalpos={pos[[dim1]],pos[[dim2]]};];
If[dimensions==3,finalpos={pos[[dim1]],pos[[dim2]],pos[[dim3]]};];
pos=Transpose[finalpos]; (* backward compatibility *)
step=5;
diagMatr=Sqrt[DiagonalMatrix[pos[[1]]]];  (* "Dimensions" largest positive eigenvalues *)
step=6;
solpts=Transpose[pos[[2]]].diagMatr;
step=7;
dMatr={};
step=8;
mdserrors=0;
step=9;
nofmdserrors=0;
step=10;
pointstress={};

(* ALL dbg *)
If[dbg==0,Print["delta=",delta//MatrixForm]];
If[dbg==0,Print["bMatr=",N[bMatr,15]//MatrixForm]];
If[dbg==0,Print["eigensys=",eigensys]];
If[0<=dbg<=1,Print["eigenvals=",eigensys[[1]]]];
If[dbg==0,Print["diagMatr=",diagMatr]];
If[0<=dbg<=1,Print["sol=",solpts]];
If[dbg==0,Print["dMatr=",dMatr//MatrixForm]];
If[0<=dbg<=1,Print[mdserrors//MatrixForm]];
If[dbg==0,Print["pointStress=",pointstress//MatrixForm]];
If[0<=dbg<=1,Print[nofmdserrors//MatrixForm]];
(* RETURN *)
Return[{solpts,mdserrors,Select[eigensys[[1]],#>0&],pointstress,nofmdserrors}]
,{step,{i,j},{i1,i2}}]];



End[]

EndPackage[]
