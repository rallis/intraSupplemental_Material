(* ::Package:: *)

BeginPackage["GenDistanceFunctions`"]


SSIMmax::usage="SSIMmax[img1_,img2_,max_] computes SSIM distance"
mydescriptor::usage="mydescriptor[img_,winsizes_,steps_,histbins_] = Returns the descriptor of a  matrix/image, given the windowsizes, windowsteps and histogram bins"
ApproxInfoDist::usage="Computes ApproxInfoDist"


Begin["`Private`"]

SSIMmax[img1_,img2_,max_]:=Module[{w, c1, c2, m1, m2, m1sq, m2sq, m1m2, sigma1sq, sigma2sq,sigma12, ssimmap, mssim},
w = GaussianMatrix[{5, 1.5}]; c1 = (0.01*max)^2; c2 = (0.03*max)^2;
m1 = ListCorrelate[w, img1]; m2 = ListCorrelate[w, img2];
m1sq = m1*m1; m2sq = m2*m2; m1m2 = m1*m2;
sigma1sq = ListCorrelate[w, img1*img1] - m1sq;
sigma2sq = ListCorrelate[w, img2*img2] - m2sq;
sigma12 = ListCorrelate[w, img1*img2] - m1m2; 
ssimmap = ((c1 + 2*m1m2)*(c2 + 2*sigma12))/((c1 + m1sq + m2sq)*(c2 + sigma1sq + sigma2sq));
mssim = Mean[Mean[ssimmap]];
Return[mssim]];


mydescriptor[img_,winsizes_,steps_,histbins_]:=
Module[{dim=Length[img],allwinds={},winds,i,j,i1,nowsize,nowstep},
If[Length@winsizes!=Length@steps,Print["Error in winsizes and steps lists."];Return[]];
For[i1=1,i1<=Length[winsizes],i1++,
nowsize=winsizes[[i1]];nowstep=steps[[i1]];
winds=Table[Take[img,{i,nowsize+i-1},{j,nowsize+j-1}],{i,1,dim-nowsize+1,nowstep},{j,1,dim-nowsize+1,nowstep}];
AppendTo[allwinds,Flatten[Table[(Length/@BinLists[Sort[Flatten[winds[[i,j]]]],{histbins}])/(nowsize^2),{i,1,Length[winds]},{j,1,Length[winds]}]]];
];
Return[Flatten[allwinds]]
];

ApproxInfoDist[i_, j_] := 
Module[{x, y, xy, yx, cgr1, cgr2},
cgr1 = i; cgr2 = j;
x = Total[Unitize[cgr1], 2]; y = Total[Unitize[cgr2], 2];
xy = Total[Unitize[cgr1 + cgr2], 2];
yx = Total[Unitize[cgr1 + cgr2], 2];(*Print[{x,y,xy,yx}];*)
Return[N[(xy - y + yx - x)/xy]];
];

End[]

EndPackage[]




