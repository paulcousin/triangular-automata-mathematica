(* ::Package:: *)

(* ::Title:: *)
(*Triangular Automata*)


(* ::Abstract:: *)
(*This package is made to compute and display cellular automata in an infinite triangular grid.*)


(* ::Author:: *)
(*Paul Cousin*)
(*https://orcid.org/0000-0002-3866-7615*)


(* ::Section:: *)
(*Begin*)


BeginPackage["TriangularAutomata`"];


(* ::Subsection:: *)
(*Utilities*)


TAConfigurationVector::usage="TAConfigurationVector[grid] returns the configuration vector of the grid.";
TANegativeGrid::usage="TANegativeGrid[grid] returns the grid with all states inverted";
TANegativeRule::usage="TANegativeRule[ruleNumber] returns the rule number that would have the same effect in a negative grid.";
TALayer::usage="TALayer[grid] returns the number of layers in the grid.";
TAPopulation::usage="TAPopulation[grid] returns the population the grid.";


(* ::Subsection:: *)
(*Grid*)


TAGrid::usage="TAGrid[stateVector] is an object representing a grid.";
TAStateVector::usage="TAStateVector[TAGrid] returns the state vector.";
TAAdjacencyMatrix::usage="TAAdjacencyMatrix[l] gives the adjacency matrix of the l layers wide grid.";
TACoordinates::usage="TACoordinates[l] gives the coordinate vector of the l layers wide grid.";


(* ::Subsection:: *)
(*Evolution*)


TAEvolve::usage="TAEvolve[grid,ruleNumber] evolves the grid with the given rule.";
TANestEvolve::usage="TANestEvolve[grid,ruleNumber,n] evolves the grid n times with the given rule and returns the final grid.";
TANestListEvolve::usage="TANestListEvolve[grid,ruleNumber,n] evolves the grid n times with the given rule and returns the list of all intermediate grids.";


(* ::Subsection:: *)
(*Plots*)


TAConfigurationPlot::usage="TAConfigurationPlot[c] plots configuration number c.";
TATransformationPlot::usage="TATransformationPlot[c,s] illustrates the transformation from configuration number c to state s.";
TARulePlot::usage="TARulePlot[ruleNumber] makes a plot illustrating the given rule.";
TAGridPlot::usage="TAGridPlot[grid] displays the grid.";
TAGridPlot3D::usage="TAGridPlot3D[grid] displays the grid in 3D.";
TAEvolutionPlot::usage="TAEvolutionPlot[grid,ruleNumber,numberOfSteps] generates an animation of the evolution of the grid.";
TAEvolutionPlot3D::usage="TAEvolutionPlot[grid,ruleNumber,numberOfSteps] generate a 3D representation of the grid evolution.";


(* ::Subsection:: *)
(*Starting Points*)


TAStartOneAlive::usage="TAStartOneAlive is a grid with only one alive cell.";
TAStartLogo::usage="TAStartLogo is a grid with alive cells in the shape of the letters TA that includes all local configurations.";
TAStartRandom::usage="TAStartRandom[n] is a grid with a random distribution of alive and dead cells on n layers.";


(* ::Subsection:: *)
(*Torus*)


TATorus::usage="TATorus[stateVector] is an object representing a torus.";
TAStateVectorT::usage="TAStateVectorT[TAGrid] returns the state vector.";
TAAdjacencyMatrixT::usage="TAAdjacencyMatrixT[l] gives the adjacency matrix of the l layers wide torus.";
TACoordinatesT::usage="TACoordinatesT[l] gives the coordinate vector of the l layers wide torus.";
TATorusPlot::usage="TATorusPlot[torus] displays the torus.";
TAStartRandomTorus::usage="TAStartRandomTorus[n] is a torus with a random distribution of alive and dead cells on n layers.";


(* ::Subsection:: *)
(*Miscellaneous *)


TAEdit::usage="TAEdit[grid] allows you to edit the grid and copy the result.";
TACenterColumn::usage="TACenterColumn[ruleNumber,n] gives the center column of the grid evolved n times.";
TAThreeCenterColumns::usage="TAThirdCenterColumn[ruleNumber,n] gives the sequence of states of the 0th, 1st and 2nd layers.";
TASlice::usage="TASlice[ruleNumber,n] returns a slice of space-time plot.";
TASlicePlot::usage="TASlicePlot[lines] returns plot of the slice.";


(* ::Section:: *)
(*Functions*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
aliveColor=RGBColor[0.5, 0, 0.5];deadColor=GrayLevel[1];unknownColor = GrayLevel[0.85];
exactCoordinates= False;
onlineDirectory="https://files.paulcousin.net/triangular-automata";
localDirectory=FileNameJoin@{Directory[],"Triangular Automata"};


(* ::Subsection:: *)
(*Utilities*)


(* ::Input::Initialization:: *)
binary[n_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(23\)]\(\(IntegerDigits[n, 10, 24]\)[\([\(-1\) - i]\)]*
\*SuperscriptBox[\(2\), \(i\)]\)\);


layerFromOrder[order_]:=Ceiling[1/6 (-3+Sqrt[3(8 order-5)])];

layerFromMatrix[matrix_] := Module[{order = Max @ Dimensions @ matrix, n, t},
 Return@layerFromOrder[order]];

layerFromCoords[coords_] := Module[{order = Length @ coords, n, t},
Return@layerFromOrder[order]];
	
orderFromLayer[layer_]:=1+3 (layer(layer+1))/2;

TALayer[TAGrid[sV_]]:=layerFromOrder@Length@sV;

TAPopulation[grid_]:=Total[If[Last@#[[1]]=={1},1-#[[1]],#[[1]]]][[1]]&@grid;


(* ::Input::Initialization:: *)
TAConfigurationVector[TAGrid[sV_]]:=Normal[4*sV+TAAdjacencyMatrix[layerFromOrder@Length@sV] . sV][[All,1]];


TANegativeGrid[TAGrid[sV_]]:=TAGrid[SparseArray@Normal@(1-sV)];

TANegativeRule[ruleNumber_]:=FromDigits[1-Reverse@IntegerDigits[ruleNumber,2,8],2];


(* ::Subsection:: *)
(*Grid*)


 TAGrid/:MakeBoxes[TAGrid[sV_],form_]:=With[
 {
 interpretation=Interpretation[
 If[layerFromOrder[Length[sV]]-2<=64,
 Deploy@TAGridPlot[TAGrid[sV],"ImageSize"->Tiny,"Frame"->(sV[[-1,-1]]==0)],
 Deploy@largeGridRepresentation[layerFromOrder[Length[sV]]-2,Total[Abs[sV[[-1,-1]]-sV][[All,1]]],If[sV[[-1,-1]]==0,White,Purple]]
 ],TAGrid[sV]]
 },
 MakeBoxes[interpretation,form]];


largeGridRepresentation[layers_,population_,universe_]:=Panel[Grid[{{
Spacer[0],
Graphics[{GrayLevel[0.8], Translate[Triangle[{{-3^Rational[-1, 2], 0}, {Rational[1, 2] 3^Rational[-1, 2], Rational[1, 2]}, {Rational[1, 2] 3^Rational[-1, 2], Rational[-1, 2]}}], {{{0.8660254037844388, -0.5}}, {{0.8660254037844388, 0.5}}, {{-0.8660254037844384, 0.5}}, {{1.7320508075688774`, 0.}}, {{1.7320508075688774`, 1.}}, {{0.8660254037844388, 1.5}}, {{-1.732050807568877, 0.}}, {{-0.8660254037844384, -1.5}}, {{0.8660254037844388, -1.5}}, {{2.5980762113533165`, 0.5}}, {{-2.598076211353315, -0.5}}, {{3.4641016151377553`, 0.}}}], Translate[Triangle[{{3^Rational[-1, 2], 0}, {Rational[-1, 2] 3^Rational[-1, 2], Rational[1, 2]}, {Rational[-1, 2] 3^Rational[-1, 2], Rational[-1, 2]}}], {{{0.5773502691896258, 0.}}, {{-0.28867513459481275`, 0.5}}, {{1.4433756729740645`, -0.5}}, {{0.577350269189626, 1.}}, {{-1.1547005383792512`, 0.}}, {{-1.1547005383792512`, -1.}}, {{0.577350269189626, -1.}}, {{2.3094010767585034`, 0.}}, {{2.3094010767585034`, 1.}}, {{1.4433756729740648`, 1.5}}, {{-2.02072594216369, -0.5}}, {{-0.28867513459481264`, -2.5}}, {{3.175426480542942, 0.5}}, {{-2.8867513459481273`, -1.}}}]}, Background -> GrayLevel[0, 0], ImageSize -> 50],
Spacer[0],
Grid[{{Text[Style["layers:",Gray,Medium]],Text[Style[ToString@layers,Medium]]},
{Text[Style["population:",Gray,Medium]],Text[Style[ToString@population,Medium]]},
{Text[Style["universe:",Gray,Medium]],universe}},
Alignment->Left,Spacings->{.6,.2}],Spacer[0]
}}],FrameMargins->4];


(* ::Subsubsection:: *)
(*State Vector & Layer*)


TAStateVector[grid_]:=grid[[1]];


(* ::Subsubsection::Closed:: *)
(*Adjacency Matrix*)


(* ::Subsubsubsection::Closed:: *)
(*Useful functions*)


(* ::Input::Initialization:: *)
cornerMerge[list_]:=Module[{i,dims, output=list[[1]]},

For[i=2,i<=Length@list,i++,
dims=Dimensions[output]+Dimensions[list[[i]]];
output=PadRight[output,dims]+PadLeft[list[[i]],dims]
];

Return@output
];

triple[m_]:=cornerMerge[{m,m,m}];

shift[m_]:=m[[RotateRight@Range@Dimensions[m][[1]],RotateRight@Range@Dimensions[m][[2]]]];

stairs[n_]:=Transpose@SparseArray[{Band[{2,2}]->1,Band[{1,2}]->1},{n+1,n+1}][[1;;n+1,2;;n+1]];

symmetrize[m_]:=Module[{dim=Dimensions[m][[2]],output},
output=PadRight[m,{dim,dim}];
 Return[output+Transpose[output]];
];


(* ::Subsubsubsection:: *)
(*Computing the Adjacency Matrix*)


(* ::Input::Initialization:: *)
expandMatrix[matrix_]:=Module[{newLayer=layerFromMatrix@matrix+1},
Which[
newLayer==0,Return@SparseArray[Automatic, {1, 1}, 0, {1, {{0, 0}, {}}, {}}],
newLayer==1, Return@SparseArray[Automatic, {1, 4}, 0, {1, {{0, 3}, {{2}, {3}, {4}}}, {1, 1, 1}}],
newLayer==2,Return@SparseArray[Automatic, {4, 10}, 0, {1, {{0, 3, 5, 7, 9}, {{2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}}}, {1, 1, 1, 1, 1, 1, 1, 1, 1}}],
newLayer==3,Return@SparseArray[Automatic, {10, 19}, 0, {1, {{0, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21}, {{2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {19}, {12}, {13}, {13}, {14}, {15}, {16}, {16}, {17}, {18}, {19}}}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}],
EvenQ[newLayer],Return@cornerMerge[{matrix,triple@cornerMerge@{stairs[newLayer/2],SparseArray@IdentityMatrix[newLayer/2-1]}}],
OddQ[newLayer],Return@cornerMerge[{matrix,shift@triple@cornerMerge@{SparseArray@IdentityMatrix[Ceiling[newLayer/2]-2],stairs[Ceiling[newLayer/2]]}}]
]]


gridMatrix=SparseArray[Automatic, {10, 19}, 0, {1, {{0, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21}, {{2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {19}, {12}, {13}, {13}, {14}, {15}, {16}, {16}, {17}, {18}, {19}}}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}];


TAAdjacencyMatrix[layer_]:=Module[{
currentLayer=layerFromOrder@Dimensions[gridMatrix][[2]]
},
If[layer>currentLayer+512,
	Module[{fileName,file,remoteFileName,toDownloadLayer=Min[2^Ceiling[Log[2,layer]],8192]},
	fileName="TAMatrix-"<>ToString@toDownloadLayer<>".wl";
	file=If[
		FileExistsQ[FileNameJoin[{localDirectory,fileName}]],
		PrintTemporary["Found matrix in local directory."];
		FileNameJoin[{localDirectory,fileName}],
		PrintTemporary["Downloading matrix..."];
		ExtractArchive[onlineDirectory<>"/matrix/"<>fileName<>".zip",Directory[]][[1]]
	];
	gridMatrix=Import@file;
	];
	currentLayer=layerFromOrder@Dimensions[gridMatrix][[2]];
];
If[currentLayer<layer,gridMatrix=Nest[expandMatrix,gridMatrix,layer-currentLayer]];
Return@symmetrize@gridMatrix[[All,;;orderFromLayer@layer]]
]


(* ::Subsubsection:: *)
(*Vertex Coordinates*)


(* ::Input::Initialization:: *)
TAExpandCoords[c_] := Module[
{i,coords = c, layer = layerFromCoords @ c +1, new,floor, ceiling},

floor=Floor[layer / 2];
ceiling =Ceiling[layer / 2];

new={coords[[Min[-3 * (layer - 1), -1]]] + If[OddQ[layer],
			{1 / Sqrt[3], 0},
			{1 / (2 Sqrt[3]), -(1/2)}
		]};

	Table[
		AppendTo[
			new
			,
			Last[new] + Which[
				i < floor,
					{0, 1}
				,
				i < floor + ceiling,
					{-(Sqrt[3] / 2), 1/2}
				,
				i < 2 floor + ceiling,
					{-(Sqrt[3] / 2), -(1/2)}
				,
				i < 2 floor + 2 ceiling,
					{0, -1}
				,
				i < 3 floor + 2 ceiling,
					{Sqrt[3] / 2, -(1/2)}
				,
				i < 3 floor + 3 ceiling,
					{Sqrt[3] / 2, 1/2}
			]
		]
		,
		{i, 0, 3 * layer - 2}
	];
	Return @Join[coords,new];
]


(* ::Input::Initialization:: *)
buildCoords[l_]:=Module[{i,coords=If[exactCoordinates,{{0,0}},{{0.,0.}}]},
Monitor[For[i=0,i<l,i++,coords=TAExpandCoords@coords],i];Return@coords];


coordinates=If[exactCoordinates,{{0,0}},{{0.,0.}}];


TACoordinates[layer_]:=Module[{
currentLayer=layerFromOrder@Length@coordinates
},
If[layer>currentLayer+512,
	Module[{fileName,file,toDownloadLayer=Min[2^Ceiling[Log[2,layer]],6000]},
	fileName=If[exactCoordinates,"TAExactCoords-","TACoords-"]<>ToString@toDownloadLayer<>".wl";
	file=If[
		FileExistsQ[FileNameJoin[{localDirectory,fileName}]],
		PrintTemporary["Found coordinates in local directory."];
		FileNameJoin[{localDirectory,fileName}],
		PrintTemporary["Downloading coordinates..."];
		ExtractArchive[onlineDirectory<>"/coordinates/"<>fileName<>".zip",Directory[]][[1]]
	];
	coordinates=Import@file;
	];
	currentLayer=layerFromOrder@Length@coordinates;
];
If[(layer>currentLayer),coordinates=Nest[TAExpandCoords,coordinates,layer-currentLayer]];
Return@coordinates[[;;orderFromLayer@layer]];
]


(* ::Subsection:: *)
(*Evolution*)


(* ::Input::Initialization:: *)
ruleState[ruleNumber_,case_]:=IntegerDigits[ruleNumber,2,8][[-case-1]];


(* ::Input::Initialization:: *)
configurationVector[aM_,sV_]:=Normal[4*sV+aM . sV];


(* ::Input::Initialization:: *)
newState[aM_,sV_,ruleNumber_]:= SparseArray[Map[ruleState[ruleNumber,#]&,configurationVector[aM,sV]]];


(* ::Input::Initialization:: *)
Options[TAEvolve]={"Shrink"->False};

TAEvolve[gr_,rN_,OptionsPattern[]]:=Module[
{grid=gr,stateVector,layer,matrix,vertexCoords,digits,lastSpecified},

digits=IntegerDigits[rN,2,8];
stateVector=TAStateVector[grid];
layer=TALayer@grid;

If[Head@grid===TAGrid,
	lastSpecified=If[#=={},0,layerFromOrder@Last[#]]&@
	If[stateVector[[-1,-1]]==0,
	Position[Normal@stateVector,1][[All,1]],
	Position[Normal@stateVector,0][[All,1]]
	];

	If[(lastSpecified==layer-2 )\[And] !OptionValue["Shrink"] ,
		layer=(layer+1);
		stateVector=ArrayReshape[stateVector,{orderFromLayer[layer],1},stateVector[[-1,1]]];
	];
	
	If[ OptionValue["Shrink"],
		layer=(layer-1);
		stateVector=stateVector[[;;orderFromLayer[layer],All]];stateVector=ReplacePart[stateVector,Table[{n,1}->Normal@stateVector[[-1-3layer,1]],{n,-3layer,-1}]];
	];
];

matrix=Which[
Head@grid===TAGrid,TAAdjacencyMatrix[layer],
Head@grid===TATorus,TAAdjacencyMatrixT[layerFromOrderT@Length@grid[[1]]]
];

stateVector=newState[matrix,stateVector,rN];

If[Head@grid===TAGrid,
stateVector=ReplacePart[stateVector,Table[{n,1}->Normal@stateVector[[-1-3layer,1]],{n,-3layer,-1}]];
];

Return@Which[
Head@grid===TAGrid,TAGrid[stateVector],
Head@grid===TATorus,TATorus[stateVector]
];
];


TANestEvolve[gr_,rN_,n_]:=Module[{i,grid=gr,rule=If[ListQ@rN,rN,{rN}]},
Monitor[For[i=1,i<=n,i++,grid=TAEvolve[grid,rule[[Mod[i,Length@rule,1]]]]],Grid@{{ProgressIndicator[Dynamic[i],{0,n}],"computing grid "<>ToString[i]<>" of "<>ToString[n]}}];
Return@grid;
];

TANestListEvolve[gr_,rN_,n_]:=Module[{i,grids={gr},rule=If[ListQ@rN,rN,{rN}]},
Monitor[For[i=1,i<=n,i++,AppendTo[grids,TAEvolve[Last@grids,rule[[Mod[i,Length@rule,1]]]]]],Grid@{{ProgressIndicator[Dynamic[i],{0,n}],"computing grid "<>ToString[i]<>" of "<>ToString[n]}}];
Return@grids;
];


(* ::Subsection:: *)
(*Plots*)


TAConfigurationPlot[c_]:=
Graphics[{
If[c<4,White,Purple], If[c<4,EdgeForm[Thickness@Small],EdgeForm[None]],
	Triangle[{{-(1/Sqrt[3]),0},{1/(2 Sqrt[3]),1/2},{1/(2 Sqrt[3]),-(1/2)}}], 
If[Mod[c,4]>0,Purple,White],If[Mod[c,4]>0,EdgeForm[None],EdgeForm[Thickness@Small]],
	Translate[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3],0}}],
If[Mod[c,4]>1,Purple,White],If[Mod[c,4]>1,EdgeForm[None],EdgeForm[Thickness@Small]],
	Translate[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3]-Sqrt[3]/2,1/2}}],
If[Mod[c,4]>2,Purple,White],If[Mod[c,4]>2,EdgeForm[None],EdgeForm[Thickness@Small]],
	Translate[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}], 1.1{{1/Sqrt[3]-Sqrt[3]/2,-(1/2)}}]
},ImageSize->{45,45}];


(* ::Input::Initialization:: *)
TATransformationPlot[before_,after_]:=
Graphics[{
Inset[TAConfigurationPlot[before],{-1.3,0},Automatic,Scaled@1],Text[Style["\[Rule]",FontSize->Scaled@.25],Scaled@{.5,.528}],Inset[
Graphics[{
If[after==0,White,Purple], If[after==0,EdgeForm[Thickness@Small],EdgeForm[None]],
	Triangle[{{-(1/Sqrt[3]),0},{1/(2 Sqrt[3]),1/2},{1/(2 Sqrt[3]),-(1/2)}}],
unknownColor,EdgeForm[None],
	Translate[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3],0}}],
	Translate[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3]-Sqrt[3]/2,1/2}}],
	Translate[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}], 1.1{{1/Sqrt[3]-Sqrt[3]/2,-(1/2)}}]
}],{1.3,0},Automatic,Scaled@1]
}];


Options[TARulePlot]={"Labeled"->False,"Portrait"->False,"Frame"->True,"ImageSize"->Automatic};

TARulePlot[rN_,OptionsPattern[]] :=Module[{ruleNumber=rN,graphicsGrid},
graphicsGrid=Graphics[
Inset[
	TATransformationPlot[#,ruleState[ruleNumber,#]],
	If[OptionValue["Portrait"]==False,{-2*Mod[#,4],.5(#-Mod[#,4])},{.5(Mod[#,4]-#),2*Mod[#,4]}],
	Automatic,
	If[OptionValue["Portrait"]==False,Scaled@.32,Scaled@.45]
]&/@Range[7,0,-1], 
Frame->OptionValue["Frame"], FrameTicks->None, FrameStyle->Thick,
Background->White,
ImagePadding->None, ImageMargins->None,
AspectRatio->If[OptionValue["Portrait"]==False,1/3.5,1.2],
PlotRange->All,
PlotRangePadding ->{1.1,1},
ImageSize->Which[
OptionValue["ImageSize"]=!=Automatic,OptionValue["ImageSize"],
OptionValue["Portrait"]==False,500,
OptionValue["Portrait"]==True,250]
];
If[OptionValue["Labeled"],
	Return@Grid[{{Text[ruleNumber " = " BaseForm[ruleNumber, 2]]},
	{graphicsGrid}}],
	Return@graphicsGrid
];
];


(* ::Input::Initialization:: *)
Options[TAGridPlot]={
"ImageSize" -> Medium,
"Time"->Null,
 "Padding"->Scaled[.03],
"Parallel"->False,
"Frame"->False,
"FrameStyle" ->Automatic,
"Rotation"->0,
"Translation"->{0,0},
"PlotRegion"->Automatic,
"AspectRatio"->Automatic,
"PlotRange"->Automatic,
"PlotRangePadding"->Automatic
};

TAGridPlot[grid_,OptionsPattern[]]:=Module[
{x,head,stateVector,displayedVertices,triangleEven, triangleOdd, trianglesEven={},trianglesOdd={},gridstate,color,graphicsList={},evenLayered,oddLayered,coords,imageSize=OptionValue["ImageSize" ]},

stateVector=grid[[1]];
gridstate=stateVector[[-1,-1]];
coords=TACoordinates@layerFromOrder@Length@stateVector;
coords=coords . RotationMatrix[-OptionValue["Rotation" ]];
coords=(#+OptionValue["Translation" ])&/@coords;

(*basic shapes*)
triangleEven=Triangle[{{-1/Sqrt[3],0},{1/(2Sqrt[3]),1/2},{1/(2Sqrt[3]),-1/2}} . RotationMatrix[-OptionValue["Rotation" ]]];
triangleOdd=Triangle[{{1/Sqrt[3],0},{-1/(2Sqrt[3]),1/2},{-(1/(2Sqrt[3])),-1/2}} . RotationMatrix[-OptionValue["Rotation" ]]];

(*select background color*)
If[gridstate==0,
color=aliveColor;
displayedVertices=Drop[ArrayRules@stateVector,-1][[;;,1,1]];
,
color= deadColor;
displayedVertices=Drop[ArrayRules[ConstantArray[1,Dimensions@stateVector]-stateVector],-1][[;;,1,1]];
];


(*compute two sets of coordinates for even and odd layers*)
evenLayered=Map[EvenQ@layerFromOrder@#&,displayedVertices];
oddLayered=Not/@evenLayered;

evenLayered = SparseArray[#->True&/@Pick[displayedVertices,evenLayered],Length@coords,False];
oddLayered= SparseArray[#->True&/@Pick[displayedVertices,oddLayered],Length@coords,False];

trianglesEven=Pick[coords,evenLayered];
trianglesOdd=Pick[coords,oddLayered];

(*construct Graphics object*)
AppendTo[graphicsList,color];
If[trianglesEven!={},AppendTo[graphicsList,Translate[triangleEven,trianglesEven]]];
If[trianglesOdd!={},AppendTo[graphicsList,Translate[triangleOdd,trianglesOdd]]];
If[OptionValue["Time"]=!=Null,AppendTo[graphicsList,Text[Style[OptionValue["Time"],FontSize->Scaled@.06],Scaled[{.98,.01}],{Right,Bottom}]]];

(*output with different Padding options*)
If[OptionValue["ImageSize" ]=="Proportional",
imageSize=128*Norm@Last[coords];
];

Return@Graphics[
graphicsList,
Background->If[gridstate==0,deadColor,aliveColor],
PlotRange->Which[
OptionValue["Padding"]===None,Automatic,
OptionValue["PlotRange"]===Automatic,Norm@coords[[-1]],
True,OptionValue["PlotRange"]],
ImageSize->imageSize,
Frame->OptionValue["Frame" ],
FrameTicks->False,
FrameStyle->OptionValue["FrameStyle" ],
PlotRegion->OptionValue["PlotRegion" ],
PlotRangePadding->OptionValue["Padding"],
AspectRatio->OptionValue["AspectRatio" ]
];

];


Options[TAEvolutionPlot]={
"ImageSize" -> Medium, 
"Timed"->True, 
"Pause"->False, 
"Animated"->True, 
"Select"->"All"};

TAEvolutionPlot[grid_, ruleNumber_, steps_, OptionsPattern[]] := 
Module[{x,i, grids = TANestListEvolve[grid, ruleNumber, steps]},
	
	grids = Transpose@{grids,Range[0,Length@grids-1]};
	If[OptionValue["Select"]=="Even", grids=grids[[1;;;;2]]];
	If[OptionValue["Select"]=="Odd", grids=grids[[2;;;;2]]];
	grids = Which[
	Head[#[[1]]]===TAGrid,
	TAGridPlot[#[[1]], "ImageSize" -> OptionValue["ImageSize"],"Time"->If[OptionValue["Timed"],#[[2]],Null]],
	Head[#[[1]]]===TATorus,
	TATorusPlot[#[[1]], "ImageSize" -> OptionValue["ImageSize"],"Time"->If[OptionValue["Timed"],#[[2]],Null]]
	]& /@ grids;
	
	
	If[OptionValue["Pause"],
	For[i=0,i<Length[grids]/16,i++,
		PrependTo[grids,First@grids];
		AppendTo[grids,Last@grids];
		AppendTo[grids,Last@grids];
	]];	
	
	Return @ If[OptionValue["Animated"],ListAnimate[grids],grids];
]


(* ::Input::Initialization:: *)
Options[TAGridPlot3D]={"Translation"->{0,0}};

TAGridPlot3D[TAGrid[sV_],time_:0, OptionsPattern[]]:=Module[{x,level=time,stateVector=sV,coords,displayedVertices,shapeEven,shapeOdd,shapesEven={},shapesOdd={},shapes={},evenLayered,oddLayered,gridstate=sV[[-1,-1]],color},

coords=TACoordinates@layerFromOrder@Length@stateVector;
coords=(#+OptionValue["Translation" ])&/@coords;

(*basic shapes*)
shapeEven=MeshRegion[{{-1/Sqrt[3],0,0},{1/(2Sqrt[3]),1/2,0},{1/(2Sqrt[3]),-(1/2),0},{-1/Sqrt[3],0,1},{1/(2Sqrt[3]),1/2,1},{1/(2Sqrt[3]),-(1/2),1}},{Polygon[{{1,2,3},{4,5,6}}],Polygon[{{1,2,5,4},{2,3,6,5},{3,1,4,6}}]}];
shapeOdd=MeshRegion[{{1/Sqrt[3],0,0},{-(1/(2Sqrt[3])),1/2,0},{-(1/(2Sqrt[3])),-(1/2),0},{1/Sqrt[3],0,1},{-(1/(2Sqrt[3])),1/2,1},{-(1/(2Sqrt[3])),-(1/2),1}},{Polygon[{{1,2,3},{4,5,6}}],Polygon[{{1,2,5,4},{2,3,6,5},{3,1,4,6}}]}];


(*select background color*)
If[gridstate==0,
color=aliveColor;
displayedVertices=Drop[ArrayRules@stateVector,-1][[;;,1,1]];
,
color= deadColor;
displayedVertices=Drop[ArrayRules[ConstantArray[1,Dimensions@stateVector]-stateVector],-1][[;;,1,1]];
];

(*compute two sets of coordinates for even and odd layers*)
evenLayered=Map[EvenQ@layerFromOrder@#&,displayedVertices];
oddLayered=Not/@evenLayered;

evenLayered = SparseArray[#->True&/@Pick[displayedVertices,evenLayered],Length@coords,False];
oddLayered= SparseArray[#->True&/@Pick[displayedVertices,oddLayered],Length@coords,False];

shapesEven=Pick[coords,evenLayered];
shapesOdd=Pick[coords,oddLayered];

shapesEven=Transpose@Append[Transpose[shapesEven],ConstantArray[-level,Length@shapesEven]];
shapesOdd=Transpose@Append[Transpose[shapesOdd],ConstantArray[-level,Length@shapesOdd]];

(*output*)
Return@If[Union[shapesEven,shapesOdd]=={},EmptyRegion[3],RegionUnion[Map[Translate[shapeEven,#]&,shapesEven]\[Union]Map[Translate[shapeOdd,#]&,shapesOdd]]
]
];


Options[TAEvolutionPlot3D]={
"Mesh" -> False,
"ImageSize"->Medium,
"ViewPoint"->{1,0,.75},
"ViewProjection"->"Orthographic",
"Translation"->{0,0}
};

TAEvolutionPlot3D[grid_,ruleNumber_,steps_, OptionsPattern[]]:=Module[
	{regions,gridPlots},
	
	regions= TANestListEvolve[grid,ruleNumber,steps];
	gridPlots=TAGridPlot3D[regions[[#]],#,"Translation"->{0,0}]&/@Range@Length@regions;

	If[OptionValue["Mesh"],
	Return@RegionUnion@gridPlots,
	Return@Graphics3D[RegionUnion@gridPlots,
		Boxed->False,
		Method->{"ShrinkWrap" -> True},
		ViewProjection->OptionValue["ViewProjection"],
		ViewPoint->OptionValue["ViewPoint"],
		ViewVertical->{0,0,1},
		ImageSize->OptionValue["ImageSize"]
]]];


(* ::Subsection:: *)
(*Starting Points*)


(* ::Input::Initialization:: *)
buildGrid[n_] :=
	Module[{matrix = buildMatrix[n]},
		Return[{matrix, SparseArray[{1, 1} -> 1, {Max @ Dimensions @ matrix,
			 1}], buildCoords[n]}]
	]

TAStartOneAlive =TAGrid[SparseArray[Automatic, {10, 1}, 0, {1, {{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {{1}}}, {1}}]];

TAStartRandom[n_,d_:.5] :=TAGrid[ArrayReshape[SparseArray@Transpose@{Table[Round[RandomReal[{-.5,.5}]+d],orderFromLayer[n]]},{orderFromLayer[n+2],1}]];

TAStartLogo=TAGrid[
SparseArray[({#,1}->1)&/@{2,3,5,6,8,11,13,16,17,19,21,22,23,27,29,31,33,34,35,41,44,49,58,67,78,88},{166,1}]];


(* ::Subsection:: *)
(*Torus*)


layerSizeT[l_]:=Max[0,12l-6];orderFromLayerT[l_]:=6*l^2;

layerFromOrderT[order_]:=Ceiling[Sqrt[order/6]];

layerFromMatrixT[matrix_] := Module[{order = Max @ Dimensions @ matrix, n, t},
 Return@layerFromOrderT[order]];

layerFromCoordsT[coords_] := Module[{order = Length @ coords, n, t},
Return@layerFromOrderT[order]];
	



 TATorus/:MakeBoxes[TATorus[sV_],form_]:=With[
 {
 interpretation=Interpretation[
 If[layerFromOrderT[Length[sV]]<=32,
 Deploy@TATorusPlot[TATorus[sV],"ImageSize"->Tiny,"Frame"->True],
 Deploy@largeTorusRepresentation[layerFromOrderT[Length[sV]],Total[sV][[1]]]
 ],TATorus[sV]]
 },
 MakeBoxes[interpretation,form]];


largeTorusRepresentation[layers_,population_]:=Panel[Grid[{{
Spacer[0],
Graphics[{{EdgeForm[{Thickness[0.02`],GrayLevel[0.8]}],GrayLevel[1, 0],InterpretationBox[PolygonBox[{{2.780837256811621`*^-16,4.541451884327381`},{-3.9330127018922187`,2.2707259421636916`},{-3.93301270189222`,-2.2707259421636894`},{-8.342511770434863`*^-16,-4.541451884327381`},{3.9330127018922183`,-2.2707259421636925`},{3.9330127018922214`,2.270725942163687`}}],RegularPolygon[{4.541451884327381`,\[Pi]/2},6]]},GrayLevel[0.8],GeometricTransformationBox[TagBox[PolygonBox[NCache[{{-(1/Sqrt[3]),0},{1/(2 Sqrt[3]),1/2},{1/(2 Sqrt[3]),-(1/2)}},{{-0.5773502691896258`,0},{0.2886751345948129`,0.5`},{0.2886751345948129`,-0.5`}}]],"Triangle"],NCache[{{{1/(4 Sqrt[3])+Sqrt[3]/4,0}},{{7/(4 Sqrt[3])-Sqrt[3]/4,-1}},{{13/(4 Sqrt[3])-Sqrt[3]/4,1/2}},{{7/(4 Sqrt[3])-Sqrt[3]/4,1}},{{-(5/(4 Sqrt[3]))-Sqrt[3]/4,1}},{{-(5/(4 Sqrt[3]))-Sqrt[3]/4,-1}},{{19/(4 Sqrt[3])-Sqrt[3]/4,1}},{{13/(4 Sqrt[3])-Sqrt[3]/4,3/2}},{{7/(4 Sqrt[3])-Sqrt[3]/4,2}},{{-(11/(4 Sqrt[3]))-Sqrt[3]/4,1/2}},{{25/(4 Sqrt[3])-Sqrt[3]/4,1/2}},{{-(17/(4 Sqrt[3]))-Sqrt[3]/4,0}}},{{{0.5773502691896257`,0}},{{0.5773502691896258`,-1}},{{1.4433756729740645`,0.5`}},{{0.5773502691896258`,1}},{{-1.1547005383792515`,1}},{{-1.1547005383792515`,-1}},{{2.3094010767585034`,1}},{{1.4433756729740645`,1.5`}},{{0.5773502691896258`,2}},{{-2.0207259421636903`,0.5`}},{{3.175426480542942`,0.5`}},{{-2.886751345948129`,0}}}]],GeometricTransformationBox[TagBox[PolygonBox[NCache[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}},{{0.5773502691896258`,0},{-0.2886751345948129`,0.5`},{-0.2886751345948129`,-0.5`}}]],"Triangle"],NCache[{{{1/(2 Sqrt[3]),-(1/2)}},{{-(1/(4 Sqrt[3]))+Sqrt[3]/4,1/2}},{{11/(4 Sqrt[3])-Sqrt[3]/4,0}},{{5/(4 Sqrt[3])-Sqrt[3]/4,3/2}},{{-(1/(4 Sqrt[3]))-Sqrt[3]/4,1}},{{-(7/(4 Sqrt[3]))-Sqrt[3]/4,1/2}},{{-(7/(4 Sqrt[3]))-Sqrt[3]/4,-(1/2)}},{{17/(4 Sqrt[3])-Sqrt[3]/4,1/2}},{{17/(4 Sqrt[3])-Sqrt[3]/4,3/2}},{{11/(4 Sqrt[3])-Sqrt[3]/4,2}},{{-(13/(4 Sqrt[3]))-Sqrt[3]/4,0}},{{-(1/(4 Sqrt[3]))-Sqrt[3]/4,-2}},{{23/(4 Sqrt[3])-Sqrt[3]/4,1}},{{-(19/(4 Sqrt[3]))-Sqrt[3]/4,-(1/2)}}},{{{0.2886751345948129`,-0.5`}},{{0.28867513459481287`,0.5`}},{{1.154700538379252`,0}},{{0.288675134594813`,1.5`}},{{-0.5773502691896257`,1}},{{-1.4433756729740645`,0.5`}},{{-1.4433756729740645`,-0.5`}},{{2.0207259421636907`,0.5`}},{{2.0207259421636907`,1.5`}},{{1.154700538379252`,2}},{{-2.3094010767585034`,0}},{{-0.5773502691896257`,-2}},{{2.8867513459481295`,1}},{{-3.1754264805429417`,-0.5`}}}]]},Background->GrayLevel[1, 0],ImageSize->40,Frame->False,FrameTicks->False,FrameStyle->Automatic,PlotRange->4.541451884327381`,PlotRegion->{{0,1},{0,1}},PlotRangePadding->Scaled[0.01`]],
Spacer[0],
Grid[{{Text[Style["layers:",Gray,Medium]],Text[Style[ToString@layers,Medium]]},
{Text[Style["population:",Gray,Medium]],Text[Style[ToString@population,Medium]]}},
Alignment->Left,Spacings->{.6,.2}],Spacer[0]
}}],FrameMargins->4];


(* ::Input::Initialization:: *)
layerBand[1]:={Band[{1,2},{5,6}]->1,{1,6}->1};
layerBand[l_]:=Catenate@{
{Band[{orderFromLayerT[l-1]+1,orderFromLayerT[l-1]+2},{orderFromLayerT[l]-1,orderFromLayerT[l]}]->1},
{{orderFromLayerT[l-1]+1,orderFromLayerT[l]}->1},
Table[Band[
{orderFromLayerT[l-2]+1,orderFromLayerT[l-1]+2}+n{1+2(l-2),1+2(l-1)},
{orderFromLayerT[l-2]+2l-3,orderFromLayerT[l-1]+2l-2}+n{1+2(l-2),1+2(l-1)},
{2,2}]->1,
{n,0,5}]
};
junctionBand[l_]:=Table[
Band[
{orderFromLayerT[l-1]+1,orderFromLayerT[l-1]+4(1+2(l-1))}+n{1+2(l-1),2l-1},
{orderFromLayerT[l-1]+2l-1,orderFromLayerT[l-1]+3(1+2(l-1))+1}+n{2(l-1)+1,2l-1},
{2,-2}]->1
,{n,0,2}];
TAAdjacencyMatrixT[l_]:=symmetrize@SparseArray@Catenate@Append[Table[layerBand@n,{n,1,l}],junctionBand[l]];


(* ::Subsubsection:: *)
(*Coordinates*)


(* ::Input::Initialization:: *)
TAExpandCoordsT[c_] := Module[
{i,coords = c, layer = layerFromCoordsT @ c +1, new,layerSize},

layerSize=layerSizeT@layer;

new={coords[[-1]] +{1/(Sqrt[3]),-1}};

	Table[
		AppendTo[
			new
			,
			Last[new] + Which[
				i < 1*layerSize/6,
					If[EvenQ@i,{1/(2 Sqrt[3]),1/2},{1/Sqrt[3],0}]
				,
				i < 2*layerSize/6,
					If[OddQ@i,{-1/(2 Sqrt[3]),1/2},{1/(2 Sqrt[3]),1/2}]
				,
				i < 3*layerSize/6,
					If[EvenQ@i,{-1/Sqrt[3],0},{-1/(2 Sqrt[3]),1/2}]
				,
				i < 4*layerSize/6,
					If[OddQ@i,{-1/(2 Sqrt[3]),-1/2},{-1/Sqrt[3],0}]
				,
				i < 5*layerSize/6,
					If[EvenQ@i,{1/(2 Sqrt[3]),-1/2},{-1/(2 Sqrt[3]),-1/2}]
				,
				i < 6*layerSize/6,
					If[OddQ@i,{1/Sqrt[3],0},{1/(2 Sqrt[3]),-1/2}]
			]
		]
		,
		{i, 0, layerSize-2}
	];
	Return @Join[coords,new];
]


coordinatesT=Table[If[exactCoordinates,{1/(2Sqrt[3]),-1/2},{1./(2Sqrt[3]),-1./2}] . RotationMatrix[x],{x,0,-5*Pi/3,-Pi/3}];


TACoordinatesT[layer_]:=Module[{
currentLayer=layerFromCoordsT@coordinatesT
},
If[(layer>currentLayer),coordinatesT=Nest[TAExpandCoordsT,coordinatesT,layer-currentLayer]];
Return@coordinatesT[[;;orderFromLayerT@layer]];
]


(* ::Subsubsection:: *)
(*Starting Points*)


TAStartRandomTorus[n_,d_:.5] :=TATorus[ArrayReshape[SparseArray@Transpose@{Table[Round[RandomReal[{-.5,.5}]+d],orderFromLayerT[n]]},{orderFromLayerT[n],1}]];


(* ::Subsubsection:: *)
(*Plot*)


(* ::Input::Initialization:: *)
Options[TATorusPlot]={
"ImageSize" -> Medium,
"Time"->Null,
"TimeColor"->aliveColor,
 "Padding"->Scaled[.03],
"Parallel"->False,
"Frame"->False,
"FrameStyle" ->Automatic,
"Rotation"->0,
"Translation"->{0,0},
"BorderSize"->4,
"Shadow"->False
};

TATorusPlot[grid_,OptionsPattern[]]:=Module[
{x,head,layers,stateVector,displayedVerticesOdd,displayedVerticesEven,triangleEven, triangleOdd, trianglesEven={},trianglesOdd={},color,graphicsList={},evenLayered,oddLayered,coords,imageSize=OptionValue["ImageSize" ]},


	head=grid[[0]];
	stateVector=grid[[1]];
	layers=layerFromOrderT@Length@stateVector;
coords=TACoordinatesT@layers;
	coords=coords . RotationMatrix[-OptionValue["Rotation" ]];
	coords=(#+OptionValue["Translation" ])&/@coords;

(*basic shapes*)
triangleEven=Triangle[{{-1/Sqrt[3],0},{1/(2Sqrt[3]),1/2},{1/(2Sqrt[3]),-1/2}} . RotationMatrix[-OptionValue["Rotation" ]]];
triangleOdd=Triangle[{{1/Sqrt[3],0},{-1/(2Sqrt[3]),1/2},{-(1/(2Sqrt[3])),-1/2}} . RotationMatrix[-OptionValue["Rotation" ]]];

(*select background color*)

displayedVerticesOdd=Select[Position[Normal[stateVector],{1}][[All,1]],OddQ];
displayedVerticesEven=Select[Position[Normal[stateVector],{1}][[All,1]],EvenQ];

(*compute two sets of coordinates for even and odd layers*)
trianglesOdd=coords[[displayedVerticesOdd]];
trianglesEven=coords[[displayedVerticesEven]];


(*construct Graphics object*)
AppendTo[graphicsList,aliveColor];
If[trianglesEven!={},AppendTo[graphicsList,Translate[triangleEven,trianglesEven]]];
If[trianglesOdd!={},AppendTo[graphicsList,Translate[triangleOdd,trianglesOdd]]];
If[OptionValue["Time"]=!=Null,AppendTo[graphicsList,{OptionValue["TimeColor"],Text[Style[OptionValue["Time"],FontSize->Scaled@.06],Scaled[{.98,.01}],{Right,Bottom}]}]];

PrependTo[graphicsList,{
If[OptionValue["Shadow" ],Prepend[DropShadowing[{0,0},10,White]],Identity]@
{Black,RegularPolygon[{(100+OptionValue["BorderSize" ])layers/100,Pi/2},6]},
White,RegularPolygon[{layers,Pi/2},6]
}];

(*output with different Padding options*)
If[OptionValue["ImageSize" ]=="Proportional",
imageSize=128*Norm@Last[coords];
];

Return@Graphics[
graphicsList,
Background->Transparent,
ImageSize->imageSize,
Frame->OptionValue["Frame" ],
FrameTicks->False,
FrameStyle->OptionValue["FrameStyle" ],
PlotRange->(100+OptionValue["BorderSize" ])layers/100{{-1,1},{-1,1}},
PlotRegion->{{0,1},{0,1}},
PlotRangePadding->OptionValue["Padding" ]

];

];


(* ::Subsection:: *)
(*Miscellaneous *)


TAEdit[grid_]:=DialogInput[DynamicModule[{
stateVector=grid[[1]],
layer,
coordinates,
gridQ=(Head[grid]===TAGrid),
torusQ=(Head[grid]===TATorus)
},
layer=layerFromOrder@Length@stateVector;
coordinates=Which[gridQ,TACoordinates[layer-2],torusQ,TACoordinatesT[layer]];
Column[{
ClickPane[
Dynamic@Which[
gridQ,TAGridPlot[TAGrid[stateVector],"ImageSize"->Large, "Padding"->0,"Frame"->True],
torusQ,TATorusPlot[TATorus[stateVector],"ImageSize"->Large, "Padding"->0,"Frame"->True]],
Module[{index=Nearest[coordinates->"Index",#][[1]]},
stateVector=ReplacePart[stateVector,{index,1}->Mod[stateVector[[index,1]]+1,2]]
]&],Button["Done",DialogReturn@Which[gridQ,TAGrid[stateVector],torusQ,TATorus[stateVector]]]}]
]];


TACenterColumn[rN_,n_,gr_:TAStartOneAlive]:=Module[
{grid=gr,rule=rN,i,nsteps=n,values={},shrinkQ=False},
Monitor[
	For[i=0,i<=nsteps,i++,
	AppendTo[values, TAStateVector[grid][[1,1]]];
	If[TALayer[grid]>1+nsteps-i,shrinkQ=True];
	grid=TAEvolve[grid,rule,"Shrink"->shrinkQ];
	],
Grid@{{ProgressIndicator[Dynamic[i],{0,n}],"computing grid "<>ToString[i+1]<>" of "<>ToString[n]}}];
Return@values;]

TAThreeCenterColumns[rN_,n_,gr_:TAStartOneAlive]:=Module[
{grid=gr,rule=rN,i,nsteps=n,values={},shrinkQ=False},
Monitor[
	For[i=0,i<=nsteps,i++,
	AppendTo[values, Normal@TAStateVector[grid][[{1,2,5},1]]];
	If[TALayer[grid]>3+nsteps-i,shrinkQ=True];
	grid=TAEvolve[grid,rule,"Shrink"->shrinkQ];
	],
Grid@{{ProgressIndicator[Dynamic[i],{0,n}],"computing grid "<>ToString[i+1]<>" of "<>ToString[n]}}];
Return@values;]


diagonalIndices[n_]:=Table[Which[
l<0,Floor[2+3/2 l^2],
l==0,1,
l>0,2+3/2 l(l-1)
],{l,-n,n}];


TASlice[rN_,n_,gr_:TAStartOneAlive]:=Module[
{
i,
grid=gr,
rule=If[ListQ@rN,rN,{rN}],
lines= {{0,1,0}},
newLine,lastLine
},

Monitor[
For[i=1,i<=n,i++,
grid=TAEvolve[grid,rule[[Mod[i,Length@rule,1]]]];
newLine=Normal[TAStateVector[grid][[diagonalIndices[TALayer@grid],1]]];
If[True,
newLine= Split[newLine];
newLine[[{1,-1}]]={First[newLine[[1]]]};
newLine=Catenate@newLine;
];
AppendTo[lines,newLine]
],i];

Return@lines
];


Options[TASlicePlot]={
"ImageSize" -> Medium,
"Antialiasing"->True
};

TASlicePlot[lines_, OptionsPattern[]]:=Module[
{
i,
array,
graphicsList={},
width=Max[Length/@lines],
height=Length@lines
},
array=ArrayPad[#,(width-Length[#])/2,Last[#]] &/@lines;
For[i=1,i<=Length@lines,i++,
If[First[array[[i]]]==0, 
If[!Equal@@array[[i]],AppendTo[graphicsList,{Purple,Translate[Rectangle[{0, 0}, {1, 1}],{#-1,-i}&/@Catenate[Position[array[[i]],1]]]}]],
AppendTo[graphicsList,{Purple,Rectangle[{0, -i}, {width, 1-i}]}];
If[!Equal@@array[[i]],AppendTo[graphicsList,{White,Translate[Rectangle[{0, 0}, {1, 1}],{#-1,-i}&/@Catenate[Position[array[[i]],0]]]}]]
]
];
Return@Graphics[graphicsList,
ImageSize->OptionValue["ImageSize"],
PlotRange->{{0,width},{-height,0}},
PlotRangePadding->None
];
];


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
