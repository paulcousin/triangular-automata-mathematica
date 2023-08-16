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
(*Evolution*)


TAEvolve::usage="TAEvolve[grid,ruleNumber] evolves the grid with the given rule.";
TANestEvolve::usage="TANestEvolve[grid,ruleNumber,n] evolves the grid n times with the given rule and returns the final grid.";
TANestListEvolve::usage="TANestListEvolve[grid,ruleNumber,n] evolves the grid n times with the given rule and returns the list of all intermediate grids.";
TANegativeGrid::usage="TANegativeGrid[grid] returns the grid with all states inverted";
TANegativeRule::usage="TANegativeRule[ruleNumber] returns the rule number that would have the same effect in a negative grid.";


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


Begin["`Private`"];


(* ::Section:: *)
(*Functions*)


(* ::Subsection::Closed:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
aliveColor=RGBColor[0.5, 0, 0.5];deadColor=GrayLevel[1];unknownColor = GrayLevel[0.85];


(* ::Subsection:: *)
(*Evolution*)


(* ::Subsubsection:: *)
(*Utilities*)


(* ::Input::Initialization:: *)
binary[n_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(23\)]\(\(IntegerDigits[n, 10, 24]\)[\([\(-1\) - i]\)]*
\*SuperscriptBox[\(2\), \(i\)]\)\);


layerFromOrder[order_]:=Ceiling[1/6 (-3+Sqrt[3(8 order-5)])];

layerFromMatrix[matrix_] := Module[{order = Max @ Dimensions @ matrix, n, t},
 Return@Ceiling@layerFromOrder[order]];

layerFromCoords[coords_] := Module[{order = Length @ coords, n, t},
Return@Ceiling@layerFromOrder[order]];
	
graphOrderFromLayer[layer_]:=1+3 (layer(layer+1))/2;


TANegativeGrid[grid_]:={grid[[1]],
SparseArray@Normal@(1-grid[[2]]),
grid[[3]]};

TANegativeRule[ruleNumber_]:=FromDigits[1-Reverse@IntegerDigits[ruleNumber,2,8],2];


(* ::Subsubsection::Closed:: *)
(*Grid Matrix and Adjacency Matrix*)


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

symmetrize[m_]:=Module[{maxDim=Max@Dimensions@m,output},
output=PadRight[m,{maxDim,maxDim}];
 Return[output+Transpose[output]];
];


(* ::Subsubsubsection::Closed:: *)
(*Expanding an existing Grid Matrix*)


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


(* ::Input::Initialization:: *)
buildMatrix[l_]:=Nest[expandMatrix,SparseArray[Automatic, {1, 1}, 0, {1, {{0, 0}, {}}, {}}],l]


(* ::Subsubsection::Closed:: *)
(*Vertex Coordinates*)


(* ::Input::Initialization:: *)
expandCoords[c_] := Module[{coords = c, layer = layerFromCoords @ c +
	 1, i},
	AppendTo[
		coords
		,
		coords[[Min[-3 * (layer - 1), -1]]] + If[OddQ[layer],
			{1 / Sqrt[3], 0}
			,
			{1 / (2 Sqrt[3]), -(1/2)}
		]
	];
	Table[
		AppendTo[
			coords
			,
			Last[coords] + Which[
				i < Floor[layer / 2],
					{0, 1}
				,
				i < Floor[layer / 2] + Ceiling[layer / 2],
					{-(Sqrt[3] / 2), 1/2}
				,
				i < 2 Floor[layer / 2] + Ceiling[layer / 2],
					{-(Sqrt[3] / 2), -(1/2)}
				,
				i < 2 Floor[layer / 2] + 2 Ceiling[layer / 2],
					{0, -1}
				,
				i < 3 Floor[layer / 2] + 2 Ceiling[layer / 2],
					{Sqrt[3] / 2, -(1/2)}
				,
				i < 3 Floor[layer / 2] + 3 Ceiling[layer / 2],
					{Sqrt[3] / 2, 1/2}
			]
		]
		,
		{i, 0, 3 * layer - 2}
	];
	Return @ coords;
]


(* ::Input::Initialization:: *)
buildCoords[l_]:=Module[{i,coords={{0,0}}},
Monitor[For[i=0,i<l,i++,coords=expandCoords@coords],i];Return@coords];


(* ::Subsubsection:: *)
(*Computing the next grid*)


(* ::Input::Initialization:: *)
ruleState[ruleNumber_,case_]:=IntegerDigits[ruleNumber,2,8][[-case-1]];


(* ::Input::Initialization:: *)
configurationVector[{aM_,sV_}]:=Normal[4*sV+aM . sV];


(* ::Input::Initialization:: *)
newState[{aM_,sV_},ruleNumber_]:= SparseArray[Map[ruleState[ruleNumber,#]&,configurationVector[{aM,sV}]]];


(* ::Input::Initialization:: *)
TAEvolve[grid_,rN_]:=Module[
{layer,
matrix=grid[[1]],
stateVector=grid[[2]],
vertexCoords=grid[[3]],
digits=IntegerDigits[rN,2,8],
lastSpecified},


layer=layerFromCoords@vertexCoords;

lastSpecified=If[#=={},0,layerFromOrder@Last[#]]&@
If[stateVector[[-1,-1]]==0,
Position[Normal@stateVector,1][[All,1]],
Position[Normal@stateVector,0][[All,1]]
];

If[lastSpecified==layer-2,
	matrix=expandMatrix@matrix;
	vertexCoords=expandCoords@vertexCoords;
	stateVector=ArrayReshape[stateVector,{Dimensions[vertexCoords][[1]],1},Normal@stateVector[[-1]]];
	layer=(layer+1);
];

stateVector=newState[{symmetrize[matrix],stateVector},rN];

stateVector=ReplacePart[stateVector,Table[{n,1}->Normal@stateVector[[-1-3layer,1]],{n,-3layer,-1}]];

Return@{matrix,stateVector, vertexCoords}
];


TANestEvolve[gr_,rN_,n_]:=Module[{i,grid=gr},
Monitor[For[i=0,i<n,i++,grid=TAEvolve[grid,rN]],Grid@{{ProgressIndicator[Dynamic[i],{0,n}],"computing grid "<>ToString[i]<>" of "<>ToString[n]}}];
Return@grid
];

TANestListEvolve[gr_,rN_,n_]:=Module[{i,grids={gr}},
Monitor[For[i=0,i<n,i++,AppendTo[grids,TAEvolve[Last@grids,rN]]],Grid@{{ProgressIndicator[Dynamic[i],{0,n}],"computing grid "<>ToString[i]<>" of "<>ToString[n]}}];
Return@grids
];


(* ::Subsection:: *)
(*Plots*)


TAConfigurationPlot[c_]:=
Graphics[{
If[c<4,White,Purple], If[c<4,EdgeForm[Thickness@Small],EdgeForm[None]],
	Triangle[{{-(1/Sqrt[3]),0},{1/(2 Sqrt[3]),1/2},{1/(2 Sqrt[3]),-(1/2)}}], 
If[Mod[c,4]>0,Purple,White],If[Mod[c,4]>0,EdgeForm[None],EdgeForm[Thickness@Small]],
	GeometricTransformation[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3],0}}],
If[Mod[c,4]>1,Purple,White],If[Mod[c,4]>1,EdgeForm[None],EdgeForm[Thickness@Small]],
	GeometricTransformation[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3]-Sqrt[3]/2,1/2}}],
If[Mod[c,4]>2,Purple,White],If[Mod[c,4]>2,EdgeForm[None],EdgeForm[Thickness@Small]],
	GeometricTransformation[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}], 1.1{{1/Sqrt[3]-Sqrt[3]/2,-(1/2)}}]
},ImageSize->{45,45}];


(* ::Input::Initialization:: *)
TATransformationPlot[before_,after_]:=
GraphicsGrid[{{
TAConfigurationPlot[before]->
Graphics[{
If[after==0,White,Purple], If[after==0,EdgeForm[Thickness@Small],EdgeForm[None]],
	Triangle[{{-(1/Sqrt[3]),0},{1/(2 Sqrt[3]),1/2},{1/(2 Sqrt[3]),-(1/2)}}],
unknownColor,EdgeForm[None],
	GeometricTransformation[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3],0}}],
	GeometricTransformation[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}],1.1{{1/Sqrt[3]-Sqrt[3]/2,1/2}}],
	GeometricTransformation[Triangle[{{1/Sqrt[3],0},{-(1/(2 Sqrt[3])),1/2},{-(1/(2 Sqrt[3])),-(1/2)}}], 1.1{{1/Sqrt[3]-Sqrt[3]/2,-(1/2)}}]
},ImageSize->{45,45}]
}},ImageSize->{128,75}];


Options[TARulePlot]={"Labeled"->False};

TARulePlot[rN_,OptionsPattern[]] :=Module[{ruleNumber=rN,graphicsGrid},
graphicsGrid=GraphicsGrid[ArrayReshape[TATransformationPlot[#,ruleState[ruleNumber,#]]&/@Range[7,0,-1],{2,4}], Frame -> True];
If[OptionValue["Labeled"],
	Return@Grid[{{Text[ruleNumber " = " BaseForm[ruleNumber, 2]]},
	{graphicsGrid}}],
	Return@graphicsGrid
];
];


(* ::Input::Initialization:: *)
Options[TAGridPlot]={"ImageSize" -> Small,"Time"->Null, "Padding"->Automatic};

TAGridPlot[grid_,OptionsPattern[]]:=Module[
{x,stateVector=grid[[2]],coords3D=grid[[3]] . Transpose[\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
FractionBox["2", 
SqrtBox["3"]], "0"},
{
RowBox[{"-", 
FractionBox["1", 
SqrtBox["3"]]}], "1"},
{
RowBox[{"-", 
FractionBox["1", 
SqrtBox["3"]]}], 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)],coords2D=grid[[3]],aliveVertices, deadVertices,triangle1, triangle2, triangles1={},triangles2={},gridstate=grid[[2]][[-1,-1]],color,graphicsList={}},

triangle1=Triangle[{{-1/Sqrt[3],0},{1/(2Sqrt[3]),1/2},{1/(2Sqrt[3]),-1/2}}];
triangle2=Triangle[{{1/Sqrt[3],0},{-1/(2Sqrt[3]),1/2},{-(1/(2Sqrt[3])),-1/2}}];

If[gridstate==0,
color=aliveColor;
aliveVertices=Drop[ArrayRules@stateVector,-1]/.({x_,1}->1)->x;
	Scan[If[
	AllTrue[coords3D[[#]],IntegerQ],
	AppendTo[triangles1,coords2D[[#]]],
	AppendTo[triangles2,coords2D[[#]]]]&,
	aliveVertices];,
color= deadColor;
deadVertices=Drop[ArrayRules[ConstantArray[1,Dimensions@stateVector]-stateVector],-1]/.({x_,1}->1)->x;
	Scan[If[
	AllTrue[coords3D[[#]],IntegerQ],
	AppendTo[triangles1,coords2D[[#]]],
	AppendTo[triangles2,coords2D[[#]]]]&,
	deadVertices];
];

triangles1=triangles1/.{x_}->{x,x,x};
triangles2=triangles2/.{x_}->{x,x,x};

AppendTo[graphicsList,color];
If[triangles1!={},AppendTo[graphicsList,Translate[triangle1,triangles1]]];
If[triangles2!={},AppendTo[graphicsList,Translate[triangle2,triangles2]]];
If[OptionValue["Time"]=!=Null,AppendTo[graphicsList,Text[Style[OptionValue["Time"],FontSize->Scaled@.06],Scaled[{.98,.01}],{Right,Bottom}]]];


Return@Which[
OptionValue["Padding"]===Automatic,
	Graphics[graphicsList,Background->If[gridstate==0,deadColor,aliveColor],
PlotRangePadding->{4,4},PlotRange->Norm@coords2D[[-1]],ImageSize->OptionValue["ImageSize"]],
OptionValue["Padding"]===None,
Graphics[graphicsList,Background->If[gridstate==0,deadColor,aliveColor],ImageSize->OptionValue["ImageSize"]]
];

]


Options[TAEvolutionPlot]={"ImageSize" -> Medium, "Timed"->True, "Pause"->True, "Animated"->True};

TAEvolutionPlot[grid_, ruleNumber_, steps_, OptionsPattern[]] := 
Module[{x,i, grids = TANestListEvolve[grid, ruleNumber, steps]},
	
	grids = TAGridPlot[grids[[#]], "ImageSize" -> OptionValue["ImageSize"],"Time"->If[OptionValue["Timed"],(#-1),Null]]& /@ Range@Length@grids;
	
	If[OptionValue["Pause"],
	For[i=0,i<Length[grids]/16,i++,
PrependTo[grids,First@grids];
AppendTo[grids,Last@grids];
AppendTo[grids,Last@grids];
]];	

	
	Return @ If[OptionValue["Animated"],ListAnimate[grids],grids];
]


(* ::Input::Initialization:: *)
TAGridPlot3D[grid_,time_]:=Module[
{x,level=time,stateVector=grid[[2]],coords3D=grid[[3]] . Transpose[\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{
FractionBox["2", 
SqrtBox["3"]], "0"},
{
RowBox[{"-", 
FractionBox["1", 
SqrtBox["3"]]}], "1"},
{
RowBox[{"-", 
FractionBox["1", 
SqrtBox["3"]]}], 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)],coords2D=grid[[3]],aliveVertices, deadVertices,shape1, shape2, shapes1={},shapes2={},shapes={},gridstate=grid[[2]][[-1,-1]],color},


shape1=MeshRegion[{{-1/Sqrt[3], 0, 0}, {1/(2Sqrt[3]), 1/2, 0}, {1/(2Sqrt[3]), -(1/2), 0}, {-1/Sqrt[3], 0, 1}, {1/(2Sqrt[3]),1/2, 1}, {1/(2Sqrt[3]), -(1/2), 1}}, {Polygon[{{1, 2, 3}, {4, 5, 6}}], Polygon[{{1, 2, 5, 4}, {2, 3, 6, 5}, {3, 1, 4, 6}}]}];

shape2=MeshRegion[{{1/Sqrt[3], 0, 0}, {-(1/(2Sqrt[3])), 1/2, 0}, {-(1/(2Sqrt[3])), -(1/2), 0}, {1/Sqrt[3], 0, 1}, {-(1/(2Sqrt[3])),1/2, 1}, {-(1/(2Sqrt[3])), -(1/2), 1}}, {Polygon[{{1, 2, 3}, {4, 5, 6}}], Polygon[{{1, 2, 5, 4}, {2, 3, 6, 5}, {3, 1, 4, 6}}]}];


If[gridstate==0,

color=aliveColor;
aliveVertices=Drop[ArrayRules@stateVector,-1]/.({x_,1}->1)->x;
Scan[If[
AllTrue[coords3D[[#]],IntegerQ],
AppendTo[shapes1,Append[coords2D[[#]],-level]],
AppendTo[shapes2,Append[coords2D[[#]],-level]]]&,
aliveVertices];

,

color= deadColor;
deadVertices=Drop[ArrayRules[ConstantArray[1,Dimensions@stateVector]-stateVector],-1]/.({x_,1}->1)->x;
Scan[If[
AllTrue[coords3D[[#]],IntegerQ],
AppendTo[shapes1,Append[coords2D[[#]],-level]],
AppendTo[shapes2,Append[coords2D[[#]],-level]]]&,
deadVertices];

];

Return@If[Union[shapes1,shapes2]=={},
EmptyRegion[3],
RegionUnion[Map[Translate[shape1,#]&,shapes1]\[Union]Map[Translate[shape2,#]&,shapes2]]
]
];


Options[TAEvolutionPlot3D]={"Mesh" -> False,"ImageSize"->Medium};

TAEvolutionPlot3D[grid_,ruleNumber_,steps_, OptionsPattern[]]:=Module[
	{regions,gridPlots},
	
	regions= TANestListEvolve[grid,ruleNumber,steps];
	gridPlots=TAGridPlot3D[regions[[#]],#]&/@Range@Length@regions;

	If[OptionValue["Mesh"],
	Return@RegionUnion@gridPlots,
	Return@Graphics3D[RegionUnion@gridPlots,
		Boxed->False,
		Method->{"ShrinkWrap" -> True},
		ViewProjection->"Orthographic",
		ViewPoint->{1,0,.75},
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

TAStartOneAlive ={SparseArray[Automatic, {4, 10}, 0, {1, {{0, 3, 5, 7, 9}, {{2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}}}, {1, 1, 1, 1, 1, 1, 1, 1, 1}}],SparseArray[Automatic, {10, 1}, 0, {1, {{0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {{1}}}, {1}}],{{0, 0}, {3^Rational[-1, 2], 0}, {3^Rational[-1, 2] + Rational[-1, 2] 3^Rational[1, 2], Rational[1, 2]}, {3^Rational[-1, 2] + Rational[-1, 2] 3^Rational[1, 2], Rational[-1, 2]}, {Rational[1, 2] 3^Rational[1, 2], Rational[-1, 2]}, {Rational[1, 2] 3^Rational[1, 2], Rational[1, 2]}, {0, 1}, {Rational[-1, 2] 3^Rational[1, 2], Rational[1, 2]}, {Rational[-1, 2] 3^Rational[1, 2], Rational[-1, 2]}, {0, -1}}};

TAStartRandom[n_] :=
	Module[{matrix = buildMatrix[n+2],coords=buildCoords[n+2]},
		Return[{matrix, 
		ArrayReshape[SparseArray@Transpose@{Table[RandomChoice[{0,1}],graphOrderFromLayer[n]]},{graphOrderFromLayer[n+2],1}],
			 coords }]
	];

TAStartLogo={buildGrid[10][[1]],
SparseArray[({#,1}->1)&/@{2,3,5,6,8,11,13,16,17,19,21,22,23,27,29,31,33,34,35,41,44,49,58,67,78,88(*,45,112*)},{166,1}],
buildGrid[10][[3]]};


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
