(* ::Package:: *)

(* ::Title:: *)
(*Triangular Automata*)


(* ::Abstract:: *)
(*This package is made to compute and display cellular automata in a triangular grid.*)


(* ::Author:: *)
(*Paul Cousin (https://orcid.org/0000-0002-3866-7615)*)


(* ::Section:: *)
(*Begin*)


BeginPackage["TriangularAutomata`"];


TAPopulation::usage="TAPopulation[grid]"
TANegative::usage="TANegative[grid]"
TANegativeRule::usage="TANegativeRule[rule]"
TADestroboscopify::usage="TADestroboscopify[rule]"

TAGrid::usage="TAGrid[states,universe,phase]"
TARandom::usage="TARandom[{x,y}] creates a x-by-y random grid.\n"<>
				"TARandom[{x,y},d] creates a random grid with density d.\n"<>
				"TARandom[{x,y},d,u] creates a random grid with the universe in state u.\n"<>
				"TARandom[{x,y},d,u,p] creates arandom grid with phase p."
				
TAPad::usage="TAPad[grid]"
TAEdit::usage="TAEdit[grid]"

TAEvolve::usage="TAEvolve[rule][grid] evolves the grid once.\n"<>
				"TAEvolve[grid,rule] evolves the grid once.\n"<>
				"TAEvolve[grid,rule,n] evolves the grid n times.\n"<>
				"TAEvolve[grid,{\!\(\*SubscriptBox[\(r\), \(1\)]\),...,\!\(\*SubscriptBox[\(r\), \(m\)]\)},n] evolves the grid n times alternating between the different rules.";

TAPlot::usage="TAPlot[grid] plots grid."
TAQuickPlot::usage="TAQuickPlot[grid] plots grid efficiently."

TAConfigurationPlot::usage="TAQuickPlot[configuration]"
TARulePlot::usage="TARulePlot[rule]"



Begin["`Private`"];


(* ::Section:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
aliveColor=RGBColor[0.5, 0, 0.5];deadColor=GrayLevel[1];unknownColor = GrayLevel[0.85];
exactCoordinates= False;
onlineDirectory="https://files.paulcousin.net/triangular-automata";
localDirectory=FileNameJoin@{Directory[],"Triangular Automata"};


(* ::Section:: *)
(*Utilities*)


TAPopulation[g_TAGrid]:=Total[Abs[g[[1]]-g[[1,1,1]]],2];
TANegative[g_TAGrid]:=TAGrid[1-g[[1]],1-g[[2]],g[[3]]];
TANegativeRule[ruleNumber_]:=FromDigits[1-Reverse@IntegerDigits[ruleNumber,2,8],2];
TADestroboscopify[n_]:={255-n,FromDigits[Reverse@IntegerDigits[n,2,8],2]};


(* ::Section:: *)
(*Grid*)


TAGrid[states_?MatrixQ]:=TAGrid[states,0];
TAGrid[states_?MatrixQ,universe_Integer]:=TAGrid[states,universe,0];
TAGrid[states_?MatrixQ,universe_Integer,phase_Integer]:=TAGrid[states,universe,phase,(Dimensions[states][[{2,1}]]-1)*{-Sqrt[3]/4,1/4}];


TAGrid[1]:=TAGrid[ArrayPad[{{1}},2,0]];
TAGrid[0]:=TAGrid[ArrayPad[{{0}},2,1],1];
TAGrid["Hexagon"]:=TAGrid[{{1,1},{1,1},{1,1}}];
TAGrid["TA"]:=TAGrid[{
	{0,0,0,0,1,1,0,0,0},
	{0,0,0,0,1,1,1,0,0},
	{0,0,1,1,1,0,1,1,0},
	{0,1,1,0,1,1,1,1,0},
	{1,1,0,0,1,1,0,0,0},
	{1,0,1,0,1,0,0,0,0},
	{0,0,1,0,1,0,0,0,0},
	{0,0,0,0,0,0,0,0,0},
	{0,0,0,1,0,0,0,0,0}}];


TARandom[{x_Integer,y_Integer},d_:0.5,u_:0,p_:0]:=TAGrid[RandomVariate[BernoulliDistribution[Abs[d-u]],{y,x}],u,p];


 TAGrid/:MakeBoxes[TAGrid[states_?MatrixQ,universe_,phase_Integer,coords_List],form_]:=With[
 {
 interpretation=Interpretation[
Which[
	Times@@Dimensions@states<=10000,
	Deploy@TAPlot[TAGrid[states,universe,phase,coords],ImageSize->Tiny,Frame->(universe===0)],
	Times@@Dimensions@states<=2000000,
	Deploy@TAQuickPlot[TAGrid[states,universe,phase,coords],ImageSize->Tiny,Frame->(universe===0)],
	True,
	Deploy@largeGridRepresentation[Dimensions@states,Total[Abs[states[[1,1]]-states],2],If[universe==0,White,Purple]]
 ],TAGrid[states,universe,phase,coords]]
 },
 MakeBoxes[interpretation,form]];


largeGridRepresentation[dimensions_,population_,universe_]:=Panel[Grid[{{
Spacer[0],
	Graphics[{GrayLevel[0.8], Translate[Triangle[{{-3^Rational[-1, 2], 0}, {Rational[1, 2] 3^Rational[-1, 2], Rational[1, 2]}, {Rational[1, 2] 3^Rational[-1, 2], Rational[-1, 2]}}], {{{0.8660254037844388, -0.5}}, {{0.8660254037844388, 0.5}}, {{-0.8660254037844384, 0.5}}, {{1.7320508075688774`, 0.}}, {{1.7320508075688774`, 1.}}, {{0.8660254037844388, 1.5}}, {{-1.732050807568877, 0.}}, {{-0.8660254037844384, -1.5}}, {{0.8660254037844388, -1.5}}, {{2.5980762113533165`, 0.5}}, {{-2.598076211353315, -0.5}}, {{3.4641016151377553`, 0.}}}], Translate[Triangle[{{3^Rational[-1, 2], 0}, {Rational[-1, 2] 3^Rational[-1, 2], Rational[1, 2]}, {Rational[-1, 2] 3^Rational[-1, 2], Rational[-1, 2]}}], {{{0.5773502691896258, 0.}}, {{-0.28867513459481275`, 0.5}}, {{1.4433756729740645`, -0.5}}, {{0.577350269189626, 1.}}, {{-1.1547005383792512`, 0.}}, {{-1.1547005383792512`, -1.}}, {{0.577350269189626, -1.}}, {{2.3094010767585034`, 0.}}, {{2.3094010767585034`, 1.}}, {{1.4433756729740648`, 1.5}}, {{-2.02072594216369, -0.5}}, {{-0.28867513459481264`, -2.5}}, {{3.175426480542942, 0.5}}, {{-2.8867513459481273`, -1.}}}]}, Background -> GrayLevel[0, 0], ImageSize -> 50],
Spacer[0],
	Grid[{{Text[Style["dimensions:",Gray,Medium]],Text[Style[ToString@dimensions,Medium]]},
		{Text[Style["population:",Gray,Medium]],Text[Style[ToString@population,Medium]]},
		{Text[Style["universe:",Gray,Medium]],universe}},
		Alignment->Left,Spacings->{.6,.2}],
Spacer[0]
}}],FrameMargins->4];


TAPad[TAGrid[s_?MatrixQ,u_Integer,p_Integer,c_List]]:=Module[{states=s,universe=u,phase=p,coords=c,pad},
pad[n_]:={Boole@*MemberQ[1-universe]/@states[[{n,-n},;;]],Boole@*MemberQ[1-universe]/@Transpose@states[[;;,{n,-n}]]};
pad=2pad[1]+pad[2]/.x_Integer:>Min[2,x];
states=ArrayPad[states,pad,universe];
phase=Mod[phase+Boole@OddQ[pad[[1,1]]+pad[[2,1]]],2];
coords+=pad[[{2,1},1]]*{-Sqrt[3]/2,1/2};
Return@TAGrid[states,universe,phase,coords];
];


Options[TAEdit]={ImageSize -> Large};

TAEdit[g_TAGrid,OptionsPattern[]]:=DialogInput[DynamicModule[
{grid=g,firstCoord=g[[4]],phase=g[[3]],dims=Dimensions@g[[1]],shift=1/(4*Sqrt[3]),coords},
coords=Transpose@CoordinateBoundsArray[{{firstCoord[[1]],firstCoord[[1]]+(dims[[2]]-1)Sqrt[3]/2},{firstCoord[[2]],firstCoord[[2]]-(dims[[1]]-1)/2}},{Sqrt[3]/2,-1/2}];
coords[[1+phase;;;;2,1;;;;2]]=Map[#+{shift,0}&,coords[[1+phase;;;;2,1;;;;2]],{2}];
coords[[2-phase;;;;2,1;;;;2]]=Map[#-{shift,0}&,coords[[2-phase;;;;2,1;;;;2]],{2}];
coords[[2-phase;;;;2,2;;;;2]]=Map[#+{shift,0}&,coords[[2-phase;;;;2,2;;;;2]],{2}];
coords[[1+phase;;;;2,2;;;;2]]=Map[#-{shift,0}&,coords[[1+phase;;;;2,2;;;;2]],{2}];
Column[{
ClickPane[
Dynamic@TAPlot[grid,ImageSize->OptionValue[ImageSize],PlotRange->Full,PlotRangePadding->None,Frame->True],
Module[{index=Nearest[Catenate[coords]->"Index",#][[1]],x,y},
x=1+Mod[index-1,dims[[2]]];y=1+Quotient[index-1,dims[[2]]];
grid[[1,y,x]]=1-grid[[1,y,x]]
]&],Button["Done",DialogReturn@grid]}]
]];


(* ::Section:: *)
(*Evolution*)


TAEvolve[r_Integer]=TAEvolve[#,r]&


TAEvolve[g_TAGrid,r_,n_Integer]:=Module[{grid=g,i,rule=If[ListQ@r,r,{r}]},
Monitor[
	Do[grid=TAEvolve[grid,rule[[Mod[i,Length[rule],1]]]],{i,n}],
	Grid[{{ProgressIndicator[i,{0,n}],"step "<>ToString[i]<>" of "<>ToString[n]}}]
];
Return@grid;
];


TAEvolve[g_TAGrid,r_,{n_Integer}]:=Module[{grids={g},i,rule=If[ListQ@r,r,{r}]},
Monitor[
	Do[AppendTo[grids,TAEvolve[Last@grids,rule[[Mod[i,Length@rule,1]]]]],{i,n}],
	Grid[{{ProgressIndicator[i,{0,n}],"step "<>ToString[i]<>" of "<>ToString[n]}}]
];
Return@grids;
];


Options[TAEvolve]={"Pad"->True};

TAEvolve[g_TAGrid,r_Integer,OptionsPattern[]]:=Module[
{grid=If[OptionValue["Pad"],TAPad[g],g],states,universe,phase,coords,digits,config,dims,padRight},
states=grid[[1]];universe=grid[[2]];phase=grid[[3]];coords=grid[[4]];

dims=Dimensions@states;
digits=Reverse@IntegerDigits[r,2,8];

If[EvenQ@dims[[2]],
	states=ArrayPad[states,{{0,0},{0,1}},universe];
	dims=Dimensions@states;
];

If[IntegerQ@universe,universe=If[universe===0,First,Last]@digits];

config=Flatten[4*states+RotateLeft@states+RotateRight@states];
states//=Flatten;
config[[1+phase;;;;2]]+=RotateLeft[states][[1+phase;;;;2]];
config[[2-phase;;;;2]]+=RotateRight[states][[2-phase;;;;2]];
states=ArrayReshape[digits[[config+1]],dims];

Return@TAGrid[states,universe,phase,coords];

];


(* ::Section:: *)
(*Plot*)


Options[TAPlot]={ImageSize->Medium,Frame->False,FrameTicks->None,PlotRange->Automatic,PlotRangePadding->Scaled[.03]};

TAPlot[TAGrid[s_?MatrixQ,u_Integer,p_Integer,c_List],OptionsPattern[]]:=Module[{
states=s,universe=If[IntegerQ@u,u,0],phase=p,
triangleLeft=Triangle[#+c&/@{{-Sqrt[3]/4,0},{Sqrt[3]/4,1/2},{Sqrt[3]/4,-1/2}}],
triangleRight=Triangle[#+c&/@{{Sqrt[3]/4,0},{-Sqrt[3]/4,1/2},{-Sqrt[3]/4,-1/2}}],
graphics,positions
},

graphics={If[universe===1,White,Purple]};
positions=Position[states,1-If[IntegerQ@universe,universe,1]]-1;

positions=Select[positions,#@*Total]&/@{EvenQ,OddQ}/.{y_Integer,x_Integer}:>{x*Sqrt[3]/2,-y/2};
If[phase==1,positions//=Reverse];

If[positions[[1]]!={},AppendTo[graphics,Translate[triangleLeft,positions[[1]]]]];
If[positions[[2]]!={},AppendTo[graphics,Translate[triangleRight,positions[[2]]]]];

Return@Graphics[graphics,
	Background->If[universe===1,Purple,White],
	ImageSize->OptionValue[ImageSize],
	Frame->OptionValue[Frame],
	FrameTicks->OptionValue[FrameTicks],
	PlotRange->If[OptionValue[PlotRange]==Full,
		Transpose@{c-{Sqrt[3]/4,-1/2},c+Dimensions[states][[{2,1}]]*{Sqrt[3]/2,-1/2}-{Sqrt[3]/4,0}},
		Evaluate@OptionValue[PlotRange]
	],
	PlotRangePadding->OptionValue[PlotRangePadding]
];
];


Options[TAQuickPlot]={ImageSize->Medium,Frame->False};

TAQuickPlot[g_TAGrid,OptionsPattern[]]:=Module[{states=g[[1]],dims=Dimensions@g[[1]]},
ArrayPlot[states,
	ColorRules->{1->Purple,0->White},
	AspectRatio->dims[[1]]/(Sqrt[3]*dims[[2]]),
	ImageSize->OptionValue[ImageSize],
	Frame->OptionValue[Frame]
]];


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
	TATransformationPlot[#,IntegerDigits[ruleNumber,2,8][[-#-1]]],
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


(* ::Section:: *)
(*End*)


End[];


EndPackage[];
