BeginPackage["MyPackage`"]

MyFunc::usage = "MyFunc is a useless function";

Begin["`Private`"]
MyFunc[x_] := x;

(* NonCommutativeMultiply is protected *)
Unprotect[NonCommutativeMultiply];

(* Remove the property flat from NonCommutativeMultiply *)
(* The flat property is associativity, we don't need to clear that? *)
ClearAttributes[NonCommutativeMultiply, Flat];

(*Linearity of addition:*)
NonCommutativeMultiply[H___,Plus[A_,B__],T___]:=NonCommutativeMultiply[H,A,T]+NonCommutativeMultiply[H,Plus[B],T]

(*Scalars come out. Define C[i] symbols to be scalars. *)
ScalarQ[f_]:=Or[NumericQ[f],Head[f]===C]
NonCommutativeMultiply[H___,Times[c_,A__],T___]:=c NonCommutativeMultiply[H,Times[A],T]/;ScalarQ[c]
NonCommutativeMultiply[H___,c_,T___]:=c NonCommutativeMultiply[H,T]/;ScalarQ[c];

(* Allow NonCommutativeMultiply[\[Psi]] to simplify to \[Psi] *)
NonCommutativeMultiply[H_]:=H/;(Head[H]===Symbol);

(*One-way flatness:*)
NonCommutativeMultiply[H___,NonCommutativeMultiply[M___],T___]:=NonCommutativeMultiply[H,M,T];

(*Canonical ordering:*)
NonCommutativeMultiply[H___,B_,A_,T___]:=-NonCommutativeMultiply[H,A,B,T]/;Not[OrderedQ[{B,A}]]

(*Squares vanish:*)
NonCommutativeMultiply[H___,A_,M___,A_,T___]:=0

(* 1**2**3 should simplify to 6, not Times[6,NonCommutativeMultiply[]] *) 
NonCommutativeMultiply[]:=Sequence[];

(* Reprotect NonCommutativeMultiply again *)
(*Protect[NonCommutativeMultiply];*)

End[]

EndPackage[]
