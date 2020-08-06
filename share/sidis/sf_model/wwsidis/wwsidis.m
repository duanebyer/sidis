(* ::Package:: *)

 


Print["Package WW-SIDIS contains the set of TMDs calculated with WW approximation and SIDIS structure functions"];
Print["Copyright: Alexei Prokudin (PSU Berks), Kemal Tezgin (UConn), Version 1.01 (06/05/2018)"];
Print["e-mail: prokudin@jlab.org, kemal.tezgin@uconn.edu"];
Print["https://github.com/prokudin/WW-SIDIS"];
Print[Style["If you use this package, please, site Bastami:2018xqd for this paper and https://github.com/prokudin/WW-SIDIS", Red]];
Print["___________________________________________________________________________"];
Print["Contains the following functions: "];
BeginPackage["wwsidis`"];
f1u::usage="f1u[x_,Q2_ ] is the unpolarised collinear PDF for u quark, bibitem{Martin:2009iq}";
f1d::usage="f1d[x_,Q2_ ] is the unpolarised collinear PDF for d quark, bibitem{Martin:2009iq}";
f1ubar::usage="f1ubar[x_,Q2_ ] is the unpolarised collinear PDF for ubar quark, bibitem{Martin:2009iq}";
f1dbar::usage="f1dbar[x_,Q2_ ] is the unpolarised collinear PDF for dbar quark, bibitem{Martin:2009iq}";
f1s::usage="f1s[x_,Q2_ ] is the unpolarised collinear PDF for s quark, bibitem{Martin:2009iq}";
f1sbar::usage="f1sbar[x_,Q2_ ] is the unpolarised collinear PDF for sbar quark, bibitem{Martin:2009iq}";
?f1u;
?f1d;
?f1ubar;
?f1dbar;
?f1s;
?f1sbar;
D1u::usage="D1u[pion_, z_, Q2_ ] is the unpolarised collinear FF for u quark, pion = pi+ or pi-, bibitem{deFlorian:2007aj}";
D1d::usage="D1d[pion_, z_, Q2_ ] is the unpolarised collinear FF for dbar quark, pion = pi+ or pi-, bibitem{deFlorian:2007aj}";
D1ubar::usage="D1ubar[pion_, z_, Q2_ ] is the unpolarised collinear FF for ubar quark, pion = pi+ or pi-, bibitem{deFlorian:2007aj}";
D1dbar::usage="D1dbar[pion_, z_, Q2_ ] is the unpolarised collinear FF for dbar quark, pion = pi+ or pi-, bibitem{deFlorian:2007aj}";
D1s::usage="D1s[pion_, z_, Q2_ ] is the unpolarised collinear FF for s quark, pion = pi+ or pi-, bibitem{deFlorian:2007aj}";
D1sbar::usage="D1sbar[pion_, z_, Q2_ ] is the unpolarised collinear FF for sbar quark, pion = pi+ or pi-, bibitem{deFlorian:2007aj}";
?D1u;
?D1d;
?D1ubar;
?D1dbar;
?D1s;
?D1sbar;
g1u::usage="g1u[x_,Q2_ ] is the helicity collinear PDF for u quark, bibitem{Gluck:1998xa}";
g1d::usage="g1d[x_,Q2_ ] is the helicity collinear PDF for d quark, bibitem{Gluck:1998xa}";
g1ubar::usage="g1ubar[x_,Q2_ ] is the helicity collinear PDF for ubar quark, bibitem{Gluck:1998xa}";
g1dbar::usage="g1dbar[x_,Q2_ ] is the helicity collinear PDF for dbar quark, bibitem{Gluck:1998xa}";
g1s::usage="g1s[x_,Q2_ ] is the helicity collinear PDF for s quark, bibitem{Gluck:1998xa}";
g1sbar::usage="g1sbar[x_,Q2_ ] is the helicity collinear PDF for sbar quark, bibitem{Gluck:1998xa}";
?g1u;
?g1d;
?g1ubar;
?g1dbar;
?g1s;
?g1sbar;
h1u::usage="h1u[x_,Q2_ ] is the transversity collinear PDF for u quark, bibitem{Anselmino:2013vqa}";
h1d::usage="h1d[x_,Q2_ ] is the transversity collinear PDF for d quark, bibitem{Anselmino:2013vqa}";
h1ubar::usage="h1ubar[x_,Q2_ ] is the transversity collinear PDF for ubar quark, bibitem{Anselmino:2013vqa}";
h1dbar::usage="h1dbar[x_,Q2_ ] is the transversity collinear PDF for dbar quark, bibitem{Anselmino:2013vqa}";
h1s::usage="h1s[x_,Q2_ ] is the transversity collinear PDF for s quark, bibitem{Anselmino:2013vqa}";
h1sbar::usage="h1sbar[x_,Q2_ ] is the transversity collinear PDF for sbar quark, bibitem{Anselmino:2013vqa}";
?h1u;
?h1d;
?h1ubar;
?h1dbar;
?h1s;
?h1sbar;
H1perpFirstMoment::usage="H1perpFirstMoment[quark_, pion_, x_, Q2_] the first moment of Collins FF favoured, bibitem{Anselmino:2013vqa}";
?H1perpFirstMoment;

f1TperpuFirstMoment::usage="f1TperpuFirstMoment[x_,Q2_] the first moment of \!\(\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)]\) Sivers PDF for u quark, bibitem{Anselmino:2011gs}";
f1TperpdFirstMoment::usage="f1TperpdFirstMoment[x_,Q2_] the first moment of \!\(\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)]\) Sivers PDF for d quark, bibitem{Anselmino:2011gs}";
f1TperpubarFirstMoment::usage="f1TperpubarFirstMoment[x_,Q2_] the first moment of \!\(\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)]\) Sivers PDF for ubar quark, bibitem{Anselmino:2011gs}";
f1TperpdbarFirstMoment::usage="f1TperpdbarFirstMoment[x_,Q2_] the first moment of \!\(\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)]\) Sivers PDF for dbar quark, bibitem{Anselmino:2011gs}";
f1TperpsFirstMoment::usage="f1TperpsFirstMoment[x_,Q2_] the first moment of \!\(\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)]\) Sivers PDF for s quark, bibitem{Anselmino:2011gs}";
f1TperpsbarFirstMoment::usage="f1TperpsbarFirstMoment[x_,Q2_] the first moment of \!\(\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)]\) Sivers PDF for sbar quark, bibitem{Anselmino:2011gs}";
?f1TperpuFirstMoment;
?f1TperpdFirstMoment;
?f1TperpubarFirstMoment;
?f1TperpdbarFirstMoment;
?f1TperpsFirstMoment;
?f1TperpsbarFirstMoment;
h1TperpuSecondMoment::usage="h1TperpuSecondMoment[x_,Q2_] the second moment of pretzelosity for u quark, bibitem{Lefky:2014eia};"
h1TperpdSecondMoment::usage="h1TperpdSecondMoment[x_,Q2_] the second moment of pretzelosity for d quark, bibitem{Lefky:2014eia};"
?h1TperpuSecondMoment;
?h1TperpdSecondMoment;
h1perpuFirstMoment::usage="h1perpuFirstMoment[x_,Q2_] the first moment of Boer-Mulders function for u quark, bibitem{Barone:2009hw};"
h1perpdFirstMoment::usage="h1perpdFirstMoment[x_,Q2_] the first moment of Boer-Mulders function for d quark, bibitem{Barone:2009hw};"
?h1perpuFirstMoment;
?h1perpdFirstMoment;
h1Lu::usage="h1Lu[x_,Q2_ ] is the \!\(\*SubsuperscriptBox[\(h\), \(1  L\), \(\[UpTee]\)]\) PDF for u quark, this work";
h1Ld::usage="h1Ld[x_,Q2_ ] is the \!\(\*SubsuperscriptBox[\(h\), \(1  L\), \(\[UpTee]\)]\) PDF for u quark, this work";
?h1Lu;
?h1Ld;
g1Tperpu::usage="g1Tperpu[x_,Q2_ ] is the \!\(\*SubsuperscriptBox[\(g\), \(1  T\), \(\[UpTee]\)]\) PDF for u quark, this work";
g1Tperpd::usage="g1Tperpd[x_,Q2_ ] is the \!\(\*SubsuperscriptBox[\(g\), \(1  T\), \(\[UpTee]\)]\) PDF for d quark, this work";
?g1Tperpu;
?g1Tperpd;
ALL::usage="ALL[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubscriptBox[\(A\), \(LL\)]\) asymmetry";
?ALL;
AUTSivers::usage="AUTSivers[pion_, x_, z_, Q2_, PT_] is the Sivers \!\(\*SubsuperscriptBox[\(A\), \(UT\), \(TraditionalForm\`sin[\*SubscriptBox[\(\[CapitalPhi]\), \(h\)] - \*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?AUTSivers;
AUTCollins::usage="AUTCollins[pion_, x_, z_, Q2_, PT_] is the Collins \!\(\*SubsuperscriptBox[\(A\), \(UT\), \(TraditionalForm\`sin[\*SubscriptBox[\(\[CapitalPhi]\), \(h\)] + \*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?AUTCollins;
AUTh1tp::usage="AUTh1tp[pion_, x_, z_, Q2_, PT_] is the Pretzelosity \!\(\*SubsuperscriptBox[\(A\), \(UT\), \(TraditionalForm\`sin[3\\\ \*SubscriptBox[\(\[CapitalPhi]\), \(h\)] - \*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?AUTh1tp;
AUUcos2phi::usage="AUUcos2phi[pion_, x_, z_, Q2_, PT_] is the Boer-Mulders \!\(\*SubsuperscriptBox[\(A\), \(UU\), \(TraditionalForm\`cos[2 \*SubscriptBox[\(\[CapitalPhi]\), \(h\)]]\)]\) asymmetry";
?AUUcos2phi;
ALT::usage="ALT[pion_, x_, z_, Q2_, PT_] is the Double Spin \!\(\*SubsuperscriptBox[\(A\), \(LT\), \(TraditionalForm\`cos[\*SubscriptBox[\(\[CapitalPhi]\), \(h\)] - \*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?ALT;
AULsin2phi::usage="AULsin2phi[pion_, x_, z_, Q2_, PT_] is the Kotzinian-Mulders \!\(\*SubsuperscriptBox[\(A\), \(UL\), \(TraditionalForm\`sin[2 \*SubscriptBox[\(\[CapitalPhi]\), \(h\)]]\)]\) asymmetry";
?AULsin2phi;
ALTcosPhi::usage="ALTcosPhi[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubsuperscriptBox[\(A\), \(LT\), \(TraditionalForm\`cos[\*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?ALTcosPhi;
AULsinphi::usage="AULsinphi[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubsuperscriptBox[\(A\), \(UL\), \(TraditionalForm\`sin[\*SubscriptBox[\(\[CapitalPhi]\), \(h\)]]\)]\) asymmetry";
?AULsinphi;
ALLcosphi::usage="ALLcosphi[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubsuperscriptBox[\(A\), \(LL\), \(TraditionalForm\`cos[\*SubscriptBox[\(\[CapitalPhi]\), \(h\)]]\)]\) asymmetry";
?ALLcosphi;
ALTcos2PhiminusPhiS::usage="ALTcos2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubsuperscriptBox[\(A\), \(LT\), \(TraditionalForm\`cos[2 \*SubscriptBox[\(\[CapitalPhi]\), \(h\)] - \*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?ALTcos2PhiminusPhiS;
AUUcosphi::usage="AUUcosphi[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubsuperscriptBox[\(A\), \(UU\), \(TraditionalForm\`cos[\*SubscriptBox[\(\[CapitalPhi]\), \(h\)]]\)]\) asymmetry";
?AUUcosphi;
AUTsinphiS::usage="AUTsinphiS[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubsuperscriptBox[\(A\), \(UT\), \(TraditionalForm\`sin[\*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?AUTsinphiS;
AUTsin2PhiminusPhiS::usage="AUTsin2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] is the \!\(\*SubsuperscriptBox[\(A\), \(UT\), \(TraditionalForm\`sin[2 \*SubscriptBox[\(\[CapitalPhi]\), \(h\)] - \*SubscriptBox[\(\[CapitalPhi]\), \(S\)]]\)]\) asymmetry";
?AUTsin2PhiminusPhiS;
(*CrossSectionUU::usage="CrossSectionUU[pion_,x_,z_,Q2_,PT_,energy_,phi_] calculates the differential cross section of unpolarized beam and target";
CrossSectionUL::usage="CrossSectionUL[pion_,x_,z_,Q2_,PT_,energy_,phi_] calculates the differential cross section of unpolarized beam and longitidunally polarized target";
CrossSectionLU::usage="CrossSectionLU[pion_,x_,z_,Q2_,PT_,energy_,phi_,helicity_] calculates the differential cross section of longitidunally polarized beam and unpolarized target";
CrossSectionLL::usage="CrossSectionLL[pion_,x_,z_,Q2_,PT_,energy_,phi_,helicity_] calculates the differential cross section of longitidunally polarized beam and longitidunally polarized target";
CrossSectionUT::usage="CrossSectionUT[pion_,x_,z_,Q2_,PT_,energy_,phih_,phiS_,helicity_] calculates the differential cross section of unpolarized beam and transversely polarized target";
CrossSectionLT::usage="CrossSectionLT[pion_,x_,z_,Q2_,PT_,energy_,phih_,phiS_,helicity_] calculates the differential cross section of longitidunally polarized beam and transversely polarized target";
*)
CrossSection::usage="Semi-Inclusive Deep Inelastic Cross Section, CrossSection[pion_,x_,z_,Q2_,PT_,energy_,phih_,phiS_,helicity_,targetpolarization_] calculates the differential cross section of specified pion_,x_,z_,Q2_,PT_,energy_,phih_,phiS_,helicity_,targetpolarization_";
?CrossSection;
Print["Contains the following constants: "];
avp::usage="avp unpolarized \!\(\*SubscriptBox[\(D\), \(1\)]\) TMD fragmentation distribution width";
?avp;
avk::usage="avk unpolarized \!\(\*SubscriptBox[\(f\), \(1\)]\) TMD  width";
?avk;
avkg::usage="avkg helicity \!\(\*SubscriptBox[\(g\), \(1\)]\)TMD width";
?avkg;
avks::usage="avks \!\(\*SubsuperscriptBox[\(f\), \(1  T\), \(\[UpTee]\)]\) (Sivers function) TMD width";
?avks;
avph::usage="avph \!\(\*SubsuperscriptBox[\(H\), \(1\), \(\[UpTee]\)]\) (Collins fragmentation function) TMD width";
?avph;
Begin["`Private`"];
DSShplus= ReadList["./Grids/fragmentationpiplus.dat",Real,RecordLists-> True];
DSShminus= ReadList["./Grids/fragmentationpiminus.dat",Real,RecordLists-> True];
uhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,3]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}]
dhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,4]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
shplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,5]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
ubhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,6]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
dbhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,7]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
sbhplus=Interpolation[Table[{{DSShplus[[i,1]],DSShplus[[i,2]]},DSShplus[[i,8]]},{i,1,Length[DSShplus]}],InterpolationOrder->{3,3}];
uhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,3]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}]
dhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,4]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
shminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,5]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
ubhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,6]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
dbhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,7]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
sbhminus=Interpolation[Table[{{DSShminus[[i,1]],DSShminus[[i,2]]},DSShminus[[i,8]]},{i,1,Length[DSShminus]}],InterpolationOrder->{3,3}];
<<mstwpdf.m;
prefix="Grids/mstw2008lo";
Timing[ReadPDFGrid[prefix,0]];
upv[x_,Q2_]:=xf[0,x,Sqrt[Q2],8]/x;
dnv[x_,Q2_]:=xf[0,x,Sqrt[Q2],7]/x;
usea[x_,Q2_]:=xf[0,x,Sqrt[Q2],-2]/x;
dsea[x_,Q2_]:=xf[0,x,Sqrt[Q2],-1]/x;
str[x_,Q2_]:=xf[0,x,Sqrt[Q2],3]/x;
sbar[x_,Q2_]:=xf[0,x,Sqrt[Q2],-3]/x;
up[x_,Q2_]:=upv[x,Sqrt[Q2]]+usea[x,Sqrt[Q2]];
dn[x_,Q2_]:=dnv[x,Sqrt[Q2]]+dsea[x,Sqrt[Q2]];
upbar[x_,Q2_]:=usea[x,Sqrt[Q2]];
dnbar[x_,Q2_]:=dsea[x,Sqrt[Q2]];
g1param= ReadList["./Grids/g1.dat",Real,RecordLists-> True];
g1u=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,3]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1d=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,4]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1s=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,5]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1ubar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,6]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1dbar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,7]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
g1sbar=Interpolation[Table[{{g1param[[i,1]],g1param[[i,2]]},g1param[[i,8]]},{i,1,Length[g1param]}],InterpolationOrder->{3,3}];
sb= ReadList["./Grids/SofferBound.dat",Real,RecordLists-> True];
sb1u=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,3]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1d=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,4]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1s=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,5]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1ubar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,6]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1dbar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,7]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
sb1sbar=Interpolation[Table[{{sb[[i,1]],sb[[i,2]]},sb[[i,8]]},{i,1,Length[sb]}],InterpolationOrder->{3,3}];
(*2005 fit Appendix A.1 [hep-ph/0501196]*)
Clear[avk];
avk=0.25; 
f1u[x_,Q2_ ]:= up[x,Q2];
f1d[x_,Q2_ ]:= dn[x,Q2] ;
f1ubar[x_,Q2_ ]:= upbar[x,Q2];  
f1dbar[x_,Q2_ ]:= dnbar[x,Q2] ; 
f1s[x_,Q2_ ]:= str[x,Q2];
f1sbar[x_,Q2_ ]:= sbar[x,Q2];


f1uTMD[x_,Q2_ ,kt_]:= up[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];
f1dTMD[x_,Q2_,kt_ ]:= dn[x,Q2] 1/(\[Pi] avk) Exp[-kt^2/avk];
f1ubarTMD[x_,Q2_,kt_ ]:= upbar[x,Q2]  1/(\[Pi] avk) Exp[-kt^2/avk];
f1dbarTMD[x_,Q2_ ,kt_]:= dnbar[x,Q2]  1/(\[Pi] avk) Exp[-kt^2/avk];
f1sTMD[x_,Q2_,kt_ ]:= str[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];
f1sbarTMD[x_,Q2_,kt_ ]:= sbar[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avk];

(*2005 fit Appendix A.1 [hep-ph/0501196]*)
Clear[avp];
avp=0.2;
 
D1uTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", uhplus[z, Q2],If[pion == "pi-", uhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1dTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", dhplus[z, Q2],If[pion == "pi-", dhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1ubarTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", ubhplus[z, Q2],If[pion == "pi-", ubhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1dbarTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", dbhplus[z, Q2],If[pion == "pi-", dbhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1sTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", shplus[z, Q2],If[pion == "pi-", shminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
D1sbarTMD[pion_, z_, Q2_, pt_] := If[pion == "pi+", sbhplus[z, Q2],If[pion == "pi-", sbhminus[z, Q2]]] 1/(\[Pi] avp) Exp[-pt^2/avp]
 
 
D1u[pion_, z_, Q2_ ] := If[pion == "pi+", uhplus[z, Q2], If[pion == "pi-", uhminus[z, Q2]]]
D1d[pion_, z_, Q2_ ] := If[pion == "pi+", dhplus[z, Q2], If[pion == "pi-", dhminus[z, Q2]]]
D1ubar[pion_, z_, Q2_ ] := If[pion == "pi+", ubhplus[z, Q2],If[pion == "pi-", ubhminus[z, Q2]]]
D1dbar[pion_, z_, Q2_ ] := If[pion == "pi+", dbhplus[z, Q2],If[pion == "pi-", dbhminus[z, Q2]]]
D1s[pion_, z_, Q2_ ] := If[pion == "pi+", shplus[z, Q2],If[pion == "pi-", shminus[z, Q2]]]
D1sbar[pion_, z_, Q2_ ] := If[pion == "pi+", sbhplus[z, Q2], If[pion == "pi-", sbhminus[z, Q2]]]
 
(*Lattice fit Appendix A. 2 [0908.1283].*)
Clear[avkg];
avkg=avk 0.76 ;


(* Helicity functions with DGLAP Evolution*)


g1uTMD[x_, Q2_, kt_] := g1u[x, Q2] 1/(\[Pi] avkg) Exp[-kt^2/avkg];
g1dTMD[x_, Q2_, kt_] := g1d[x, Q2] 1/(\[Pi] avkg) Exp[-kt^2/avkg];
 
(*2013 fit Appendix A.4 [1303.3822].*)
Clear[NuT,NdT,alphaT,betaT];
NuT=0.46;
NdT=-1.000;
alphaT=1.11;
betaT=3.64;

(* Transversity function with DGLAP Evolution*)

h1u[x_, Q2_] := (NuT x^alphaT (1 - x)^betaT (alphaT + betaT)^(alphaT + betaT)/((alphaT^alphaT) (betaT^betaT))) sb1u[x, Q2];
h1d[x_, Q2_] := (NdT x^alphaT (1 - x)^betaT (alphaT + betaT)^(alphaT + betaT)/((alphaT^alphaT) (betaT^betaT))) sb1d[x, Q2];
h1ubar[x_, Q2_] := 0.;
h1dbar[x_, Q2_] := 0.;
h1s[x_, Q2_] := 0.;
h1sbar[x_, Q2_] := 0.;

h1uTMD[x_, Q2_, kt_] := h1u[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avk];
h1dTMD[x_, Q2_, kt_] := h1d[x, Q2] 1/(\[Pi] avk) Exp[-kt^2/avk];

(*2013 fit Appendix A.4 [1303.3822].*)
Clear[MC,NCfav,NCunf,alphaC,betaC, Mh];
MC=Sqrt[1.5];
NCfav=0.49;
NCunf=-1.000;
alphaC=1.06;
betaC=0.07;
Mh = 0.135;

avph = avp MC^2/(avp+MC^2); (*collins FF width*)

(*Collins function with DGLAP Evolution*)

H1perpFavHalfMoment[x_, Q2_] := Sqrt[2 E] (NCfav x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1u["pi+", x, Q2]/4 Sqrt[(\[Pi] avp)/(MC^2 + avp)^3];
H1perpUnfHalfMoment[x_, Q2_] := Sqrt[2 E] (NCunf x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1d["pi+", x, Q2]/4 Sqrt[(\[Pi] avp)/(MC^2 + avp)^3];

H1perpFavFirstMoment[x_, Q2_] := Sqrt[E/2] (NCfav x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1u["pi+", x, Q2] (MC^3 avp)/(Mh x (MC^2 + avp)^2);
H1perpUnfFirstMoment[x_, Q2_] := Sqrt[E/2] (NCunf x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1d["pi+", x, Q2] (MC^3 avp)/(Mh x (MC^2 + avp)^2);

H1perpFavTMD[x_, Q2_, kt_] := (x Mh)/MC Exp[-kt^2/MC^2] Sqrt[2 E] (NCfav x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1u["pi+", x, Q2] 1/(\[Pi] avp) Exp[-kt^2/avp];
H1perpUnfTMD[x_, Q2_, kt_] := (x Mh)/MC Exp[-kt^2/MC^2] Sqrt[2 E] (NCunf x^alphaC (1 - x)^betaC (alphaC + betaC)^(alphaC + betaC)/((alphaC^alphaC) (betaC^betaC))) D1d["pi+", x, Q2] 1/(\[Pi] avp) Exp[-kt^2/avp];

H1perpFirstMoment[quark_, pion_, x_, Q2_] := If[(quark == "u" && pion == "pi+") || (quark == "bard" && pion == "pi+") || (quark == "d" && pion == "pi-") || (quark == "baru" && pion == "pi-"),H1perpFavFirstMoment[x, Q2], H1perpUnfFirstMoment[x, Q2]]
  
(*2011 fit Appendix A. 3 1107.4446.*)
Clear[Ms,avks,Nu,Nusea,Nd,Ndsea,Nst,Nstbar,alphauv,alphadv,asea,beta, Mp];
Ms=Sqrt[0.19];
avks=avk Ms^2/(avk +Ms^2);
Nu=0.4;
Nd=-0.97;
alphauv=0.35;
alphadv=0.44;
betauv=2.6;
betadv=0.9;
avPTs[z_]:= Sqrt[avp^2+avks^2 z^2]
Mp=0.938;

(* Sivers function with DGLAP Evolution*)
usiv[x_,Q2_]:=(Nu x^alphauv (1-x)^betauv (alphauv+ betauv)^(alphauv+ betauv)/((alphauv^alphauv)( betauv^betauv))) up[x,Q2]
dsiv[x_,Q2_]:=(Nd x^alphadv (1-x)^betadv (alphadv+ betadv)^(alphadv+ betadv)/((alphadv^alphadv)( betadv^betadv)))dn[x,Q2]
ubarsiv[x_,Q2_]:=0.
dbarsiv[x_,Q2_]:=0.
ssiv[x_,Q2_]:=0.
sbarsiv[x_,Q2_]:=0.

 (*Sivers function with DGLAP Evolution*)
f1TperpuFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp Ms) usiv[x,Q2] avks^2/avk ;
f1TperpdFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp Ms) dsiv[x,Q2] avks^2/avk;
f1TperpubarFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp Ms) ubarsiv[x,Q2] avks^2/avk ;
f1TperpdbarFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp Ms) dbarsiv[x,Q2] avks^2/avk;
f1TperpsFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp Ms) ssiv[x,Q2] avks^2/avk ;
f1TperpsbarFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp Ms) sbarsiv[x,Q2] avks^2/avk;


f1TperpuTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]usiv[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpdTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]dsiv[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpubarTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]ubarsiv[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpdbarTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]dbarsiv[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpsTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]ssiv[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avks];
f1TperpsbarTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]sbarsiv[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avks];



(*2015 fit Appendix A. 6 [1411.0580].*)
Clear[MTT,avkTT,NuTT,NdTT,alphaTT,betaTT];
MTT=Sqrt[0.18];
avkTT=avk MTT^2/(avk+MTT^2);
NuTT=1;
NdTT=-1;
alphaTT=2.5;
betaTT=2.;
(* Pretzelosity function with DGLAP Evolution*)
huTT[x_,Q2_]:=E(NuTT x^alphaTT (1-x)^betaTT (alphaTT+ betaTT)^(alphaTT+ betaTT)/((alphaTT^alphaTT)( betaTT^betaTT))) (up[x,Q2]-g1u[x,Q2])
hdTT[x_,Q2_]:=E(NdTT x^alphaTT(1-x)^betaTT (alphaTT+ betaTT)^(alphaTT+ betaTT)/((alphaTT^alphaTT)( betaTT^betaTT)))(dn[x,Q2]-g1d[x,Q2])


 (*Pretzelosity function with DGLAP Evolution*)
h1TperpuFirstMoment[x_,Q2_]:=1/(2 MTT^2) huTT[x,Q2] avkTT^2/avk ;
h1TperpdFirstMoment[x_,Q2_]:=1/(2 MTT^2) hdTT[x,Q2] avkTT^2/avk;



h1TperpuTMD[x_,Q2_,kt_]:=Mp^2/MTT^2  huTT[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avkTT];
h1TperpdTMD[x_,Q2_,kt_]:=Mp^2/MTT^2  hdTT[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avkTT];


h1TperpuSecondMoment[x_,Q2_]:=1/(2 Mp^2 MTT^2) huTT[x,Q2] avkTT^3/avk ;
h1TperpdSecondMoment[x_,Q2_]:=1/(2 Mp^2 MTT^2) hdTT[x,Q2] avkTT^3/avk;


(*2010 fit Barone et al Appendix A. 5 [0912.5194]. https://arxiv.org/pdf/0912.5194.pdf*)
MsBM=Sqrt[0.34];
avksBM=avk MsBM^2/(avk +MsBM^2);
avkBM =avksBM;
NuBM=2.1*0.35;
NuseaBM=-1.*0.04;
NdBM=(-1.111)*(-0.9000);
NdseaBM=-1.*0.4;
NstBM=0;
NstbarBM=0;
alphauvBM=0.73;
alphadvBM=1.08;
aseaBM=0.79;
betaBM=3.46;
avPTsBM[z_]:= Sqrt[avp^2+avksBM^2 z^2] 



 (*Boer-Mulders function with DGLAP Evolution*)
usivBM[x_,Q2_]:=(NuBM x^alphauvBM (1-x)^betaBM (alphauvBM+ betaBM)^(alphauvBM+ betaBM)/((alphauvBM^alphauvBM)( betaBM^betaBM))) up[x,Q2]
dsivBM[x_,Q2_]:=(NdBM x^alphadvBM (1-x)^betaBM (alphadvBM+ betaBM)^(alphadvBM+ betaBM)/((alphadvBM^alphadvBM)( betaBM^betaBM)))dn[x,Q2]
ubarsivBM[x_,Q2_]:=(NuseaBM x^aseaBM (1-x)^betaBM (aseaBM+ betaBM)^(aseaBM+ betaBM)/((aseaBM^aseaBM)( betaBM^betaBM))) upbar[x,Q2]
dbarsivBM[x_,Q2_]:=(NdseaBM x^aseaBM (1-x)^betaBM(aseaBM+ betaBM)^(aseaBM+ betaBM)/((aseaBM^aseaBM)( betaBM^betaBM))) dnbar[x,Q2]
ssivBM[x_,Q2_]:=(NstBM x^aseaBM (1-x)^betaBM (aseaBM+ betaBM)^(aseaBM+ betaBM)/((aseaBM^aseaBM)( betaBM^betaBM))) str[x,Q2]
sbarsivBM[x_,Q2_]:=(NstbarBM x^aseaBM (1-x)^betaBM (aseaBM+ betaBM)^(aseaBM+ betaBM)/((aseaBM^aseaBM)( betaBM^betaBM))) sbar[x,Q2]

 (*Sivers function with DGLAP Evolution*)
h1perpuFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp MsBM) usivBM[x,Q2] avksBM^2/avk ;
h1perpdFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp MsBM) dsivBM[x,Q2] avksBM^2/avk;
h1perpubarFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp MsBM) ubarsivBM[x,Q2] avksBM^2/avk ;
h1perpdbarFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp MsBM) dbarsivBM[x,Q2] avksBM^2/avk;
h1perpsFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp MsBM) ssivBM[x,Q2] avksBM^2/avk ;
h1perpsbFirstMoment[x_,Q2_]:=-Sqrt[ E/2] 1/(Mp MsBM) sbarsivBM[x,Q2] avksBM^2/avk;


h1perpuTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]usivBM[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avksBM];
h1perpdTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]dsivBM[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avksBM];
h1perpubarTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]ubarsivBM[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avksBM];
h1perpdbarTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]dbarsivBM[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avksBM];
h1perpsTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]ssivBM[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avksBM];
h1perpsbTMD[x_,Q2_,kt_]:=-(Mp/Ms) Sqrt[2 E]sbarsivBM[x,Q2]1/(\[Pi] avk) Exp[-kt^2/avksBM];


(*WW-type relations*)
(*use grid to speed up*)
Clear[grid,TableU,TableD,xh1Lp1u,xh1Lp1d,h1Lu,h1Ld];
grid=ReadList["./Grids/xh1Lperp_u_d.dat",Real,RecordLists-> True];
TableU= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,3]]},{i,1,Length[grid]}];
TableD= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,4]]},{i,1,Length[grid]}];
xh1Lp1u=Interpolation[TableU];
xh1Lp1d=Interpolation[TableD];

h1Lu[x_,Q2_]:=xh1Lp1u[x,Q2]/x  ; 
h1Ld[x_,Q2_]:=xh1Lp1d[x,Q2]/x   ; 
(*
h1Lu[x_, Q_] := -x^2 NIntegrate[h1u[y, Q]/y^2, {y, x, 1.}];
h1Ld[x_, Q_] := -x^2  NIntegrate[h1d[y, Q]/y^2, {y, x, 1.}];
*)


Clear[grid,TableU,TableD,TableUBAR,TableDBAR, TableS,TableSBAR];
grid=ReadList["./Grids/gT_u_d_ubar_dbar_s_sbar.dat",Real,RecordLists-> True];
TableU= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,3]]},{i,1,Length[grid]}];
TableD= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,4]]},{i,1,Length[grid]}];
TableUBAR= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,5]]},{i,1,Length[grid]}];
TableDBAR= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,6]]},{i,1,Length[grid]}];
TableS= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,5]]},{i,1,Length[grid]}];
TableSBAR= Table[{{grid[[i,1]],grid[[i,2]]},grid[[i,6]]},{i,1,Length[grid]}];
xgTu=Interpolation[TableU];
xgTd=Interpolation[TableD];
xgTubar=Interpolation[TableUBAR];
xgTdbar=Interpolation[TableDBAR];
xgTs=Interpolation[TableUBAR];
xgTsbar=Interpolation[TableDBAR];

gTu[x_,Q2_]:=xgTu[x,Q2]/x ;  
gTd[x_,Q2_]:=xgTd[x,Q2]/x ;  
 gTubar[x_,Q2_] := xgTubar[x,Q2]/x ;
gTdbar[x_,Q2_] := xgTdbar[x,Q2]/x ;
gTs[x_,Q2_] := xgTs[x,Q2]/x ;
gTsbar[x_,Q2_] := xgTsbar[x,Q2]/x ;
(*Eq. 3.2 (a)*)
(*
gTu[x_,Q2_] := NIntegrate[g1u[y,Q2]/y,{y,x,1.}];
gTd[x_,Q2_] := NIntegrate[g1d[y,Q2]/y,{y,x,1.}];
gTubar[x_,Q2_] := NIntegrate[g1ubar[y,Q2]/y,{y,x,1.}];
gTdbar[x_,Q2_] := NIntegrate[g1dbar[y,Q2]/y,{y,x,1.}];
gTs[x_,Q2_] := NIntegrate[g1s[y,Q2]/y,{y,x,1.}];
gTsbar[x_,Q2_] := NIntegrate[g1sbar[y,Q2]/y,{y,x,1.}];
*)


g1Tperpu[x_,Q2_] := x gTu[x,Q2];
g1Tperpd[x_,Q2_] := x gTd[x,Q2];
g1Tperpubar[x_,Q2_] := x gTubar[x,Q2];
g1Tperpdbar[x_,Q2_] := x gTdbar[x,Q2];
g1Tperps[x_,Q2_] := x gTs[x,Q2];
g1Tperpsbar[x_,Q2_] := x gTsbar[x,Q2];

(*Eq. 3.6 (a)*)
(*
g1Tperpu[x_,Q2_] := x NIntegrate[g1u[y,Q2]/y,{y,x,1.}];
g1Tperpd[x_,Q2_] := x NIntegrate[g1d[y,Q2]/y,{y,x,1.}];
g1Tperpubar[x_,Q2_] := x NIntegrate[g1ubar[y,Q2]/y,{y,x,1.}];
g1Tperpdbar[x_,Q2_] := x NIntegrate[g1dbar[y,Q2]/y,{y,x,1.}];
g1Tperps[x_,Q2_] := x NIntegrate[g1s[y,Q2]/y,{y,x,1.}];
g1Tperpsbar[x_,Q2_] := x NIntegrate[g1sbar[y,Q2]/y,{y,x,1.}];
*)

 



(*Leading Structure Functions and Asymmetries*)
(*F_{UU}*)
avPT[z_] := avp + avk z^2;

FUU[pion_, x_, z_, Q2_, PT_] := ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x,Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sbar[x, Q2] D1sbar[pion, z, Q2]) Exp[(-PT^2/avPT[z])]/(\[Pi] avPT[z])

FUUIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x,Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sbar[x, Q2] D1sbar[pion, z, Q2])
 
(*F_{LL}*)
avLLPT[z_] := avp + avkg z^2;


FLL[pion_, x_, z_, Q2_, PT_] := ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2] ) Exp[(-PT^2/avLLPT[z])]/(\[Pi] avLLPT[z])


FLLIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2])
ALL[pion_,x_,z_,Q2_,PT_] := FLL[pion,x,z,Q2,PT]/FUU[pion,x,z,Q2,PT]  ;

(*F_{UT}^{Sin\phi_h-Sin\phi_S}*)
avUTPTf1tperp[z_] := avp + avks z^2;

GUTf1tperp[z_, PT_] := -(Mp^2/(\[Pi] avUTPTf1tperp[z]^2)) Exp[(-PT^2/avUTPTf1tperp[z])]
AUTf1tperp[z_] := -(( Mp z Sqrt[\[Pi]])/ Sqrt[avUTPTf1tperp[z]])

FUTf1tperp[pion_, x_, z_, Q2_, PT_] := (2 z PT)/Mp ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbarFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) GUTf1tperp[z, PT]


FUTf1tperpIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbarFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) AUTf1tperp[z]

AUTSivers[pion_, x_, z_, Q2_, PT_] := FUTf1tperp[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTSiversIntegrated[pion_, x_, z_, Q2_] := FUTf1tperpIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];

(*F_{UT}^{Sin\phi_h+Sin\phi_S}*)
avUTPTh1[z_] := avp MC^2/(avp + MC^2) + avk z^2;

GUTh1[z_, PT_] := Mp^2/(\[Pi] avUTPTh1[z]^2) Exp[(-PT^2/avUTPTh1[z])]
AUTh1[z_] := ( Mh z Sqrt[\[Pi]])/ Sqrt[avUTPTh1[z]]

FUTh1[pion_, x_, z_, Q2_, PT_] := (2 z PT Mh)/Mp^2 ((4.0/9.0) h1u[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1d[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GUTh1[z, PT]


FUTh1Integrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1u[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1d[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AUTh1[z]

AUTCollins[pion_, x_, z_, Q2_, PT_] := FUTh1[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTCollinsIntegrated[pion_, x_, z_, Q2_] := FUTh1Integrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
 
(*F_{UT}^{Sin[3\phi_h-\phi_S]}*)
avkP = avk MTT^2/(avk + MTT^2);
avUTPTh1tp[z_] := avp MC^2/(avp + MC^2) + avkP z^2;

GUTh1tp[z_, PT_] := (Mp^4 Mh)/(\[Pi] avUTPTh1tp[z]^4) Exp[(-PT^2/avUTPTh1tp[z])]
AUTh1tp[z_] := ( 3 Mh Mp^2 z^3 Sqrt[\[Pi]])/( 2 Sqrt[avUTPTh1tp[z]^3])


FUTh1tp[pion_, x_, z_, Q2_, PT_] := (2 (z^3) (PT^3) )/Mp^2 ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GUTh1tp[z, PT]


FUTh1tpIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AUTh1tp[z]

AUTh1tp[pion_, x_, z_, Q2_, PT_] := FUTh1tp[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTh1tpIntegrated[pion_, x_, z_, Q2_] := FUTh1tpIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UU}^{Cos[2\phi_h]}*)
avc = MC^2 avp/(MC^2 + avp);
avUUcos2phiPT[z_] := avc + avkBM z^2;

GfactorUUcos2phi[z_, PT_] := (4 Mp^3 Mh)/(\[Pi] avUUcos2phiPT[z]^3) Exp[(-PT^2/avUUcos2phiPT[z])]
AfactorUUcos2phi[z_] := ( 4 Mh Mp (z^2)  )/ avUUcos2phiPT[z]

FUUcos2phi[pion_, x_, z_, Q2_, PT_] := (z^2 PT^2)/Mp^2 ((4.0/9.0) h1perpuFirstMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1perpdFirstMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorUUcos2phi[z, PT]


FUUcos2phiIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1perpuFirstMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1perpdFirstMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorUUcos2phi[z]

AUUcos2phi[pion_, x_, z_, Q2_, PT_] := FUUcos2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUUcos2phiIntegrated[pion_, x_, z_, Q2_] := FUUcos2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
 
(*F_{LT}^{Cos[\phi_h-\phi_S]}*)
avkg1T = avkg; (*assumption*)
avLLPT[z_] := avp + avkg1T z^2;

GLT[z_, PT_] := (2 Mp^2)/(\[Pi] avLLPT[z]^2) Exp[(-PT^2/avLLPT[z])]
ALT[z_] := ( Mp z Sqrt[\[Pi]])/ Sqrt[avLLPT[z]]

FLT[pion_, x_, z_, Q2_, PT_] := (z PT)/Mp ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) GLT[z, PT]


FLTIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) ALT[z]

ALT[pion_, x_, z_, Q2_, PT_] := FLT[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALTIntegrated[pion_, x_, z_, Q2_] := FLTIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
 
(*F_{UL}^{Sin[2\phi_h]}*)
avc = MC^2 avp/(MC^2 + avp);
avkh1l = avk; (* assumption *)
avULsin2phiPT[z_] := avc + avkh1l z^2;

GfactorULsin2phi[z_, PT_] := (4 Mp^3 Mh)/(\[Pi] avULsin2phiPT[z]^3) Exp[(-PT^2/avULsin2phiPT[z])]
AfactorULsin2phi[z_] := ( 4 Mh Mp (z^2)  )/ avULsin2phiPT[z]

FULsin2phi[pion_, x_, z_, Q2_, PT_] := (z^2 PT^2)/Mp^2 ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorULsin2phi[z, PT]


FULsin2phiIntegrated[pion_, x_, z_, Q2_] := ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorULsin2phi[z]

AULsin2phi[pion_, x_, z_, Q2_, PT_] := FULsin2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AULsin2phiIntegrated[pion_, x_, z_, Q2_] := FULsin2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*Sub-leading Structure Functions and Asymmetries*)

(*F_{LT}^{Cos[\phi_S]}*)
avLLPT[z_] := avp + avkg z^2;

GLTCOSPhi[z_, PT_] := 1 /(\[Pi] avLLPT[z]) Exp[(-PT^2/avLLPT[z])]
ALTCOSPhi[z_] := 1.

FLTcosPhi[pion_, x_, z_, Q2_, PT_] := (-2 Mp x)/Sqrt[Q2] ((4.0/9.0) gTu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) gTd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) gTs[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) gTubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) gTdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) gTsbar[x, Q2] D1sbar[pion, z, Q2] ) GLTCOSPhi[z, PT];


FLTcosPhiIntegrated[pion_, x_, z_, Q2_] := (-2 Mp x)/Sqrt[Q2] ((4.0/9.0) gTu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) gTd[x,Q2] D1d[pion, z, Q2] + (1.0/9.0) gTs[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) gTubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) gTdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) gTsbar[x, Q2] D1sbar[pion, z, Q2] ) ALTCOSPhi[z];

ALTcosPhi[pion_, x_, z_, Q2_, PT_] := FLTcosPhi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALTcosPhiIntegrated[pion_, x_, z_, Q2_] := FLTcosPhiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UL}^{Sin[\phi_h]}*)
avc = MC^2 avp/(MC^2 + avp);
avkh1l = avk; (* assumption *)
avULsinphiPT[z_] := avc + avkh1l z^2;
(* x hL = h_1L^perp(1) *)

GfactorULsinphi[z_, PT_] := (2 (Mp^2) )/(\[Pi] avULsinphiPT[z]^2) Exp[(-PT^2/avULsinphiPT[z])]
AfactorULsinphi[z_] := ( Mh z Sqrt[Pi]  )/ (avULsinphiPT[z]^(1/2))

FULsinphi[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (z PT Mh)/Mp^2 (-2) ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorULsinphi[z, PT]


FULsinphiIntegrated[pion_, x_, z_, Q2_] := ((2 Mp)/Sqrt[Q2]) (-2) ((4.0/9.0) h1Lu[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1Ld[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorULsinphi[z]

AULsinphi[pion_, x_, z_, Q2_, PT_] := FULsinphi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AULsinphiIntegrated[pion_, x_, z_, Q2_] := FULsinphiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{LL}^{Cos[\phi_h]}*)
avLLcosphiPT[z_] := avp + avkg z^2;

GfactorLLcosphi[z_, PT_] := avkg/(\[Pi] avLLcosphiPT[z]^2) Exp[(-PT^2/avLLcosphiPT[z])]
AfactorLLcosphi[z_] := ( z Sqrt[Pi] avkg)/(2 Mp Sqrt[avLLcosphiPT[z]]) ;

FLLcosphi[pion_, x_, z_, Q2_, PT_] := (-((2 Mp)/Sqrt[Q2])) (z PT)/Mp ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2] ) GfactorLLcosphi[z, PT]


FLLcosphiIntegrated[pion_, x_, z_, Q2_] := (-((2 Mp)/Sqrt[Q2])) ((4.0/9.0) g1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1sbar[x, Q2] D1sbar[pion, z, Q2] ) AfactorLLcosphi[z]

ALLcosphi[pion_, x_, z_, Q2_, PT_] := FLLcosphi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALLcosphiIntegrated[pion_, x_, z_, Q2_] := FLLcosphiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{LT}^{Cos[2\phi_h-\phi_S]}*)
avkg1T = avkg; (*assumption*)
avLTcos2phiPT[z_] := avp + avkg1T z^2;

GLTcos2phi[z_, PT_] := (Mp^2 avkg1T)/(\[Pi] avLTcos2phiPT[z]^3) Exp[(-PT^2/avLTcos2phiPT[z])]
AfactorLTcos2phi[z_] := -(( z^2 avkg1T)/ avLTcos2phiPT[z])

FLTcos2phi[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (-((z^2 PT^2)/Mp^2)) ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) GLTcos2phi[z, PT]


FLTcos2phiIntegrated[pion_, x_, z_, Q2_] := ((2 Mp)/Sqrt[Q2]) ((4.0/9.0) g1Tperpu[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) g1Tperpd[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) g1Tperps[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) g1Tperpubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) g1Tperpdbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) g1Tperpsbar[x, Q2] D1sbar[pion, z, Q2] ) AfactorLTcos2phi[z]

ALTcos2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] := FLTcos2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

ALTcos2PhiminusPhiSIntegrated[pion_, x_, z_, Q2_] := FLTcos2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UU}^{cos\phi_h}*)
avUUcosphi2PT[z_] := avp + avk z^2;
GfactorUUcosphi2[z_, PT_] := 1/(\[Pi] avUUcosphi2PT[z]) Exp[(-PT^2/avUUcosphi2PT[z])]
AfactorUUcosphi2[z_] := (  z Sqrt[Pi] avk)/(2 Mp  Sqrt[avUUcosphi2PT[z]]) ;

FUUcosphi[pion_, x_, z_, Q2_, PT_] := (-((2 Mp)/Sqrt[Q2])) (z PT 2 Mp)/avUUcosphi2PT[z] avk/(2 Mp^2) ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sbar[x, Q2] D1sbar[pion, z, Q2]) GfactorUUcosphi2[z, PT]


FUUcosphiIntegrated[pion_, x_, z_, Q2_] := (-((2 Mp)/Sqrt[Q2])) ((4.0/9.0) f1u[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1d[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1s[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1ubar[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1dbar[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1sbar[x, Q2] D1sbar[pion, z, Q2]) AfactorUUcosphi2[z]

AUUcosphi[pion_, x_, z_, Q2_, PT_] := FUUcosphi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUUcosphiIntegrated[pion_, x_, z_, Q2_] := FUUcosphiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UT}^{Sin[\phi_S]}*)
avkh1 = avk; (*asumption*)
avc = MC^2 avp/(MC^2 + avp);
avUTsinphi2PT[z_] := avc + avkh1 z^2;
GfactorUTsinphi2[z_, PT_] :=  (1 - PT^2/avUTsinphi2PT[z])/(\[Pi] avUTsinphi2PT[z]) Exp[(-PT^2/avUTsinphi2PT[z])]
AfactorUTsinphi2[z_] := 0. ;

FUTsinphiS[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (4 z^2 Mh Mp)/avUTsinphi2PT[z] avkh1/(2 Mp^2) ((4.0/9.0) h1u[x, Q2] If[pion == "pi+", H1perpFavFirstMoment[z, Q2], H1perpUnfFirstMoment[z, Q2]] + (1.0/9.0) h1d[x, Q2] If[pion == "pi+", H1perpUnfFirstMoment[z, Q2], H1perpFavFirstMoment[z, Q2]]) GfactorUTsinphi2[z, PT]


FUTsinphiSIntegrated[pion_, x_, z_, Q2_] := 0

AUTsinphiS[pion_, x_, z_, Q2_, PT_] := FUTsinphiS[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTsinphiSIntegrated[pion_, x_, z_, Q2_] := FUTsinphiSIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(*F_{UT}^{Sin[2\phi_h-\phi_S]}*)
avUTsin2phi1PT[z_] := avp + avks z^2;

GfactorUTsin2phi1[z_, PT_] := (Mp^2)  /(\[Pi] avUTsin2phi1PT[z]) Exp[(-PT^2/avUTsin2phi1PT[z])]
AfactorUTsin2phi1[z_] := (Mp^2 z^2)/ avUTsin2phi1PT[z] ;

avkP = avk MTT^2/(avk + MTT^2);
avc = MC^2 avp/(MC^2 + avp);
avUTsin2phi2PT[z_] := avc + avkP z^2;
GfactorUTsin2phi2[z_, PT_] := 1/(\[Pi] avUTsin2phi2PT[z]) Exp[(-PT^2/avUTsin2phi2PT[z])]
AfactorUTsin2phi2[z_] := (4 Mh Mp z^2)/ avUTsin2phi1PT[z] ;

FUTsin2phi[pion_, x_, z_, Q2_, PT_] := ((2 Mp)/Sqrt[Q2]) (z^2 PT^2)/avUTsin2phi1PT[z]^2 avks/Mp^2 ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbarFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) GfactorUTsin2phi1[z, PT] + ((2 Mp)/Sqrt[Q2]) (z^2 PT^2 4 Mp Mh)/avUTsin2phi2PT[z]^2 (-1) ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) GfactorUTsin2phi2[z, PT]


FUTsin2phiIntegrated[pion_, x_, z_, Q2_] := ((2 Mp)/Sqrt[Q2]) avks/Mp^2 ((4.0/9.0) f1TperpuFirstMoment[x, Q2] D1u[pion, z, Q2] + (1.0/9.0) f1TperpdFirstMoment[x, Q2] D1d[pion, z, Q2] + (1.0/9.0) f1TperpsFirstMoment[x, Q2] D1s[pion, z, Q2]  + (4.0/9.0) f1TperpubarFirstMoment[x, Q2] D1ubar[pion, z, Q2] + (1.0/9.0) f1TperpdbarFirstMoment[x, Q2] D1dbar[pion, z, Q2] + (1.0/9.0) f1TperpsbarFirstMoment[x, Q2] D1sbar[pion, z, Q2] ) AfactorUTsin2phi1[z] + ((2 Mp)/Sqrt[Q2]) (-1) ((4.0/9.0) h1TperpuSecondMoment[x, Q2] H1perpFirstMoment["u", pion, z, Q2] + (1.0/9.0) h1TperpdSecondMoment[x, Q2] H1perpFirstMoment["d", pion, z, Q2]) AfactorUTsin2phi2[z]

AUTsin2PhiminusPhiS[pion_, x_, z_, Q2_, PT_] := FUTsin2phi[pion, x, z, Q2, PT]/FUU[pion, x, z, Q2, PT];

AUTsin2PhiminusPhiSIntegrated[pion_, x_, z_, Q2_] := FUTsin2phiIntegrated[pion, x, z, Q2]/FUUIntegrated[pion, x, z, Q2];
 
(* Calculation of Unpolarized Cross Section *)
p1[y_]:=(1-y)/(1-y+y^2/2);
p2[y_]:=y (1-y/2)/(1-y+y^2/2);
p3[y_]:=(2-y) Sqrt[1-y]/(1-y+y^2/2);
p4[y_]:=y Sqrt[1-y]/(1-y+y^2/2);
coupling=1./137.;
y[x_,Q2_,energy_]:=Q2/(2.*Mp*energy*x);
msg::nnarg="Invalid input.";
 
(* Cross Section of unpolarized beam and unpolarized target *)
CrossSectionUU[pion_,x_,z_,Q2_,PT_,energy_,phi_]:=coupling^2/(x y[x,Q2,energy] Q2) (1-y[x,Q2,energy]+y[x,Q2,energy]^2/2)* x *(FUU[pion,x,z,Q2,PT]+Cos[phi] p3[y[x,Q2,energy]] FUUcosphi[pion,x,z,Q2,PT]+Cos[2 phi] p1[y[x,Q2,energy]] FUUcos2phi[pion,x,z,Q2,PT])
 
(* Cross Section of unpolarized beam and longitidunally polarized target *)
CrossSectionUL[pion_,x_,z_,Q2_,PT_,energy_,phi_]:=CrossSectionUU[pion,x,z,Q2,PT,energy,phi]+coupling^2/(x y[x,Q2,energy] Q2) (1-y[x,Q2,energy]+y[x,Q2,energy]^2/2) * x *(Sin[2 phi] p1[y[x,Q2,energy]] FULsin2phi[pion, x, z, Q2, PT]+Sin[phi] p3[y[x,Q2,energy]] FULsinphi[pion, x, z, Q2, PT])
                                                                                                                                        
(* Cross Section of longitidunally polarized beam and unpolarized target *)
CrossSectionLU[pion_,x_,z_,Q2_,PT_,energy_,phi_,helicity_]:=CrossSectionUU[pion,x,z,Q2,PT,energy,phi]+0.0
 
(* Cross Section of longitidunally polarized beam and longitidunally polarized target *)
CrossSectionLL[pion_,x_,z_,Q2_,PT_,energy_,phi_,helicity_]:=CrossSectionUU[pion,x,z,Q2,PT,energy,phi]+coupling^2/(x y[x,Q2,energy] Q2) (1-y[x,Q2,energy]+y[x,Q2,energy]^2/2) * x *(helicity p2[y[x,Q2,energy]] FLL[pion, x, z, Q2, PT] + helicity Cos[phi] p4[y[x,Q2,energy]] FLLcosphi[pion, x, z, Q2, PT])
 
(* Cross Section of unpolarized beam and transversely polarized target *)
CrossSectionUT[pion_,x_,z_,Q2_,PT_,energy_,phih_,phiS_]:=CrossSectionUU[pion,x,z,Q2,PT,energy,phih]+coupling^2/(x y[x,Q2,energy] Q2) (1-y[x,Q2,energy]+y[x,Q2,energy]^2/2) * x *(Sin[phih - phiS] FUTf1tperp[pion, x, z, Q2, PT] + Sin[phih + phiS] p1[y[x,Q2,energy]] FUTh1[pion, x, z, Q2, PT] + Sin[3*phih - phiS] p1[y[x,Q2,energy]] FUTh1tp[pion, x, z, Q2, PT] + Sin[phiS] p3[y[x,Q2,energy]] FUTsinphiS[pion, x, z, Q2, PT] + Sin[2*phih - phiS] p3[y[x,Q2,energy]] FUTsin2phi[pion, x, z, Q2, PT])
 
(* Cross Section of longitudinally polarized beam and transversely polarized target *)
CrossSectionLT[pion_,x_,z_,Q2_,PT_,energy_,phih_,phiS_,helicity_]:=CrossSectionUU[pion,x,z,Q2,PT,energy,phih]+coupling^2/(x y[x,Q2,energy] Q2) (1-y[x,Q2,energy]+y[x,Q2,energy]^2/2) * x *(helicity Cos[phih - phiS] p2[y[x,Q2,energy]] FLT[pion, x, z, Q2, PT] + helicity Cos[phiS] p4[y[x,Q2,energy]] FLTcosPhi[pion, x, z, Q2, PT] + helicity Cos[2*phih - phiS] p4[y[x,Q2,energy]] FLTcos2phi[pion, x, z, Q2, PT])
 
CrossSection[pion_,x_,z_,Q2_,PT_,energy_,phih_,phiS_,helicity_,targetpolarization_]:=Which[helicity==0 && targetpolarization=="U", CrossSectionUU[pion,x,z,Q2,PT,energy,phih],helicity==0 && targetpolarization=="L", CrossSectionUL[pion,x,z,Q2,PT,energy,phih], (helicity==1 || helicity==-1) && targetpolarization=="U", 0.0, (helicity==1 || helicity==-1) && targetpolarization=="L", CrossSectionLL[pion,x,z,Q2,PT,energy,phih,helicity], helicity==0 && targetpolarization=="T", CrossSectionUT[pion,x,z,Q2,PT,energy,phih,phiS], (helicity==1 || helicity==-1) && targetpolarization=="T", CrossSectionLT[pion,x,z,Q2,PT,energy,phih,phiS,helicity], True, Message[msg::nnarg]]

End[];
EndPackage[];
