(* ::Package:: *)

Print["Mathematica package for MSTW PDFs\nby Graeme Watt <Graeme.Watt(at)cern.ch>.\nFor information on usage, see ?ReadPDFGrid and ?xf."];
BeginPackage["mstwpdf`"]
ReadPDFGrid::usage="The function ReadPDFGrid[prefix,ih] reads the PDF grid corresponding to eigenvector set ih, where ih=0 is the central set, into memory from the location specified by prefix."
xf::usage="The function xf[ih,x,q,f] returns x times the parton distribution of flavour f corresponding to eigenvector set ih, where ih=0 is the central set, at a given momentum fraction x and scale q in GeV.  Note that ih and f must be integers, and x and q must be numeric quantities.  The PDG convention is used for the flavour f (apart from the gluon has f=0, not 21), that is, f = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 corresponds to tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t.  The valence quark distributions can also be obtained directly with f = 7, 8, 9, 10, 11, 12 corresponding to dv, uv, sv, cv, bv, tv.  The photon distribution is obtained with f = 13."
nhess::usage="The number of eigenvector PDF sets."
distance::usage="Distance along a particular eigenvector direction (t)."
tolerance::usage="Square root of the change in the global chi-squared (T) obtained by moving a distance t along a particular eigenvector direction.  Under ideal quadratic behaviour, t = T."
mCharm::usage="Charm quark mass in GeV."
mBottom::usage="Bottom quark mass in GeV."
alphaSQ0::usage="Value of alphaS at Q0 = 1 GeV."
alphaSMZ::usage="Value of alphaS at the Z mass."
alphaSorder::usage="Perturbative order of alphaS:  0=LO, 1=NLO, 2=NNLO."
alphaSnfmax::usage="Maximum number of flavours in alphaS evolution."
Begin["`Private`"]
xmin=1.*^-6;xmax=1.0;qsqmin=1.0;qsqmax=1.*^9;eps=1.*^-6;
np=12;nx=64;nq=48;nhess=2*20;nqc0=0;nqb0=0;
xx={1.*^-6,2.*^-6,4.*^-6,6.*^-6,8.*^-6,1.*^-5,2.*^-5,4.*^-5,6.*^-5,8.*^-5, 1.*^-4,2.*^-4,4.*^-4,6.*^-4,8.*^-4, 1.*^-3,2.*^-3,4.*^-3,6.*^-3,8.*^-3,1.*^-2,1.4*^-2,2.*^-2,3.*^-2,4.*^-2,6.*^-2,8.*^-2,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,.325,.35,.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0};
qq={1.*^0,1.25*^0,1.5*^0,0.,0.,2.5*^0,3.2*^0,4.*^0,5.*^0,6.4*^0,8.*^0,1.*^1,1.2*^1,0.,0.,2.6*^1,4.*^1,6.4*^1,1.*^2,1.6*^2,2.4*^2,4.*^2,6.4*^2,1.*^3,1.8*^3,3.2*^3,5.6*^3,1.*^4,1.8*^4,3.2*^4,5.6*^4,1.*^5,1.8*^5,3.2*^5,5.6*^5,1.*^6,1.8*^6,3.2*^6,5.6*^6,1.*^7,1.8*^7,3.2*^7,5.6*^7,1.*^8,1.8*^8,3.2*^8,5.6*^8,1.*^9};
qqorig=qq;
cc=Table[0.,{ip,1,np},{ih,1,nhess+1},{ix,1,nx},{iq,1,nq},{4},{4}];
ReadPDFGrid[prefix_String,ih_Integer]:=Module[{pdffile,pdfdat,mc2,mb2},
pdffile=prefix<>"."<>If[ih<10,"0",""]<>ToString[ih]<>".dat";
pdfdat=OpenRead[pdffile];
ReadList[pdfdat,String,2];
ReadList[pdfdat,Word,3];distance=Read[pdfdat,Number];tolerance=Read[pdfdat,Number];
ReadList[pdfdat,Word,3];mCharm=Read[pdfdat,Number];
ReadList[pdfdat,Word,3];mBottom=Read[pdfdat,Number];
ReadList[pdfdat,Word,3];alphaSQ0=Read[pdfdat,Number];
ReadList[pdfdat,Word,3];alphaSMZ=Read[pdfdat,Number];
ReadList[pdfdat,Word,3];alphaSorder=Read[pdfdat,Number];comma=Read[pdfdat,Character];alphaSnfmax=Read[pdfdat,Number];
ReadList[pdfdat,Word,3];nExtraFlavours=Read[pdfdat,Number];
(* Next line ensures compatibility with both Mathematica 5 and 6. *)
While[StringLength[ReadList[pdfdat,String,1]][[1]]==0,Null];
mc2=mCharm^2;
mb2=mBottom^2;
qq=qqorig;
Which[mc2<=qq[[1]]||mc2+eps>=qq[[8]],{Print["Error in ReadPDFGrid: invalid mCharm = ",mCharm];Abort[];},mc2<qq[[2]],{nqc0=2;qq[[4]]=qq[[2]];qq[[5]]=qq[[3]];},mc2<qq[[3]],{nqc0=3;qq[[5]]=qq[[3]];},mc2<qq[[6]],{nqc0=4;},mc2<qq[[7]],{nqc0=5;qq[[4]]=qq[[6]];},True,{nqc0=6;qq[[4]]=qq[[6]];qq[[5]]=qq[[7]];}];
Which[mb2<=qq[[12]]||mb2+eps>=qq[[17]],{Print["Error in ReadPDFGrid: invalid mBottom = ",mBottom];Abort[];},mb2<qq[[13]],{nqb0=13;qq[[15]]=qq[[13]];},mb2<qq[[16]],{nqb0=14;},True,{nqb0=15;qq[[14]]=qq[[16]];}];
qq[[nqc0]]=mc2;
qq[[nqc0+1]]=mc2+eps;
qq[[nqb0]]=mb2;
qq[[nqb0+1]]=mb2+eps;
If[nExtraFlavours<0||nExtraFlavours>1,{Print["Error in ReadPDFGrid: invalid nExtraFlavours = ",nExtraFlavours];Abort[];}];
xxlog=Log[10,xx];
qqlog=Log[10,qq];
grid=Table[0.,{ip,1,np},{ixq,1,nx*nq}];
Do[grid[[ip,(ix-1)*nq+iq]]={{xxlog[[ix]],qqlog[[iq]]},If[ix==nx||(alphaSorder<2&&(ip==10||ip==11))||(nExtraFlavours<1&&ip==12),0.,Read[pdfdat,Number]]},{ix,1,nx},{iq,1,nq},{ip,1,np}];
If[Read[pdfdat]==EndOfFile,Null,Null,{Print["Error in ReadPDFGrid: not at end of ",pdffile];Abort[];}];
Close[pdfdat];ff=Table[0.,{ip,1,np},{ix,1,nx},{iq,1,nq}];
Do[ff[[ip,ix,iq]]=grid[[ip,(ix-1)*nq+iq,2]],{ip,1,np},{ix,1,nx},{iq,1,nq}];
Do[InitialisePDF[p,ih],{p,1,np}];
Print["PDF grid read from ",pdffile];
]
polderiv1[x1_,x2_,x3_,y1_,y2_,y3_]:=(x3*x3*(y1-y2)+2.*x1*(x3*(-y1+y2)+x2*(y1-y3))+x2*x2*(-y1+y3)+x1*x1*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3))
polderiv2[x1_,x2_,x3_,y1_,y2_,y3_]:=
(x3*x3*(y1-y2)-2.*x2*(x3*(y1-y2)+x1*(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
polderiv3[x1_,x2_,x3_,y1_,y2_,y3_]:=(x3*x3*(-y1+y2)+2.*x2*x3*(y1-y3)+x1*x1*(y2-y3)+x2*x2*(-y1+y3)+2.*x1*x3*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3))
InitialisePDF[ip_Integer,ih_Integer]:=Module[{j,k,l,n,m,ff1,ff2,ff12,ff21,yy0,yy1,yy2,yy12,z,cl,iwt,d1,d2,d1d2,xxd},
ff1=Table[0.,{ix,1,nx},{iq,1,nq}];ff2=Table[0.,{ix,1,nx},{iq,1,nq}];
ff12=Table[0.,{ix,1,nx},{iq,1,nq}];ff21=Table[0.,{ix,1,nx},{iq,1,nq}];
yy0=Table[0.,{4}];yy1=Table[0.,{4}];yy2=Table[0.,{4}];yy12=Table[0.,{4}];
z=Table[0.,{16}];cl=Table[0.,{16}];
iwt={{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{ 0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},{-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},{2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},{ 0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},{ 0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},{-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},{9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
{-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2}, {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},{-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},{4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};
Do[{
ff1[[1,m]]=polderiv1[xxlog[[1]],xxlog[[2]],xxlog[[3]],ff[[ip,1,m]],ff[[ip,2,m]],ff[[ip,3,m]]];
ff1[[nx,m]]=polderiv3[xxlog[[nx-2]],xxlog[[nx-1]],xxlog[[nx]],ff[[ip,nx-2,m]],ff[[ip,nx-1,m]],ff[[ip,nx,m]]];
Do[ff1[[n,m]]=polderiv2[xxlog[[n-1]],xxlog[[n]],xxlog[[n+1]],ff[[ip,n-1,m]],ff[[ip,n,m]],ff[[ip,n+1,m]]],{n,2,nx-1}];
},{m,1,nq}];
Do[
Which[
nqc0==2&&m==1,
ff2[[n,m]]=(ff[[ip,n,m+1]]-ff[[ip,n,m]])/(qqlog[[m+1]]-qqlog[[m]]),
nqc0==2&&m==2,
ff2[[n,m]]=(ff[[ip,n,m]]-ff[[ip,n,m-1]])/(qqlog[[m]]-qqlog[[m-1]]),
m==1||m==nqc0+1||m==nqb0+1,
ff2[[n,m]]=polderiv1[qqlog[[m]],qqlog[[m+1]],qqlog[[m+2]],ff[[ip,n,m]],ff[[ip,n,m+1]],ff[[ip,n,m+2]]],
m==nq||m==nqc0||m==nqb0,
ff2[[n,m]]=polderiv3[qqlog[[m-2]],qqlog[[m-1]],qqlog[[m]],ff[[ip,n,m-2]],ff[[ip,n,m-1]],ff[[ip,n,m]]],
True,
ff2[[n,m]]=polderiv2[qqlog[[m-1]],qqlog[[m]],qqlog[[m+1]],ff[[ip,n,m-1]],ff[[ip,n,m]],ff[[ip,n,m+1]]]
],{n,1,nx},{m,1,nq}];
Do[{ff12[[1,m]]=polderiv1[xxlog[[1]],xxlog[[2]],xxlog[[3]],ff2[[1,m]],ff2[[2,m]],ff2[[3,m]]];
ff12[[nx,m]]=polderiv3[xxlog[[nx-2]],xxlog[[nx-1]],xxlog[[nx]],ff2[[nx-2,m]],ff2[[nx-1,m]],ff2[[nx,m]]];
Do[ff12[[n,m]]=polderiv2[xxlog[[n-1]],xxlog[[n]],xxlog[[n+1]],ff2[[n-1,m]],ff2[[n,m]],ff2[[n+1,m]]],{n,2,nx-1}]
},{m,1,nq}];
Do[
Which[
nqc0==2&&m==1,
ff21[[n,m]]=(ff1[[n,m+1]]-ff1[[n,m]])/(qqlog[[m+1]]-qqlog[[m]]),
nqc0==2&&m==2,
ff21[[n,m]]=(ff1[[n,m]]-ff1[[n,m-1]])/(qqlog[[m]]-qqlog[[m-1]]),
m==1||m==nqc0+1||m==nqb0+1,
ff21[[n,m]]=polderiv1[qqlog[[m]],qqlog[[m+1]],qqlog[[m+2]],ff1[[n,m]],ff1[[n,m+1]],ff1[[n,m+2]]],
m==nq||m==nqc0||m==nqb0,
ff21[[n,m]]=polderiv3[qqlog[[m-2]],qqlog[[m-1]],qqlog[[m]],ff1[[n,m-2]],ff1[[n,m-1]],ff1[[n,m]]],
True,
ff21[[n,m]]=polderiv2[qqlog[[m-1]],qqlog[[m]],qqlog[[m+1]],ff1[[n,m-1]],ff1[[n,m]],ff1[[n,m+1]]]
],{n,1,nx},{m,1,nq}];
Do[ff12[[n,m]]=0.5*(ff12[[n,m]]+ff21[[n,m]]),{n,1,nx},{m,1,nq}];
Do[
{
d1=xxlog[[n+1]]-xxlog[[n]];
d2=qqlog[[m+1]]-qqlog[[m]];
d1d2=d1*d2;
yy0[[1]]=ff[[ip,n,m]];
yy0[[2]]=ff[[ip,n+1,m]];
yy0[[3]]=ff[[ip,n+1,m+1]];
yy0[[4]]=ff[[ip,n,m+1]];
yy1[[1]]=ff1[[n,m]];
yy1[[2]]=ff1[[n+1,m]];
yy1[[3]]=ff1[[n+1,m+1]];
yy1[[4]]=ff1[[n,m+1]];
yy2[[1]]=ff2[[n,m]];
yy2[[2]]=ff2[[n+1,m]];
yy2[[3]]=ff2[[n+1,m+1]];
yy2[[4]]=ff2[[n,m+1]];
yy12[[1]]=ff12[[n,m]];
yy12[[2]]=ff12[[n+1,m]];
yy12[[3]]=ff12[[n+1,m+1]];
yy12[[4]]=ff12[[n,m+1]];
Do[{
z[[k]]=yy0[[k]];
z[[k+4]]=yy1[[k]]*d1;
z[[k+8]]=yy2[[k]]*d2;
z[[k+12]]=yy12[[k]]*d1d2;
},{k,1,4}];
(* Equivalent of Fortran code: vector = vector * matrix. *)
(*
Do[{
xxd=0.;
Do[xxd+=iwt[[l,k]]*z[[k]],{k,1,16}];
cl[[l]]=xxd;
},{l,1,16}];
*)
cl=iwt.z; (* Faster to use Mathematica routine instead. *)
l=0;
Do[{
l++;cc[[ip,ih+1,n,m,k,j]]=cl[[l]];
},{k,1,4},{j,1,4}];
},{n,1,nx-1},{m,1,nq-1}];
]
locx[x_Real]:=Module[{ju,jl,jm},
If[x==xxlog[[1]],Return[1]];
If[x==xxlog[[nx]],Return[nx-1]];
ju=nx+1;jl=0;
While[(ju-jl)>1,{jm=IntegerPart[0.5*(ju+jl)];If[x>=xxlog[[jm]],jl=jm,ju=jm]}];
Return[jl];
]
locq[y_Real]:=Module[{ju,jl,jm},
If[y==qqlog[[1]],Return[1]];
If[y==qqlog[[nq]],Return[nq-1]];
ju=nq+1;jl=0;
While[(ju-jl)>1,{jm=IntegerPart[0.5*(ju+jl)];If[y>=qqlog[[jm]],jl=jm,ju=jm]}];
Return[jl];
]
InterpolatePDF[ip_Integer,ih_Integer,x_Real,y_Real]:=
Module[{n,m,t,u,z},
n=locx[x];
m=locq[y];
t=(x-xxlog[[n]])/(xxlog[[n+1]]-xxlog[[n]]);
u=(y-qqlog[[m]])/(qqlog[[m+1]]-qqlog[[m]]);
z=0.;
Do[z=t*z+((cc[[ip,ih+1,n,m,l,4]]*u+cc[[ip,ih+1,n,m,l,3]])*u+cc[[ip,ih+1,n,m,l,2]])*u+cc[[ip,ih+1,n,m,l,1]],{l,4,1,-1}];
Return[z];
]
ExtrapolatePDF[ip_Integer,ih_Integer,x_Real,y_Real]:=
Module[{n,m,f0,f1,z,z0,z1},
n=locx[x];
m=locq[y];
Which[
n==0&&m>0&&m<nq,{
f0=InterpolatePDF[ip,ih,xxlog[[1]],y];
f1=InterpolatePDF[ip,ih,xxlog[[2]],y];
If[f0>=1.*^-3&&f1>=1.*^-3,z=Exp[Log[f0]+(Log[f1]-Log[f0])/(xxlog[[2]]-xxlog[[1]])*(x-xxlog[[1]])],z=f0+(f1-f0)/(xxlog[[2]]-xxlog[[1]])*(x-xxlog[[1]])];
},
n>0&&m==nq,{
f0=InterpolatePDF[ip,ih,x,qqlog[[nq]]];
f1=InterpolatePDF[ip,ih,x,qqlog[[nq-1]]];
If[f0>=1.*^-3&&f1>=1.*^-3,z=Exp[Log[f0]+(Log[f0]-Log[f1])/(qqlog[[nq]]-qqlog[[nq-1]])*(y-qqlog[[nq]])],
z=f0+(f0-f1)/(qqlog[[nq]]-qqlog[[nq-1]])*(y-qqlog[[nq]])];
},
n==0&&m==nq,{
f0=InterpolatePDF[ip,ih,xxlog[[1]],qqlog[[nq]]];
f1=InterpolatePDF[ip,ih,xxlog[[1]],qqlog[[nq-1]]];
If[f0>=1.*^-3&&f1>=1.*^-3,
z0=Exp[Log[f0]+(Log[f0]-Log[f1])/(qqlog[[nq]]-qqlog[[nq-1]])*(y-qqlog[[nq]])],
z0=f0+(f0-f1)/(qqlog[[nq]]-qqlog[[nq-1]])*(y-qqlog[[nq]])];
f0=InterpolatePDF[ip,ih,xxlog[[2]],qqlog[[nq]]];
f1=InterpolatePDF[ip,ih,xxlog[[2]],qqlog[[nq-1]]];
If[f0>=1.*^-3&&f1>=1.*^-3,
z1=Exp[Log[f0]+(Log[f0]-Log[f1])/(qqlog[[nq]]-qqlog[[nq-1]])*(y-qqlog[[nq]])],
z1=f0+(f0-f1)/(qqlog[[nq]]-qqlog[[nq-1]])*(y-qqlog[[nq]])];
If[z0>=1.*^-3&&z1>=1.*^-3,
z=Exp[Log[z0]+(Log[z1]-Log[z0])/(xxlog[[2]]-xxlog[[1]])*(x-xxlog[[1]])],
z=z0+(z1-z0)/(xxlog[[2]]-xxlog[[1]])*(x-xxlog[[1]])];
},
True,{
Print["Error in ExtrapolatePDF"];
Abort[];
}
];
Return[z];
]
xf[ih_Integer,xin_?NumericQ,qin_?NumericQ,f_Integer]:=Module[{qsq,xlog,qlog,res=0.,res1,anom,ip},
If[Im[xin]!=0||Im[qin]!=0,{Print["Error in xf: invalid (x,q) = (",xin,",",qin,")"];Abort[];}];
x=N[xin];q=N[qin];qsq=q*q;
If[Im[x]!=0||Im[q]!=0,{}];
If[qsq>qq[[nqc0]]&&qsq<qq[[nqc0+1]],qsq=qq[[nqc0+1]]];
If[qsq>qq[[nqb0]]&&qsq<qq[[nqb0+1]],qsq=qq[[nqb0+1]]];
xlog=Log[10,x];qlog=Log[10,qsq];
If[ih<0||ih>nhess,{Print["Error in xf: invalid ih = ",ih];Abort[];}];
Which[f==0,ip=1,1<=f<=5,ip=f+1,-5<=f<=-1,ip=-f+1,7<=f<=11,ip=f,f==13,ip=12,Abs[f]!=6&&f!=12,{Print["Error in xf: invalid f = ",f];Abort[];}];
Which[x<=0.||x>xmax||q<=0.,{Print["Error in xf: invalid (x,q) = (",x,",",q,")"];Abort[];},Abs[f]==6||f==12,res=0.,qsq<qsqmin,{If[x<xmin,{res=ExtrapolatePDF[ip,ih,xlog,Log[10,qsqmin]];res1=ExtrapolatePDF[ip,ih,xlog,Log[10,1.01*qsqmin]];If[-5<=f<=-1,{res-=ExtrapolatePDF[ip+5,ih,xlog,Log[10,qsqmin]];res1-=ExtrapolatePDF[ip+5,ih,xlog,Log[10,1.01*qsqmin]];}];},{res=InterpolatePDF[ip,ih,xlog,Log[10,qsqmin]];res1=InterpolatePDF[ip,ih,xlog,Log[10,1.01*qsqmin]];If[-5<=f<=-1,{res-=InterpolatePDF[ip+5,ih,xlog,Log[10,qsqmin]];res1-=InterpolatePDF[ip+5,ih,xlog,Log[10,1.01*qsqmin]];}];}];anom=If[Abs[res]>=1.*^-5,Max[-2.5,(res1-res)/res/0.01],1.];res*=(qsq/qsqmin)^(anom*qsq/qsqmin+1.-qsq/qsqmin);},x<xmin||qsq>qsqmax,{res=ExtrapolatePDF[ip,ih,xlog,qlog];If[-5<=f<=-1,res-=ExtrapolatePDF[ip+5,ih,xlog,qlog]];},Abs[f]!=6&&f!=12,{res=InterpolatePDF[ip,ih,xlog,qlog];If[-5<=f<=-1,res-=InterpolatePDF[ip+5,ih,xlog,qlog]];}];
Return[res];
]
End[]
EndPackage[]
