/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "SwitchingFunction.h"
#include "Tools.h"
#include "Keywords.h"
#include <vector>
#include <limits>

#ifdef __PLUMED_HAS_MATHEVAL
#include <matheval.h>
#endif


using namespace std;
namespace PLMD{

//+PLUMEDOC INTERNAL switchingfunction 
/*
Functions that measure whether values are less than a certain quantity.

Switching functions \f$s(r)\f$ take a minimum of one input parameter \f$d_0\f$.
For \f$r \le d_0 \quad s(r)=1.0\f$ while for \f$r > d_0\f$ the function decays smoothly to 0.
The various switching functions available in plumed differ in terms of how this decay is performed.

Where there is an accepted convention in the literature (e.g. \ref COORDINATION) on the form of the 
switching function we use the convention as the default.  However, the flexibility to use different
switching functions is always present generally through a single keyword. This keyword generally 
takes an input with the following form:

\verbatim
KEYWORD={TYPE <list of parameters>}
\endverbatim  

The following table contains a list of the various switching functions that are available in plumed 2
together with an example input.

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> 
<td> TYPE </td> <td> FUNCTION </td> <td> EXAMPLE INPUT </td> <td> DEFAULT PARAMETERS </td>
</tr> <tr> <td>RATIONAL </td> <td>
\f$
s(r)=\frac{ 1 - \left(\frac{ r - d_0 }{ r_0 }\right)^{n} }{ 1 - \left(\frac{ r - d_0 }{ r_0 }\right)^{m} } 
\f$
</td> <td>
{RATIONAL R_0=\f$r_0\f$ D_0=\f$d_0\f$ NN=\f$n\f$ MM=\f$m\f$}
</td> <td> \f$d_0=0.0\f$, \f$n=6\f$, \f$m=2n\f$ </td>
</tr> <tr>
<td> EXP </td> <td>
\f$
s(r)=\exp\left(-\frac{ r - d_0 }{ r_0 }\right)
\f$
</td> <td> 
{EXP  R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> \f$d_0=0.0\f$ </td>
</tr> <tr>
<td> GAUSSIAN </td> <td>
\f$
s(r)=\exp\left(-\frac{ (r - d_0)^2 }{ 2r_0^2 }\right)
\f$
</td> <td>
{GAUSSIAN R_0=\f$r_0\f$ D_0=\f$d_0\f$} 
</td> <td> \f$d_0=0.0\f$ </td>
</tr> <tr> 
<td> SMAP </td> <td>
\f$
s(r) = \left[ 1 + ( 2^{a/b} -1 )\left( \frac{r-d_0}{r_0} \right)^a \right]^{-b/a}
\f$
</td> <td>
{SMAP R_0=\f$r_0\f$ D_0=\f$d_0\f$ A=\f$a\f$ B=\f$b\f$}
</td> <td> \f$d_0=0.0\f$ </td>
</tr> <tr> 
<td> Q </td> <td>
\f$
s(r) = \frac{1}{1 + \exp(\beta(r_{ij} - \lambda r_{ij}^0))}
\f$
</td> <td>
{Q REF=\f$r_{ij}^0\f$ BETA=\f$\beta\f$ LAMBDA=\f$\lambda\f$ }
</td> <td> \f$\lambda=1.8\f$,  \f$\beta=50 nm^-1\f$ (all-atom)<br/>\f$\lambda=1.5\f$,  \f$\beta=50 nm^-1\f$ (coarse-grained)  </td>
</tr> <tr> 
<td> CUBIC </td> <td>
\f$
s(r) = (y-1)^2(1+2y) \qquad \textrm{where} \quad y = \frac{r - r_1}{r_0-r_1}
\f$
</td> <td>
{CUBIC D_0=\f$r_1\f$ D_MAX=\f$r_0\f$}
</td> <td> </td>
</tr> <tr> 
<td> TANH </td> <td>
\f$
s(r) = 1 - \tanh\left( \frac{ r - d_0 }{ r_0 } \right) 
\f$
</td> <td>
{TANH R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> </td>
</tr> <tr> 
<td> MATHEVAL </td> <td>
\f$
s(r) = FUNC
\f$
</td> <td>
{MATHEVAL FUNC=1/(1+x^6) R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> </td>
</tr>
</table>

\attention
Similarly to the \ref MATHEVAL function, the MATHEVAL switching function 
only works if libmatheval is installed on the system and
PLUMED has been linked to it
Also notice that using MATHEVAL is much slower than using e.g. RATIONAL.
Thus, the MATHEVAL switching function is useful to perform quick
tests on switching functions with arbitrary form before proceeding to their
implementation in C++.

For all the switching functions in the above table one can also specify a further (optional) parameter using the parameter
keyword D_MAX to assert that for \f$r>d_{\textrm{max}}\f$ the switching function can be assumed equal to zero. 
In this case the function is brought smoothly to zero by stretching and shifting it.
\verbatim
KEYWORD={RATIONAL R_0=1 D_MAX=3}
\endverbatim
the resulting switching function will be
\f$
s(r) = \frac{s'(r)-s'(d_{max})}{s'(0)-s'(d_{max})}
\f$
where
\f$
s'(r)=\frac{1-r^6}{1-r^{12}}
\f$
Since PLUMED 2.2 this is the default. The old behavior (no stretching) can be obtained with the
NOSTRETCH flag. The NOSTRETCH keyword is only provided for backward compatibility and might be
removed in the future. Similarly, the STRETCH keyword is still allowed but has no effect.

Notice that switching functions defined with the simplified syntax are never stretched
for backward compatibility. This might change in the future.

*/
//+ENDPLUMEDOC

void SwitchingFunction::registerKeywords( Keywords& keys ){
  keys.add("compulsory","R_0","the value of R_0 in the switching function");
  keys.add("compulsory","D_0","0.0","the value of D_0 in the switching function");
  keys.add("optional","D_MAX","the value at which the switching function can be assumed equal to zero");
  keys.add("compulsory","NN","6","the value of n in the switching function (only needed for TYPE=RATIONAL)");
  keys.add("compulsory","MM","0","the value of m in the switching function (only needed for TYPE=RATIONAL); 0 implies 2*NN");
  keys.add("compulsory","A","the value of a in the switching funciton (only needed for TYPE=SMAP)");
  keys.add("compulsory","B","the value of b in the switching funciton (only needed for TYPE=SMAP)"); 
  keys.add("compulsory","RLOWER","the lower cutoff in the switching function (only needed for MODEXP)");
  keys.add("compulsory","RUPPER","the upper cutoff in the switching function (only needed for MODEXP)");
  keys.addFlag("NOTRUNC",false,"switching function is still computed every distance <= 0.0");
}

void SwitchingFunction::set(const std::string & definition,std::string& errormsg){
  vector<string> data=Tools::getWords(definition);
  if( data.size()<1 ) errormsg="missing all input for switching function"; 
  string name=data[0];
  data.erase(data.begin());
  invr0=0.0;
  invr0_2=0.0;
  d0=0.0;
  dmax=std::numeric_limits<double>::max();
  dmax_2=std::numeric_limits<double>::max();
  stretch=1.0;
  shift=0.0;
  init=true;

  bool present;

  present=Tools::findKeyword(data,"D_0");
  if(present && !Tools::parse(data,"D_0",d0)) errormsg="could not parse D_0";

  present=Tools::findKeyword(data,"D_MAX");
  if(present && !Tools::parse(data,"D_MAX",dmax)) errormsg="could not parse D_MAX";
  dmax_2=dmax*dmax;
  bool dostretch=false;
  Tools::parseFlag(data,"STRETCH",dostretch); // this is ignored now
  dostretch=true;
  bool dontstretch=false;
  Tools::parseFlag(data,"NOSTRETCH",dontstretch); // this is ignored now
  if(dontstretch) dostretch=false;
  do_trunc = true;
  bool notrunc; Tools::parseFlag(data,"NOTRUNC",notrunc);
  if(notrunc) do_trunc = false;
  double r0;
  if(name=="CUBIC"){
     r0 = dmax - d0;
  } else {
     bool found_r0=Tools::parse(data,"R_0",r0);
     if(!found_r0) errormsg="R_0 is required";
  }
  invr0=1.0/r0;
  invr0_2=invr0*invr0;

  if(name=="RATIONAL"){
    type=rational;
    nn=6;
    mm=0;
    present=Tools::findKeyword(data,"NN");
    if(present && !Tools::parse(data,"NN",nn)) errormsg="could not parse NN";
    present=Tools::findKeyword(data,"MM");
    if(present && !Tools::parse(data,"MM",mm)) errormsg="could not parse MM";
    if(mm==0) mm=2*nn;
  } else if(name=="SMAP"){
    type=smap;
    present=Tools::findKeyword(data,"A");
    if(present && !Tools::parse(data,"A",a)) errormsg="could not parse A";
    present=Tools::findKeyword(data,"B");
    if(present && !Tools::parse(data,"B",b)) errormsg="could not parse B";
    c=pow(2., static_cast<double>(a)/static_cast<double>(b) ) - 1; 
    d = -static_cast<double>(b) / static_cast<double>(a);
  } 
  else if(name=="Q") {
    type=nativeq; 
    beta = 50.0;  // nm-1
    lambda = 1.8; // unitless
    present=Tools::findKeyword(data,"BETA");
    if(present && !Tools::parse(data, "BETA", beta)) errormsg="could not parse BETA";
    present=Tools::findKeyword(data,"LAMBDA");
    if(present && !Tools::parse(data, "LAMBDA", lambda)) errormsg="could not parse LAMBDA";
    bool found_ref=Tools::parse(data,"REF",ref); // nm
    if(!found_ref) errormsg="REF (reference disatance) is required for native Q";

  }
  else if(name=="EXP") type=exponential;
  else if(name=="GAUSSIAN") type=gaussian;
  else if(name=="CUBIC") type=cubic;
  else if(name=="TANH") type=tanh;
  else if(name=="STEP") type=step;
  else if(name=="MODEXP") {
     type=modulated_exponential;
     do_division = false;
     r1 = 0.0; r2 = 0.0;
     if(!Tools::parse(data, "RLOWER", r1)) 
        errormsg="RLOWER is required for MODEXP";
     if(!Tools::parse(data, "RUPPER", r2)) 
        errormsg="RUPPER is required for MODEXP";
     if(fabs(r1 - r2) < 1e-6) errormsg = "RLOWER and RUPPER are too close";
     if(r1 <= 0.0 || r2 <= 0.0) 
        errormsg = "RLOWER and RUPPER should be postive";
     if(r1 > r2) errormsg = "RUPPER should be larger than RLOWER";
     //invr0 is the kappa
     double common = exp(-invr0 * (r1-d0)) / (r2 - r1) / (r2 - r1) / (r2 - r1);
     a3 = common * (2 + invr0 * r1 - invr0 * r2);
     a2 = -common * (
           3 * r1 + invr0 * r1 * r1 + 3 * r2 + 
           invr0 * r1 * r2 - 2 * invr0 * r2 * r2);
     a1 = common * r2 * (6 * r1 + 2 * invr0 * r1 * r1 
           - invr0 * r1 * r2 - invr0 * r2 * r2);
     a0 = -common * r2 * r2 * (3 * r1 + invr0 * r1 * r1 - r2 - invr0 * r1 * r2);
     dmax = r2;
     dmax_2 = r2 * r2;
  }
  else if(name=="POLY3") {
     type = polynomial3;
     do_division = false;
     dmax = r0 + d0; dmax_2 = dmax * dmax;
  }
  else if(name=="POLY3PLT") {
     type = polynomial3_plateau;
     do_division = false;
     dmax = r0 + d0; dmax_2 = dmax * dmax;
  }
  else if(name == "CONSTANT") {
     type = constant;
     do_division = false;
     invr0 = r0 = invr0_2 = 1.0;
  }
#ifdef __PLUMED_HAS_MATHEVAL
  else if(name=="MATHEVAL"){
    type=matheval;
    std::string func;
    Tools::parse(data,"FUNC",func);
    evaluator=evaluator_create(const_cast<char*>(func.c_str()));
    char **check_names;
    int    check_count;
    evaluator_get_variables(evaluator,&check_names,&check_count);
    if(check_count!=1){
      errormsg="wrong number of arguments in MATHEVAL switching function";
      return;
    } 
    if(std::string(check_names[0])!="x"){
      errormsg ="argument should be named 'x'";
      return;
    }
    evaluator_deriv=evaluator_derivative(evaluator,const_cast<char*>("x"));
  }
#endif
  else errormsg="cannot understand switching function type '"+name+"'";
  if( !data.empty() ){
      errormsg="found the following rogue keywords in switching function input : ";
      for(unsigned i=0;i<data.size();++i) errormsg = errormsg + data[i] + " "; 
  }

  if(dostretch && dmax!=std::numeric_limits<double>::max()){
    double dummy;
    double s0=calculate(0.0,dummy);
    double sd=calculate(dmax,dummy);
    stretch=1.0/(s0-sd);
    shift=-sd*stretch;
  }
}

std::string SwitchingFunction::description() const {
  std::ostringstream ostr;
  ostr<<1./invr0<<".  Using ";
  if(type==rational){
     ostr<<"rational";
  } else if(type==exponential){
     ostr<<"exponential";
  } else if(type==nativeq){
     ostr<<"nativeq";     
  } else if(type==gaussian){
     ostr<<"gaussian";
  } else if(type==smap){
     ostr<<"smap";
  } else if(type==cubic){
     ostr<<"cubic";
  } else if(type==tanh){
     ostr<<"tanh";
  } else if(type==step){
     ostr<<"step";
  } else if(type==modulated_exponential){
     ostr<<"modulated_exponential";
  } else if(type==polynomial3){
     ostr<<"polynomial3";
  } else if(type==polynomial3_plateau){
     ostr<<"polynomial3_plateau";
  } else if(type==constant) {
     ostr << "constant";
#ifdef __PLUMED_HAS_MATHEVAL
  } else if(type==matheval){
     ostr<<"matheval";
#endif
  } else{
     plumed_merror("Unknown switching function type");
  }
  ostr<<" swiching function with parameters d0="<<d0;
  if(type==rational){
    ostr<<" nn="<<nn<<" mm="<<mm;
  } else if(type==nativeq){
    ostr<<" beta="<<beta<<" lambda="<<lambda<<" ref="<<ref;
  } else if(type==smap){
    ostr<<" a="<<a<<" b="<<b;
  } else if(type==cubic){
    ostr<<" dmax="<<dmax;
#ifdef __PLUMED_HAS_MATHEVAL
  } else if(type==matheval){
     ostr<<" func="<<evaluator_get_string(evaluator);
#endif

  }
  return ostr.str(); 
}

double SwitchingFunction::do_rational(double rdist,double&dfunc,int nn,int mm)const{
      double result;
      if(2*nn==mm){
// if 2*N==M, then (1.0-rdist^N)/(1.0-rdist^M) = 1.0/(1.0+rdist^N)
        double rNdist=Tools::fastpow(rdist,nn-1);
        double iden=1.0/(1+rNdist*rdist);
        dfunc = -nn*rNdist*iden*iden;
        result = iden;
      } else {
        if(rdist>(1.-100.0*epsilon) && rdist<(1+100.0*epsilon)){
           result=nn/mm;
           dfunc=0.5*nn*(nn-mm)/mm;
        }else{
           double rNdist=Tools::fastpow(rdist,nn-1);
           double rMdist=Tools::fastpow(rdist,mm-1);
           double num = 1.-rNdist*rdist;
           double iden = 1./(1.-rMdist*rdist);
           double func = num*iden;
           result = func;
           dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist));
        }
      }
    return result;
}

double SwitchingFunction::calculateSqr(double distance2,double&dfunc)const{
  if(type==rational && nn%2==0 && mm%2==0 && d0==0.0){
    if(distance2>dmax_2){
      dfunc=0.0;
      return 0.0;
    }
    const double rdist_2 = distance2*invr0_2;
    double result=do_rational(rdist_2,dfunc,nn/2,mm/2);
// chain rule:
    dfunc*=2*invr0_2;
// stretch:
    result=result*stretch+shift;
    dfunc*=stretch;
    return result;
  } else {
    double distance=std::sqrt(distance2);
    return calculate(distance,dfunc);
  }
}

double SwitchingFunction::calculate(double distance,double&dfunc)const{
  plumed_massert(init,"you are trying to use an unset SwitchingFunction");
  if(distance>dmax){
    dfunc=0.0;
    return 0.0;
  }
  const double rdist = (distance-d0)*invr0;
  double result;

  if(rdist<=0. && do_trunc){
     result=1.;
     dfunc=0.0;
  }else{
    if(type==smap){
      double sx=c*pow( rdist, a ); 
      result=pow( 1.0 + sx, d ); 
      dfunc=-b*sx/rdist*result/(1.0+sx); 
    } else if(type == constant) {
       result = 1.0;
       dfunc = 0.0;
    } else if(type==rational){
      result=do_rational(rdist,dfunc,nn,mm);
    }else if(type==exponential){
      result=exp(-rdist);
      dfunc=-result;
    }else if(type==nativeq){
      double rdist2 = beta*(distance - lambda * ref);
      double exprdist=exp(rdist2);
      result=1./(1.+exprdist);
      dfunc=-exprdist/(1.+exprdist)/(1.+exprdist);
    }else if(type==gaussian){
      result=exp(-0.5*rdist*rdist);
      dfunc=-rdist*result;
    }else if(type==cubic){
      double tmp1=rdist-1, tmp2=(1+2*rdist);
      result=tmp1*tmp1*tmp2;
      dfunc=2*tmp1*tmp2 + 2*tmp1*tmp1;
    }else if(type==tanh){
      double tmp1=std::tanh(rdist);
      result = 1.0 - tmp1;
      dfunc=-(1-tmp1*tmp1);
    }else if(type==step){
       result = 0.0;
       dfunc = 0.0;
    }else if(type==modulated_exponential){
       if(distance > r2) {
          result = 0.0; dfunc = 0.0;
       } else if(distance > r1) {
          result = distance*(distance*(distance*a3+a2)+a1)+a0;
          dfunc = (distance*(3*distance*a3+2*a2)+a1)/invr0;

//printf("poly: d0 = %f distance = %f f = %f\n",d0, distance, result);  //@@@

       } else {
          result = exp(-rdist);
          dfunc = -result;

//printf("expo: d0 = %f distance = %f rdist = %f f = %f\n",d0,distance, rdist, result);  //@@@

       }
    } else if(type == polynomial3) {
       double t = fabs(rdist) - 0.5;
       double t2 = t * t;
       result = (t*(2 * t2 - 1.5) + 0.5);
       dfunc = static_cast<int>((rdist > 0.0) - (rdist < 0.0)) * (6 * t2 - 1.5);
    } else if(type == polynomial3_plateau) {
       if(distance > d0 or distance < -d0) {
          double rdist_;
          if(distance < -d0)
            rdist_ = (distance + d0) * invr0;
          else rdist_ = rdist; // rdist = (distance + d0) * invr0
          double t = fabs(rdist_) - 0.5;
          double t2 = t * t;
          result = (t*(2 * t2 - 1.5) + 0.5);
          dfunc = static_cast<int>((rdist_ > 0.0) - (rdist_ < 0.0))
                  * (6 * t2 - 1.5);
       } else {
          result = 1.0; dfunc = 0.0;
       }
#ifdef __PLUMED_HAS_MATHEVAL
    }else if(type==matheval){
      result=evaluator_evaluate_x(evaluator,rdist);
      dfunc=evaluator_evaluate_x(evaluator_deriv,rdist);
#endif
    }else plumed_merror("Unknown switching function type");
// this is for the chain rule:
    dfunc*=invr0;
// this is because calculate() sets dfunc to the derivative divided times the distance.
// (I think this is misleading and I would like to modify it - GB)
    if(do_division) dfunc/=distance;
  }

  result=result*stretch+shift;
  dfunc*=stretch;

  return result;
}

SwitchingFunction::SwitchingFunction():
  init(false),
  do_division(true),
  type(rational),
  invr0(0.0),
  d0(0.0),
  dmax(0.0),
  nn(6),
  mm(0),
  a(0.0),
  b(0.0),
  c(0.0),
  d(0.0),
  lambda(0.0),
  beta(0.0),
  ref(0.0),
  invr0_2(0.0),
  dmax_2(0.0),
  stretch(1.0),
  shift(0.0),
  evaluator(NULL),
  evaluator_deriv(NULL)
{
}

SwitchingFunction::SwitchingFunction(const SwitchingFunction&sf):
  init(sf.init),
  do_division(sf.do_division),
  type(sf.type),
  invr0(sf.invr0),
  d0(sf.d0),
  dmax(sf.dmax),
  nn(sf.nn),
  mm(sf.mm),
  a(sf.a),
  b(sf.b),
  c(sf.c),
  d(sf.d),
  lambda(sf.lambda),
  beta(sf.beta),
  ref(sf.ref),
  invr0_2(sf.invr0_2),
  dmax_2(sf.dmax_2),
  stretch(sf.stretch),
  shift(sf.shift),
  evaluator(NULL),
  evaluator_deriv(NULL),
  r1(sf.r1),
  r2(sf.r2),
  a3(sf.a3),
  a2(sf.a2),
  a1(sf.a1),
  a0(sf.a0),
  do_trunc(sf.do_trunc)
{
#ifdef __PLUMED_HAS_MATHEVAL
  if(sf.evaluator) evaluator=evaluator_create(evaluator_get_string(sf.evaluator));
  if(sf.evaluator_deriv) evaluator_deriv=evaluator_create(evaluator_get_string(sf.evaluator_deriv));
#endif
}

void SwitchingFunction::set(int nn,int mm,double r0,double d0){
  init=true;
  do_trunc = true;
  type=rational;
  if(mm==0) mm=2*nn;
  this->nn=nn;
  this->mm=mm;
  this->invr0=1.0/r0;
  this->invr0_2=this->invr0*this->invr0;
  this->d0=d0;
  this->dmax=d0+r0*pow(0.00001,1./(nn-mm));
  this->dmax_2=this->dmax*this->dmax;
  // do stretch by default
  if(true){
    double dummy;
    double s0=calculate(0.0,dummy);
    double sd=calculate(dmax,dummy);
    stretch=1.0/(s0-sd);
    shift=-sd*stretch;
  }
}

void SwitchingFunction::set_d0(double d_0){
  d0 = d_0;
}

double SwitchingFunction::get_r0() const {
  return 1./invr0;
}

double SwitchingFunction::get_d0() const {
  return d0;
}

double SwitchingFunction::get_dmax() const {
  return dmax;
}

double SwitchingFunction::get_dmax2() const {
  return dmax_2;
}

SwitchingFunction::~SwitchingFunction(){
#ifdef __PLUMED_HAS_MATHEVAL
  if(evaluator) evaluator_destroy(evaluator);
  if(evaluator_deriv) evaluator_destroy(evaluator_deriv);
#endif
}


}



