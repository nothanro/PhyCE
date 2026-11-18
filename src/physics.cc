// physics.cc - Phy CE (PHYSCE): the physics problems of the forms (focus_menu.cc), their answers and
// their steps. It takes the place of FlowCE's calculus steps (steps.cc) behind the same interface
// (steps.h): answer_before (answer.cc) asks steps_app for a problem's answer, the steps' page
// (focus_steps.cc) asks steps_of for its steps. A problem is a call whose arguments are named
// (torque(r=0.35,F=120,phi=40)) or in the form's order (torque(0.35,120,40)); an argument left
// out is unknown. Lengths in m, masses in kg, times in s, angles turned in rad, angular speeds in
// rad/s, the angle of a force in degrees; any of them can carry a unit instead (1200_rpm, 30_cm, 2_rev).
// The steps write each line as giac reads it (parsed from text, never evaluated): the equation,
// the unknown alone, the numbers in, the value with its unit. Numbers to 4 significant digits.
#if defined PHYSCE
#include "steps.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
using namespace giac;
#define NOINL __attribute__((noinline))

static const double G = 9.8, PI_ = 3.141592653589793;

// ---------------------------------------------------------------- rules: title, formula (0: none)
enum {
  R_CONT, R_GIVEN, R_K1, R_K2, R_K3, R_K4, R_K5, R_REV,
  R_VT, R_AT, R_AC, R_ATOT, R_OMEGA, R_PERIOD, R_FREQ,
  R_ANGLE, R_TORQUE, R_TSIGN, R_TNET,
  R_XCOM, R_YCOM, R_TIP,
  R_BWEIGHT, R_PIVOT, R_BTORQ, R_BFORCE,
  R_INERT, R_PAX, R_TIA, R_PMASS, R_PPUL, R_PLINK, R_PSOLVE, R_PTENS, R_PFALL, R_AMASS, R_ATENS,
  R_KROT, R_WORK, R_TOSTOP, R_ROLL, R_RSOLVE, R_ROMEGA, R_RSLIDE, R_RKE,
  R_LCONS, R_LSTICK, R_KBEFORE, R_KAFTER, R_KLOST, R_KGAIN,
  R_CONV, R_REF, R_PALPHA,
  R_NRULES
};
struct prule { const char * t, * f; };
static const prule PR[] = {
  {"", 0},
  {"What you know", 0},
  {"The equation without the angle", "omega=omega0+alpha*t"},
  {"The equation without the final speed", "Dtheta=omega0*t+1/2*alpha*t^2"},
  {"The equation without the time", "omega^2=omega0^2+2*alpha*Dtheta"},
  {"The equation without the acceleration", "Dtheta=1/2*(omega0+omega)*t"},
  {"The equation without the starting speed", "Dtheta=omega*t-1/2*alpha*t^2"},
  {"In revolutions: divide by 2 pi", "1_rev=2*pi*rad"},
  {"Tangential speed", "v_t=omega*r"},
  {"Tangential acceleration", "a_t=alpha*r"},
  {"Centripetal acceleration", "a_c=omega^2*r"},
  {"Total acceleration", "a=sqrt(a_t^2+a_c^2)"},
  {"Angular speed", "omega=2*pi*f"},
  {"Period", "T=1/f"},
  {"Frequency", "f=1/T"},
  {"The angle between r and F", 0},
  {"Torque", "tau=r*F*sin(phi)"},
  {"Counterclockwise +, clockwise -", 0},
  {"Net torque", "tau_net=tau_ccw-tau_cw"},
  {"Center of mass", "x_com=(sum(m*x))/(sum(m))"},
  {"Center of mass, up and down", "y_com=(sum(m*y))/(sum(m))"},
  {"It tips once its center of mass passes the edge", "tan(theta)=(w/2)/(h/2)"},
  {"The beam's weight acts at its middle", "W=M*g"},
  {"Pivot at the first support: its force drops out", 0},
  {"Torques balance", "tau_ccw=tau_cw"},
  {"Forces balance", "F_up=F_down"},
  {"Moment of inertia", "I=sum(m*r^2)"},
  {"Parallel axis theorem", "I=I_com+M*h^2"},
  {"Newton's second law for rotation", "tau_net=I*alpha"},
  {"The hanging mass: down is +", "m*g-T=m*a"},
  {"The pulley", "T*R=I*alpha"},
  {"The rope does not slip", "a=alpha*R"},
  {"Put them together and solve for a", 0},
  {"The tension", "T=m*(g-a)"},
  {"Its speed after falling h", "v^2=2*a*h"},
  {"Each mass, the heavier one going down", "m2*g-T2=m2*a"},
  {"The tensions", "T1=m1*(g+a)"},
  {"Rotational kinetic energy", "K=1/2*I*omega^2"},
  {"Work done by a torque", "W=tau*Dtheta"},
  {"The torque's work takes all of K", "tau*Dtheta=K"},
  {"Energy: mgh becomes both kinetic energies", "M*g*h=1/2*M*v^2+1/2*I*omega^2"},
  {"Put in I and omega, then solve for v", "M*g*h=1/2*(1+k)*M*v^2"},
  {"Rolling without slipping", "omega=v/R"},
  {"Sliding with no friction would give", "v=sqrt(2*g*h)"},
  {"Kinetic energy at the bottom", "K=M*g*h"},
  {"Angular momentum is conserved", "I1*omega1=I2*omega2"},
  {"They grab: angular momentum is conserved", "I1*omega1+I2*omega2=(I1+I2)*omega"},
  {"Kinetic energy before", "K=1/2*I*omega^2"},
  {"Kinetic energy after", 0},
  {"Kinetic energy lost (to the grab)", 0},
  {"Kinetic energy gained (work done changing shape)", 0},
  {"Multiply by the conversion factor", 0},
  {"", 0},
  {"The pulley's angular acceleration", "alpha=a/R"},
};

// ---------------------------------------------------------------- state and helpers
static std::vector<step_t> * out_; // the steps (0: the answer only)
static const context * ctx_;
static const char * msg_;          // why a problem cannot be answered
static vecteur * res_;             // the answer: its values

static NOINL gen P(const std::string & s){ return gen(s,ctx_); }
// a-b as a+(-b): giac prints a parsed binary minus with its left product in parentheses ((m*g)-T)
static NOINL gen nominus(const gen & g){
  if (g.type!=_SYMB) return g;
  const gen & f=g._SYMBptr->feuille;
  if (g.is_symb_of_sommet(at_binary_minus) && f.type==_VECT && f._VECTptr->size()==2)
    return symbolic(at_plus,makesequence(nominus(f._VECTptr->front()),symbolic(at_neg,nominus(f._VECTptr->back()))));
  if (f.type!=_VECT) return symbolic(g._SYMBptr->sommet,nominus(f));
  vecteur v(*f._VECTptr);
  for (unsigned i=0;i<v.size();++i) v[i]=nominus(v[i]);
  return symbolic(g._SYMBptr->sommet,gen(v,f.subtype));
}
static NOINL void S1(int r,const std::string & line){
  step_t s; s.rule=(unsigned char)r; s.depth=0; s.rel=0; s.a=nominus(P(line)); s.b=0; s.c=0; s.v=0;
  out_->push_back(s);
}
// a step's line; a=b=c as two lines, a=b then a=c (giac prints a=(b=c))
static NOINL void S(int r,const std::string & line){
  if (!out_) return;
  int d=0,e1=-1,e2=-1;
  for (unsigned i=0;i<line.size();++i){
    const char c=line[i];
    if (c=='(' || c=='[') ++d; else if (c==')' || c==']') --d;
    else if (c==',' && !d){ S1(r,line); return; } // (a list of givens: r=0.35, F=120)
    else if (c=='=' && !d && (!i || (line[i-1]!='<' && line[i-1]!='>' && line[i-1]!='!'))){ if (e1<0) e1=i; else if (e2<0) e2=i; }
  }
  if (e2<0){ S1(r,line); return; }
  const std::string lhs=line.substr(0,e1);
  S1(r,line.substr(0,e2));
  S1(R_CONT,lhs+"="+line.substr(e2+1,line.size()-e2-1));
}
// x with 4 significant digits, as giac reads it ("-0.1047", "2.5e-05", "125.7")
static NOINL std::string fs(double x){
  if (!(x==x)) return "undef";
  if (fabs(x)<1e-12) return "0";
  const bool neg=x<0;
  double a=fabs(x);
  int e=(int)floor(log10(a));
  double m=a/pow(10.0,(double)e);
  m=floor(m*1000+0.5)/1000; // 4 digits of mantissa
  if (m>=10){ m/=10; ++e; }
  char b[32];
  std::string r;
  if (e>=6 || e<=-5){ // m*10^e
    long q=(long)floor(m*1000+0.5);
    sprintf(b,"%ld",q);
    r=std::string(b,1);
    std::string fr(b+1);
    { size_t L=fr.size(); while (L && fr[L-1]=='0') --L; fr=fr.substr(0,L); } // (uSTL's erase(pos) and substr(pos) take npos as is)
    if (!fr.empty()) r+="."+fr;
    sprintf(b,"%d",e);
    r+="e"; r+=b;
  }
  else {
    const int dec=3-e>0?3-e:0;
    long q=(long)floor(a*pow(10.0,(double)dec)+0.5);
    sprintf(b,"%ld",q);
    std::string d(b);
    if (dec>0){
      while ((int)d.size()<=dec) d="0"+d;
      std::string ip=d.substr(0,d.size()-dec),fr=d.substr(d.size()-dec,dec);
      { size_t L=fr.size(); while (L && fr[L-1]=='0') --L; fr=fr.substr(0,L); } // (uSTL's erase(pos) and substr(pos) take npos as is)
      r=fr.empty()?ip:ip+"."+fr;
    }
    else r=d;
  }
  return neg?"-"+r:r;
}
static NOINL std::string fp(double x){ return x<0?"("+fs(x)+")":fs(x); } // in a product or a difference
static NOINL std::string fu(double x,const char * u){ return u && *u?fs(x)+"_"+u:fs(x); } // with its unit
static NOINL void ans(const char * name,double x,const char * u){ if (res_) res_->push_back(P(std::string(name)+"="+fu(x,u))); }

// ---- the arguments: named (r=0.35) or in the problem's order
static vecteur * av_;   // the call's arguments
static const char * nm_; // the problem's names, separated by spaces
static NOINL bool named(const gen & g,const char * s,int n){ return g.type==_IDNT && (int)strlen(g._IDNTptr->id_name)==n && !strncmp(g._IDNTptr->id_name,s,n); }
// the value of argument `name` (0: left out)
static NOINL const gen * arg_(const char * name){
  const int n=strlen(name);
  int pos=0,k=0; // the name's position in the order
  for (const char * p=nm_;*p;++k){
    const char * q=strchr(p,' ');
    const int l=q?q-p:strlen(p);
    if (l==n && !strncmp(p,name,n)){ pos=k; break; }
    p+=l; if (*p) ++p;
    pos=-1;
  }
  int i=0;
  for (const_iterateur it=av_->begin();it!=av_->end();++it){
    if (it->is_symb_of_sommet(at_equal) && it->_SYMBptr->feuille.type==_VECT && it->_SYMBptr->feuille._VECTptr->size()==2){
      if (named(it->_SYMBptr->feuille._VECTptr->front(),name,n)) return &it->_SYMBptr->feuille._VECTptr->back();
    }
    else if (i++==pos) return &*it;
  }
  return 0;
}
static NOINL bool unitis(const gen & u,const char * s){ return u.type==_IDNT && !strcmp(u._IDNTptr->id_name,s); }
// v as a number in SI units. kind: 'a' an angle turned (rad), 'w' an angular speed (rad/s), 'd' a
// force's angle (degrees), 0 anything else (giac's units: to mksa)
static NOINL bool num(gen v,double & d,char kind){
  if (v.is_symb_of_sommet(at_unit) && v._SYMBptr->feuille.type==_VECT && v._SYMBptr->feuille._VECTptr->size()==2){
    const gen & u=v._SYMBptr->feuille._VECTptr->back();
    double f=0;
    if (unitis(u,"_rpm")) f=PI_/30;
    else if (unitis(u,"_tr") || unitis(u,"_rev")) f=2*PI_; // (giac's revolution, and rev)
    else if (unitis(u,"_rps")) f=2*PI_;
    else if (unitis(u,"_deg")) f=kind=='d'?1:PI_/180;
    else if (unitis(u,"_rad")) f=kind=='d'?180/PI_:1;
    if (f){
      if (!num(v._SYMBptr->feuille._VECTptr->front(),d,0)) return false;
      d*=f;
      return true;
    }
    v=mksa_reduce(v,ctx_);
  }
  v=evalf(v,1,ctx_);
  if (v.is_symb_of_sommet(at_unit) && v._SYMBptr->feuille.type==_VECT) v=v._SYMBptr->feuille._VECTptr->front();
  if (v.is_symb_of_sommet(at_prod) || v.is_symb_of_sommet(at_pow)) v=mksa_reduce(v,ctx_); // (2*_N)
  if (v.is_symb_of_sommet(at_unit) && v._SYMBptr->feuille.type==_VECT) v=evalf(v._SYMBptr->feuille._VECTptr->front(),1,ctx_);
  if (v.type==_INT_){ d=v.val; return true; }
  if (v.type==_DOUBLE_){ d=v._DOUBLE_val; return true; }
  if (v.type==_FRAC){ gen e=evalf_double(v,1,ctx_); if (e.type==_DOUBLE_){ d=e._DOUBLE_val; return true; } }
  return false;
}
// argument name as a number: 1, 0 left out (or empty), -1 not a number (msg_ set)
static NOINL int get(const char * name,double & d,char kind=0){
  const gen * g=arg_(name);
  if (!g || (g->type==_VECT && g->_VECTptr->empty())) return 0;
  if (num(*g,d,kind)) return 1;
  msg_="A value is not a number";
  return -1;
}
// argument name as a list of numbers (one number: a list of one); its size, 0 left out, -1 bad
static NOINL int getl(const char * name,std::vector<double> & v,char kind=0){
  v.clear();
  const gen * g=arg_(name);
  if (!g) return 0;
  if (g->type!=_VECT){ double d; if (!num(*g,d,kind)){ msg_="A value is not a number"; return -1; } v.push_back(d); return 1; }
  for (const_iterateur it=g->_VECTptr->begin();it!=g->_VECTptr->end();++it){
    double d;
    if (!num(*it,d,kind)){ msg_="A list holds something that is not a number"; return -1; }
    v.push_back(d);
  }
  return v.size();
}
static NOINL bool need(int k,const char * why){ if (k<=0 && !msg_) msg_=why; return k>0; }
static NOINL std::string vlist(const std::vector<double> & v,const char * u=0){
  std::string s="[";
  for (unsigned i=0;i<v.size();++i){ if (i) s+=","; s+=u?fu(v[i],u):fs(v[i]); }
  return s+"]";
}

// ---------------------------------------------------------------- the problems
// 1. constant angular acceleration: three of Dtheta, omega0, omega, alpha, t give the other two
enum { KT, KW0, KW, KA, KTI }; // theta, omega0, omega, alpha, t
static const char * const KN[5]={"Dtheta","omega0","omega","alpha","t"};
static const char * const KU[5]={"rad","(rad/s)","(rad/s)","(rad/s^2)","s"};
// equation e (KT..KTI: the one without that variable) solved for variable x
static NOINL bool kin_solve(int e,int x,double * k,bool * have){
  const double th=k[KT],w0=k[KW0],w=k[KW],a=k[KA],t=k[KTI];
  double r=0;
  std::string f; // the equation solved for x
  switch (e){
  case KT: // omega=omega0+alpha*t
    if (x==KW){ r=w0+a*t; f="omega=omega0+alpha*t"; S(R_K1,f); S(R_CONT,"omega="+fp(w0)+"+"+fp(a)+"*"+fp(t)); }
    else if (x==KW0){ r=w-a*t; S(R_K1,"omega0=omega-alpha*t"); S(R_CONT,"omega0="+fp(w)+"-"+fp(a)+"*"+fp(t)); }
    else if (x==KA){ if (t==0){ msg_="t is 0: no acceleration to find"; return false; } r=(w-w0)/t; S(R_K1,"alpha=(omega-omega0)/t"); S(R_CONT,"alpha=("+fp(w)+"-"+fp(w0)+")/"+fp(t)); }
    else { if (a==0){ msg_="alpha is 0: the speed never changes"; return false; } r=(w-w0)/a; S(R_K1,"t=(omega-omega0)/alpha"); S(R_CONT,"t=("+fp(w)+"-"+fp(w0)+")/"+fp(a)); }
    break;
  case KW: // Dtheta=omega0*t+1/2*alpha*t^2
    if (x==KT){ r=w0*t+a*t*t/2; S(R_K2,"Dtheta=omega0*t+1/2*alpha*t^2"); S(R_CONT,"Dtheta="+fp(w0)+"*"+fp(t)+"+1/2*"+fp(a)+"*"+fp(t)+"^2"); }
    else if (x==KW0){ if (t==0){ msg_="t is 0"; return false; } r=(th-a*t*t/2)/t; S(R_K2,"omega0=(Dtheta-1/2*alpha*t^2)/t"); S(R_CONT,"omega0=("+fp(th)+"-1/2*"+fp(a)+"*"+fp(t)+"^2)/"+fp(t)); }
    else if (x==KA){ if (t==0){ msg_="t is 0"; return false; } r=2*(th-w0*t)/(t*t); S(R_K2,"alpha=(2*(Dtheta-omega0*t))/t^2"); S(R_CONT,"alpha=(2*("+fp(th)+"-"+fp(w0)+"*"+fp(t)+"))/"+fp(t)+"^2"); }
    else { // 1/2 alpha t^2 + omega0 t - Dtheta = 0
      if (a==0){ if (w0==0){ msg_="No motion: omega0 and alpha are 0"; return false; } r=th/w0; S(R_K2,"t=Dtheta/omega0"); S(R_CONT,"t="+fp(th)+"/"+fp(w0)); break; }
      const double D=w0*w0+2*a*th;
      if (D<0){ msg_="It never turns that far"; return false; }
      r=(-w0+sqrt(D))/a;
      if (r<0) r=(-w0-sqrt(D))/a;
      S(R_K2,"1/2*alpha*t^2+omega0*t-Dtheta=0");
      S(R_CONT,"t=(-omega0+sqrt(omega0^2+2*alpha*Dtheta))/alpha");
      S(R_CONT,"t=(-"+fp(w0)+"+sqrt("+fp(w0)+"^2+2*"+fp(a)+"*"+fp(th)+"))/"+fp(a));
    }
    break;
  case KTI: // omega^2=omega0^2+2*alpha*Dtheta
    if (x==KW){ const double D=w0*w0+2*a*th; if (D<0){ msg_="It stops before turning that far"; return false; } r=sqrt(D); if (w0<0 || (w0==0 && a<0)) r=-r; S(R_K3,"omega=sqrt(omega0^2+2*alpha*Dtheta)"); S(R_CONT,"omega=sqrt("+fp(w0)+"^2+2*"+fp(a)+"*"+fp(th)+")"); }
    else if (x==KW0){ const double D=w*w-2*a*th; if (D<0){ msg_="No starting speed fits"; return false; } r=sqrt(D); if (w<0) r=-r; S(R_K3,"omega0=sqrt(omega^2-2*alpha*Dtheta)"); S(R_CONT,"omega0=sqrt("+fp(w)+"^2-2*"+fp(a)+"*"+fp(th)+")"); }
    else if (x==KA){ if (th==0){ msg_="Dtheta is 0"; return false; } r=(w*w-w0*w0)/(2*th); S(R_K3,"alpha=(omega^2-omega0^2)/(2*Dtheta)"); S(R_CONT,"alpha=("+fp(w)+"^2-"+fp(w0)+"^2)/(2*"+fp(th)+")"); }
    else { if (a==0){ msg_="alpha is 0: the speed never changes"; return false; } r=(w*w-w0*w0)/(2*a); S(R_K3,"Dtheta=(omega^2-omega0^2)/(2*alpha)"); S(R_CONT,"Dtheta=("+fp(w)+"^2-"+fp(w0)+"^2)/(2*"+fp(a)+")"); }
    break;
  case KA: // Dtheta=1/2*(omega0+omega)*t
    if (x==KT){ r=(w0+w)*t/2; S(R_K4,"Dtheta=1/2*(omega0+omega)*t"); S(R_CONT,"Dtheta=1/2*("+fp(w0)+"+"+fp(w)+")*"+fp(t)); }
    else if (x==KW0){ if (t==0){ msg_="t is 0"; return false; } r=2*th/t-w; S(R_K4,"omega0=(2*Dtheta)/t-omega"); S(R_CONT,"omega0=(2*"+fp(th)+")/"+fp(t)+"-"+fp(w)); }
    else if (x==KW){ if (t==0){ msg_="t is 0"; return false; } r=2*th/t-w0; S(R_K4,"omega=(2*Dtheta)/t-omega0"); S(R_CONT,"omega=(2*"+fp(th)+")/"+fp(t)+"-"+fp(w0)); }
    else { if (w0+w==0){ msg_="omega0 + omega is 0"; return false; } r=2*th/(w0+w); S(R_K4,"t=(2*Dtheta)/(omega0+omega)"); S(R_CONT,"t=(2*"+fp(th)+")/("+fp(w0)+"+"+fp(w)+")"); }
    break;
  default: // KW0: Dtheta=omega*t-1/2*alpha*t^2
    if (x==KT){ r=w*t-a*t*t/2; S(R_K5,"Dtheta=omega*t-1/2*alpha*t^2"); S(R_CONT,"Dtheta="+fp(w)+"*"+fp(t)+"-1/2*"+fp(a)+"*"+fp(t)+"^2"); }
    else if (x==KW){ if (t==0){ msg_="t is 0"; return false; } r=(th+a*t*t/2)/t; S(R_K5,"omega=(Dtheta+1/2*alpha*t^2)/t"); S(R_CONT,"omega=("+fp(th)+"+1/2*"+fp(a)+"*"+fp(t)+"^2)/"+fp(t)); }
    else if (x==KA){ if (t==0){ msg_="t is 0"; return false; } r=2*(w*t-th)/(t*t); S(R_K5,"alpha=(2*(omega*t-Dtheta))/t^2"); S(R_CONT,"alpha=(2*("+fp(w)+"*"+fp(t)+"-"+fp(th)+"))/"+fp(t)+"^2"); }
    else { // 1/2 alpha t^2 - omega t + Dtheta = 0
      if (a==0){ if (w==0){ msg_="No motion"; return false; } r=th/w; S(R_K5,"t=Dtheta/omega"); S(R_CONT,"t="+fp(th)+"/"+fp(w)); break; }
      const double D=w*w-2*a*th;
      if (D<0){ msg_="It never turns that far"; return false; }
      r=(w-sqrt(D))/a;
      if (r<0) r=(w+sqrt(D))/a;
      S(R_K5,"1/2*alpha*t^2-omega*t+Dtheta=0");
      S(R_CONT,"t=(omega-sqrt(omega^2-2*alpha*Dtheta))/alpha");
      S(R_CONT,"t=("+fp(w)+"-sqrt("+fp(w)+"^2-2*"+fp(a)+"*"+fp(th)+"))/"+fp(a));
    }
  }
  k[x]=r; have[x]=true;
  S(R_CONT,std::string(KN[x])+"="+fu(r,KU[x]));
  if (x==KT && fabs(r)>1e-9) S(R_REV,fs(r)+"/(2*pi)="+fs(r/(2*PI_))+"_rev");
  ans(KN[x],r,KU[x]);
  return true;
}
static NOINL bool p_kin(){
  double k[5]={0,0,0,0,0};
  bool have[5];
  int n=0;
  std::string given;
  for (int i=0;i<5;++i){
    const int r=get(KN[i],k[i],i==KT?'a':i==KW0 || i==KW?'w':0);
    if (r<0) return false;
    have[i]=r>0; n+=have[i];
    if (have[i]){ if (!given.empty()) given+=","; given+=std::string(KN[i])+"="+fu(k[i],KU[i]); }
  }
  if (n<3){ msg_="Fill in three of the five"; return false; }
  if (n==5){ msg_="All five are filled in: leave out what to find"; return false; }
  S(R_GIVEN,given);
  int u[2],nu=0;
  for (int i=0;i<5;++i) if (!have[i]) u[nu++]=i;
  bool hv[5]; for (int i=0;i<5;++i) hv[i]=have[i]; // (each unknown from the knowns alone)
  double kk[5]; for (int i=0;i<5;++i) kk[i]=k[i];
  if (nu==1){ // any equation with it and none of... the formula sheet's first that has it
    static const signed char order[4]={KT,KW,KTI,KA}; // (their missing variable)
    for (int j=0;j<4;++j) if (order[j]!=u[0]) return kin_solve(order[j],u[0],kk,hv);
  }
  if (!kin_solve(u[1],u[0],kk,hv)) return false; // the first unknown: the equation without the second
  for (int i=0;i<5;++i){ kk[i]=k[i]; hv[i]=have[i]; }
  return kin_solve(u[0],u[1],kk,hv);
}
// 2. a point at r on a spinning object
static NOINL bool p_spin(){
  double r,w,a=0;
  if (!need(get("r",r),"Fill in r") || !need(get("omega",w,'w'),"Fill in omega")) return false;
  const int ha=get("alpha",a);
  if (ha<0) return false;
  S(R_GIVEN,"r="+fu(r,"m")+",omega="+fu(w,"(rad/s)")+(ha?",alpha="+fu(a,"(rad/s^2)"):""));
  const double vt=w*r,at=a*r,ac=w*w*r;
  S(R_VT,"v_t="+fp(w)+"*"+fp(r)); S(R_CONT,"v_t="+fu(vt,"(m/s)")); ans("v_t",vt,"(m/s)");
  S(R_AC,"a_c="+fp(w)+"^2*"+fp(r)); S(R_CONT,"a_c="+fu(ac,"(m/s^2)")); ans("a_c",ac,"(m/s^2)");
  if (ha){
    S(R_AT,"a_t="+fp(a)+"*"+fp(r)); S(R_CONT,"a_t="+fu(at,"(m/s^2)")); ans("a_t",at,"(m/s^2)");
    const double tot=sqrt(at*at+ac*ac);
    S(R_ATOT,"a=sqrt("+fp(at)+"^2+"+fp(ac)+"^2)"); S(R_CONT,"a="+fu(tot,"(m/s^2)")); ans("a",tot,"(m/s^2)");
  }
  return true;
}
// 3. period, frequency, angular speed: one gives the others
static NOINL bool p_period(){
  double T,f,w;
  const int hT=get("T",T),hf=get("f",f),hw=get("omega",w,'w');
  if (hT<0 || hf<0 || hw<0) return false;
  if (hT+hf+hw!=1){ msg_="Fill in exactly one of T, f and omega"; return false; }
  if (hT){ if (T==0){ msg_="T is 0"; return false; } S(R_GIVEN,"T="+fu(T,"s")); f=1/T; S(R_FREQ,"f=1/"+fp(T)); S(R_CONT,"f="+fu(f,"Hz")); w=2*PI_*f; S(R_OMEGA,"omega=2*pi*"+fp(f)); S(R_CONT,"omega="+fu(w,"(rad/s)")); }
  else if (hf){ if (f==0){ msg_="f is 0"; return false; } S(R_GIVEN,"f="+fu(f,"Hz")); T=1/f; S(R_PERIOD,"T=1/"+fp(f)); S(R_CONT,"T="+fu(T,"s")); w=2*PI_*f; S(R_OMEGA,"omega=2*pi*"+fp(f)); S(R_CONT,"omega="+fu(w,"(rad/s)")); }
  else { if (w==0){ msg_="omega is 0"; return false; } S(R_GIVEN,"omega="+fu(w,"(rad/s)")); f=w/(2*PI_); S(R_OMEGA,"f=omega/(2*pi)"); S(R_CONT,"f="+fp(w)+"/(2*pi)"); S(R_CONT,"f="+fu(f,"Hz")); T=1/f; S(R_PERIOD,"T=1/"+fp(f)); S(R_CONT,"T="+fu(T,"s")); }
  if (!hT) ans("T",T,"s");
  if (!hf) ans("f",f,"Hz");
  if (!hw) ans("omega",w,"(rad/s)");
  return true;
}
// 4. the torque of a force: phi, the angle between r and F, or perp, the force's angle from the perpendicular to r
static NOINL bool p_torque(){
  double r,F,phi,pp;
  if (!need(get("r",r),"Fill in r") || !need(get("F",F),"Fill in F")) return false;
  const int hp=get("phi",phi,'d'),hq=get("phi_perp",pp,'d');
  if (hp<0 || hq<0) return false;
  if (hp+hq!=1){ msg_="Fill in one angle: to r, or from perpendicular"; return false; }
  S(R_GIVEN,"r="+fu(r,"m")+",F="+fu(F,"N")+(hp?",phi="+fu(phi,"deg"):",phi_perp="+fu(pp,"deg")));
  if (hq){ phi=90-pp; S(R_ANGLE,"phi=90-phi_perp"); S(R_CONT,"phi=90-"+fp(pp)); S(R_CONT,"phi="+fu(phi,"deg")); }
  const double tau=r*F*sin(phi*PI_/180);
  S(R_TORQUE,"tau="+fp(r)+"*"+fp(F)+"*sin("+fu(phi,"deg")+")");
  S(R_CONT,"tau="+fu(tau,"(N*m)"));
  ans("tau",tau,"(N*m)");
  return true;
}
// 5. net torque: forces F (+ counterclockwise, - clockwise) at distances r, at angles phi (90 left out)
static NOINL bool p_nettorque(){
  std::vector<double> F,r,ph;
  const int n=getl("F",F),nr=getl("r",r),np=getl("phi",ph,'d');
  if (n<0 || nr<0 || np<0) return false;
  if (!n || !nr){ msg_="Fill in the forces and their distances"; return false; }
  if (n!=nr || (np && np!=n)){ msg_="The lists need one value per force"; return false; }
  S(R_GIVEN,"F="+vlist(F,"N")+",r="+vlist(r,"m")+(np?",phi="+vlist(ph,"deg"):""));
  double net=0,ccw=0,cw=0;
  std::string sum;
  for (int i=0;i<n;++i){
    const double p=np?ph[i]:90,t=r[i]*F[i]*sin(p*PI_/180);
    std::string li="tau"+std::string(1,'1'+i%9)+"="+fp(r[i])+"*"+fp(fabs(F[i]))+"*sin("+fu(p,"deg")+")";
    S(i?R_CONT:R_TORQUE,li);
    S(R_CONT,"tau"+std::string(1,'1'+i%9)+"="+fu(fabs(t),"(N*m)"));
    if (t>=0) ccw+=t; else cw-=t;
    net+=t;
  }
  S(R_TNET,"tau_net="+fs(ccw)+"-"+fs(cw));
  S(R_CONT,"tau_net="+fu(net,"(N*m)"));
  S(R_TSIGN,net>0?"tau_net>0":net<0?"tau_net<0":"tau_net=0");
  ans("tau_net",net,"(N*m)");
  return true;
}
// 6. center of mass of masses m at x (and y)
static NOINL bool p_com(){
  std::vector<double> m,x,y;
  const int n=getl("m",m),nx=getl("x",x),ny=getl("y",y);
  if (n<0 || nx<0 || ny<0) return false;
  if (!n || !nx){ msg_="Fill in the masses and their positions"; return false; }
  if (n!=nx || (ny && ny!=n)){ msg_="The lists need one value per mass"; return false; }
  S(R_GIVEN,"m="+vlist(m,"kg")+",x="+vlist(x,"m")+(ny?",y="+vlist(y,"m"):""));
  double M=0,sx=0,sy=0;
  std::string tx,ty,tm;
  for (int i=0;i<n;++i){
    M+=m[i]; sx+=m[i]*x[i]; if (ny) sy+=m[i]*y[i];
    if (i){ tx+="+"; ty+="+"; tm+="+"; }
    tx+=fp(m[i])+"*"+fp(x[i]); if (ny) ty+=fp(m[i])+"*"+fp(y[i]); tm+=fp(m[i]);
  }
  if (M==0){ msg_="The masses add to 0"; return false; }
  S(R_XCOM,"x_com=("+tx+")/("+tm+")"); S(R_CONT,"x_com="+fu(sx/M,"m")); ans("x_com",sx/M,"m");
  if (ny){ S(R_YCOM,"y_com=("+ty+")/("+tm+")"); S(R_CONT,"y_com="+fu(sy/M,"m")); ans("y_com",sy/M,"m"); }
  return true;
}
// 7. a beam on two supports s (positions from its left end) with loads m at x: the supports' forces
static NOINL bool p_beam(){
  double M=0,L=0;
  std::vector<double> s,m,x;
  const int hM=get("M",M),hL=get("L",L),ns=getl("s",s),nm=getl("m",m),nx=getl("x",x);
  if (hM<0 || hL<0 || ns<0 || nm<0 || nx<0) return false;
  if (ns!=2){ msg_="Give the two supports' positions: [x1,x2]"; return false; }
  if (nm!=nx){ msg_="One position per load"; return false; }
  if (hM && M!=0 && !hL){ msg_="Fill in the beam's length (its weight acts at L/2)"; return false; }
  if (s[0]==s[1]){ msg_="The supports are at the same place"; return false; }
  S(R_GIVEN,"s="+vlist(s,"m")+(hM?",M="+fu(M,"kg"):"")+(hL?",L="+fu(L,"m"):"")+(nm?",m="+vlist(m,"kg")+",x="+vlist(x,"m"):""));
  // the weights and where they act (the beam's at its middle)
  std::vector<double> W,X;
  if (hM && M!=0){ W.push_back(M*G); X.push_back(L/2); S(R_BWEIGHT,"W="+fp(M)+"*9.8="+fu(M*G,"N")); S(R_CONT,"x=L/2="+fu(L/2,"m")); }
  for (int i=0;i<nm;++i){ W.push_back(m[i]*G); X.push_back(x[i]); }
  if (W.empty()){ msg_="No loads: fill in the beam's mass or a load"; return false; }
  // torques about support 1: F2*(s2-s1) = sum W*(X-s1)
  const double d2=s[1]-s[0];
  double sum=0,tot=0;
  std::string t;
  for (unsigned i=0;i<W.size();++i){ sum+=W[i]*(X[i]-s[0]); tot+=W[i]; if (i) t+="+"; t+=fp(W[i])+"*"+fp(X[i]-s[0]); }
  S(R_PIVOT,"x1="+fu(s[0],"m"));
  S(R_BTORQ,"F2*"+fp(d2)+"="+t);
  const double F2=sum/d2,F1=tot-F2;
  S(R_CONT,"F2="+fu(F2,"N"));
  S(R_BFORCE,"F1+F2="+fs(tot));
  S(R_CONT,"F1="+fs(tot)+"-"+fp(F2));
  S(R_CONT,"F1="+fu(F1,"N"));
  ans("F1",F1,"N"); ans("F2",F2,"N");
  return true;
}
// 8. a uniform box of height h and width w tipped onto an edge: the angle where it falls over
static NOINL bool p_tip(){
  double h,w;
  if (!need(get("h",h),"Fill in the height") || !need(get("w",w),"Fill in the width")) return false;
  if (h<=0){ msg_="The height must be more than 0"; return false; }
  S(R_GIVEN,"h="+fu(h,"m")+",w="+fu(w,"m"));
  const double th=atan(w/h)*180/PI_;
  S(R_TIP,"tan(theta)="+fp(w/2)+"/"+fp(h/2));
  S(R_CONT,"theta=atan("+fs(w/h)+")");
  S(R_CONT,"theta="+fu(th,"deg"));
  ans("theta",th,"deg");
  return true;
}
// 9. moment of inertia of point masses m at distances r from the axis
static NOINL bool p_inertia(){
  std::vector<double> m,r;
  const int n=getl("m",m),nr=getl("r",r);
  if (n<0 || nr<0) return false;
  if (!n || !nr){ msg_="Fill in the masses and their distances"; return false; }
  if (n!=nr){ msg_="One distance per mass"; return false; }
  S(R_GIVEN,"m="+vlist(m,"kg")+",r="+vlist(r,"m"));
  double I=0;
  std::string t;
  for (int i=0;i<n;++i){ I+=m[i]*r[i]*r[i]; if (i) t+="+"; t+=fp(m[i])+"*"+fp(r[i])+"^2"; }
  S(R_INERT,"I="+t);
  S(R_CONT,"I="+fu(I,"(kg*m^2)"));
  ans("I",I,"(kg*m^2)");
  return true;
}
// 10. parallel axis: three of I_com, M, h, I give the fourth
static NOINL bool p_pax(){
  double Ic,M,h,I;
  const int a=get("I_com",Ic),b=get("M",M),c=get("h",h),d=get("I",I);
  if (a<0 || b<0 || c<0 || d<0) return false;
  if (a+b+c+d!=3){ msg_="Fill in three of the four"; return false; }
  std::string g;
  if (a) g+="I_com="+fu(Ic,"(kg*m^2)");
  if (b) g+=std::string(g.empty()?"":",")+"M="+fu(M,"kg");
  if (c) g+=std::string(g.empty()?"":",")+"h="+fu(h,"m");
  if (d) g+=std::string(g.empty()?"":",")+"I="+fu(I,"(kg*m^2)");
  S(R_GIVEN,g);
  if (!d){ I=Ic+M*h*h; S(R_PAX,"I="+fp(Ic)+"+"+fp(M)+"*"+fp(h)+"^2"); S(R_CONT,"I="+fu(I,"(kg*m^2)")); ans("I",I,"(kg*m^2)"); }
  else if (!a){ Ic=I-M*h*h; S(R_PAX,"I_com=I-M*h^2"); S(R_CONT,"I_com="+fp(I)+"-"+fp(M)+"*"+fp(h)+"^2"); S(R_CONT,"I_com="+fu(Ic,"(kg*m^2)")); ans("I_com",Ic,"(kg*m^2)"); }
  else if (!b){ if (h==0){ msg_="h is 0"; return false; } M=(I-Ic)/(h*h); S(R_PAX,"M=(I-I_com)/h^2"); S(R_CONT,"M=("+fp(I)+"-"+fp(Ic)+")/"+fp(h)+"^2"); S(R_CONT,"M="+fu(M,"kg")); ans("M",M,"kg"); }
  else { if (M==0 || I<Ic){ msg_="I must be at least I_com"; return false; } h=sqrt((I-Ic)/M); S(R_PAX,"h=sqrt((I-I_com)/M)"); S(R_CONT,"h=sqrt(("+fp(I)+"-"+fp(Ic)+")/"+fp(M)+")"); S(R_CONT,"h="+fu(h,"m")); ans("h",h,"m"); }
  return true;
}
// 11. tau = I alpha: two give the third
static NOINL bool p_taualpha(){
  double t,I,a;
  const int ht=get("tau",t),hI=get("I",I),ha=get("alpha",a);
  if (ht<0 || hI<0 || ha<0) return false;
  if (ht+hI+ha!=2){ msg_="Fill in two of the three"; return false; }
  std::string g;
  if (ht) g+="tau="+fu(t,"(N*m)");
  if (hI) g+=std::string(g.empty()?"":",")+"I="+fu(I,"(kg*m^2)");
  if (ha) g+=std::string(g.empty()?"":",")+"alpha="+fu(a,"(rad/s^2)");
  S(R_GIVEN,g);
  if (!ht){ t=I*a; S(R_TIA,"tau="+fp(I)+"*"+fp(a)); S(R_CONT,"tau="+fu(t,"(N*m)")); ans("tau",t,"(N*m)"); }
  else if (!hI){ if (a==0){ msg_="alpha is 0"; return false; } I=t/a; S(R_TIA,"I=tau/alpha"); S(R_CONT,"I="+fp(t)+"/"+fp(a)); S(R_CONT,"I="+fu(I,"(kg*m^2)")); ans("I",I,"(kg*m^2)"); }
  else { if (I==0){ msg_="I is 0"; return false; } a=t/I; S(R_TIA,"alpha=tau/I"); S(R_CONT,"alpha="+fp(t)+"/"+fp(I)); S(R_CONT,"alpha="+fu(a,"(rad/s^2)")); ans("alpha",a,"(rad/s^2)"); }
  return true;
}
// the pulley's I: given, or a solid disk's 1/2 M R^2 (k: I/R^2, the pulley's share of the mass)
static NOINL bool pulley_I(double M,double R,int hI,double I,double & k){
  if (hI){
    if (R<=0){ msg_="Fill in R (I is given)"; return false; }
    k=I/(R*R);
    S(R_CONT,"I/R^2="+fp(I)+"/"+fp(R)+"^2="+fu(k,"kg"));
  }
  else { k=M/2; S(R_CONT,"I=1/2*M*R^2"); }
  return true;
}
// 12. a mass m hanging from a rope around a pulley (mass M, radius R, I: a solid disk's if left out); h: how far it falls
static NOINL bool p_pulley(){
  double m,M=0,R=0,I=0,h=0;
  if (!need(get("m",m),"Fill in the hanging mass")) return false;
  const int hM=get("M",M),hR=get("R",R),hI=get("I",I),hh=get("h",h);
  if (hM<0 || hR<0 || hI<0 || hh<0) return false;
  if (!hM && !hI){ msg_="Fill in the pulley's mass (or its I and R)"; return false; }
  S(R_GIVEN,"m="+fu(m,"kg")+(hM?",M="+fu(M,"kg"):"")+(hR?",R="+fu(R,"m"):"")+(hI?",I="+fu(I,"(kg*m^2)"):"")+(hh?",h="+fu(h,"m"):""));
  S(R_PMASS,"m*g-T=m*a");
  S(R_PPUL,"T*R=I*alpha");
  S(R_PLINK,"alpha=a/R");
  double k;
  if (!pulley_I(M,R,hI,I,k)) return false;
  S(R_CONT,hI?"T=(I/R^2)*a":"T=1/2*M*a");
  const double a=m*G/(m+k),T=k*a;
  S(R_PSOLVE,hI?"a=(m*g)/(m+I/R^2)":"a=(m*g)/(m+1/2*M)");
  S(R_CONT,"a=("+fp(m)+"*9.8)/("+fp(m)+"+"+fp(k)+")");
  S(R_CONT,"a="+fu(a,"(m/s^2)"));
  S(R_PTENS,"T="+fp(m)+"*(9.8-"+fp(a)+")");
  S(R_CONT,"T="+fu(T,"N"));
  ans("a",a,"(m/s^2)"); ans("T",T,"N");
  if (hR && R>0){ S(R_PALPHA,"alpha="+fp(a)+"/"+fp(R)+"="+fu(a/R,"(rad/s^2)")); ans("alpha",a/R,"(rad/s^2)"); }
  if (hh){ const double v=sqrt(2*a*h); S(R_PFALL,"v=sqrt(2*"+fp(a)+"*"+fp(h)+")"); S(R_CONT,"v="+fu(v,"(m/s)")); ans("v",v,"(m/s)"); }
  return true;
}
// 13. an Atwood machine: masses m1 and m2 over a pulley of mass M (radius R, I: a solid disk's if left out)
static NOINL bool p_atwood(){
  double m1,m2,M=0,R=0,I=0;
  if (!need(get("m1",m1),"Fill in m1") || !need(get("m2",m2),"Fill in m2")) return false;
  const int hM=get("M",M),hR=get("R",R),hI=get("I",I);
  if (hM<0 || hR<0 || hI<0) return false;
  S(R_GIVEN,"m1="+fu(m1,"kg")+",m2="+fu(m2,"kg")+(hM?",M="+fu(M,"kg"):"")+(hR?",R="+fu(R,"m"):"")+(hI?",I="+fu(I,"(kg*m^2)"):""));
  const bool sw=m1>m2; // the heavier one, m2, goes down
  const double lo=sw?m2:m1,hi=sw?m1:m2;
  S(R_AMASS,sw?"m1*g-T1=m1*a":"m2*g-T2=m2*a");
  S(R_CONT,sw?"T2-m2*g=m2*a":"T1-m1*g=m1*a");
  S(R_PPUL,"(T2-T1)*R=I*alpha");
  S(R_PLINK,"alpha=a/R");
  double k=0;
  if (hM || hI){ if (!pulley_I(M,R,hI,I,k)) return false; }
  const double a=(hi-lo)*G/(lo+hi+k),Tl=lo*(G+a),Th=hi*(G-a);
  S(R_PSOLVE,sw?"a=((m1-m2)*g)/(m1+m2+I/R^2)":"a=((m2-m1)*g)/(m1+m2+I/R^2)");
  S(R_CONT,"a=(("+fp(hi)+"-"+fp(lo)+")*9.8)/("+fp(m1)+"+"+fp(m2)+"+"+fp(k)+")");
  S(R_CONT,"a="+fu(a,"(m/s^2)"));
  S(R_ATENS,sw?"T2="+fp(lo)+"*(9.8+"+fp(a)+")":"T1="+fp(lo)+"*(9.8+"+fp(a)+")");
  S(R_CONT,sw?"T2="+fu(Tl,"N"):"T1="+fu(Tl,"N"));
  S(R_CONT,sw?"T1="+fp(hi)+"*(9.8-"+fp(a)+")":"T2="+fp(hi)+"*(9.8-"+fp(a)+")");
  S(R_CONT,sw?"T1="+fu(Th,"N"):"T2="+fu(Th,"N"));
  ans("a",a,"(m/s^2)");
  ans("T1",sw?Th:Tl,"N"); ans("T2",sw?Tl:Th,"N");
  return true;
}
// 14. rotational kinetic energy (I, omega); with a torque tau: the work over Dtheta, or the angle that takes all of K
static NOINL bool p_rotke(){
  double I,w,t=0,th=0;
  if (!need(get("I",I),"Fill in I") || !need(get("omega",w,'w'),"Fill in omega")) return false;
  const int ht=get("tau",t),hth=get("Dtheta",th,'a');
  if (ht<0 || hth<0) return false;
  S(R_GIVEN,"I="+fu(I,"(kg*m^2)")+",omega="+fu(w,"(rad/s)")+(ht?",tau="+fu(t,"(N*m)"):"")+(hth?",Dtheta="+fu(th,"rad"):""));
  const double K=I*w*w/2;
  S(R_KROT,"K=1/2*"+fp(I)+"*"+fp(w)+"^2"); S(R_CONT,"K="+fu(K,"J")); ans("K",K,"J");
  if (ht && hth){ const double W=t*th; S(R_WORK,"W="+fp(t)+"*"+fp(th)); S(R_CONT,"W="+fu(W,"J")); ans("W",W,"J"); }
  else if (ht){
    if (t==0){ msg_="tau is 0"; return false; }
    const double a=K/fabs(t);
    S(R_TOSTOP,"Dtheta=K/tau="+fs(K)+"/"+fp(fabs(t))); S(R_CONT,"Dtheta="+fu(a,"rad")); S(R_REV,fs(a)+"/(2*pi)="+fs(a/(2*PI_))+"_rev");
    ans("Dtheta",a,"rad");
  }
  else if (hth){ if (th==0){ msg_="Dtheta is 0"; return false; } const double tq=K/th; S(R_TOSTOP,"tau=K/Dtheta="+fs(K)+"/"+fp(th)); S(R_CONT,"tau="+fu(tq,"(N*m)")); ans("tau",tq,"(N*m)"); }
  return true;
}
// 15. rolling without slipping down a height h: I = k M R^2 (sphere 2/5, disk 1/2, hoop 1, hollow sphere 2/3)
static NOINL bool p_roll(){
  double h,k,R=0,M=0;
  if (!need(get("h",h),"Fill in the height") || !need(get("k",k),"Fill in k: sphere 2/5, disk 1/2, hoop 1")) return false;
  const int hR=get("R",R),hM=get("M",M);
  if (hR<0 || hM<0) return false;
  if (h<0 || k<0){ msg_="h and k must be positive"; return false; }
  S(R_GIVEN,"h="+fu(h,"m")+",k="+fs(k)+(hR?",R="+fu(R,"m"):"")+(hM?",M="+fu(M,"kg"):""));
  S(R_ROLL,"M*g*h=1/2*M*v^2+1/2*(k*M*R^2)*(v/R)^2");
  const double v=sqrt(2*G*h/(1+k));
  S(R_RSOLVE,"v=sqrt((2*g*h)/(1+k))");
  S(R_CONT,"v=sqrt((2*9.8*"+fp(h)+")/(1+"+fs(k)+"))");
  S(R_CONT,"v="+fu(v,"(m/s)"));
  ans("v",v,"(m/s)");
  if (hR && R>0){ S(R_ROMEGA,"omega="+fp(v)+"/"+fp(R)); S(R_CONT,"omega="+fu(v/R,"(rad/s)")); ans("omega",v/R,"(rad/s)"); }
  if (hM){
    const double K=M*G*h,Kt=M*v*v/2;
    S(R_RKE,"K="+fp(M)+"*9.8*"+fp(h)+"="+fu(K,"J"));
    S(R_CONT,"K_trans=1/2*M*v^2="+fu(Kt,"J"));
    S(R_CONT,"K_rot=K-K_trans="+fu(K-Kt,"J"));
    ans("K_trans",Kt,"J"); ans("K_rot",K-Kt,"J");
  }
  S(R_RSLIDE,"v=sqrt(2*9.8*"+fp(h)+")="+fu(sqrt(2*G*h),"(m/s)"));
  return true;
}
// 16. a shape change: I1 omega1 = I2 omega2, three give the fourth; the kinetic energies
static NOINL bool p_shape(){
  double v[4];
  static const char * const n[4]={"I1","omega1","I2","omega2"};
  static const char * const u[4]={"(kg*m^2)","(rad/s)","(kg*m^2)","(rad/s)"};
  int h[4],c=0,x=-1;
  std::string g;
  for (int i=0;i<4;++i){
    h[i]=get(n[i],v[i],i&1?'w':0);
    if (h[i]<0) return false;
    if (h[i]){ ++c; g+=std::string(g.empty()?"":",")+n[i]+"="+fu(v[i],u[i]); } else x=i;
  }
  if (c!=3){ msg_="Fill in three of the four"; return false; }
  S(R_GIVEN,g);
  S(R_LCONS,"I1*omega1=I2*omega2");
  const int o=x^1,p=(x+2)&3,q=p^1; // its partner, and the other side's pair
  if (v[o]==0){ msg_="Division by 0"; return false; }
  v[x]=v[p]*v[q]/v[o];
  S(R_CONT,std::string(n[x])+"=("+n[p]+"*"+n[q]+")/"+n[o]);
  S(R_CONT,std::string(n[x])+"=("+fp(v[p])+"*"+fp(v[q])+")/"+fp(v[o]));
  S(R_CONT,std::string(n[x])+"="+fu(v[x],u[x]));
  ans(n[x],v[x],u[x]);
  const double K1=v[0]*v[1]*v[1]/2,K2=v[2]*v[3]*v[3]/2;
  S(R_KBEFORE,"K1=1/2*"+fp(v[0])+"*"+fp(v[1])+"^2="+fu(K1,"J"));
  S(R_KAFTER,"K2=1/2*"+fp(v[2])+"*"+fp(v[3])+"^2="+fu(K2,"J"));
  if (fabs(K2-K1)>1e-9*fabs(K1)) S(K2>K1?R_KGAIN:R_KLOST,"K2-K1="+fu(K2-K1,"J"));
  return true;
}
// 17. two disks that grab and spin together: I1 at omega1, I2 at omega2 (0 left out)
static NOINL bool p_disks(){
  double I1,w1,I2,w2=0;
  if (!need(get("I1",I1),"Fill in I1") || !need(get("omega1",w1,'w'),"Fill in omega1") || !need(get("I2",I2),"Fill in I2")) return false;
  const int h2=get("omega2",w2,'w');
  if (h2<0) return false;
  if (I1+I2==0){ msg_="I1 + I2 is 0"; return false; }
  S(R_GIVEN,"I1="+fu(I1,"(kg*m^2)")+",omega1="+fu(w1,"(rad/s)")+",I2="+fu(I2,"(kg*m^2)")+",omega2="+fu(w2,"(rad/s)"));
  const double w=(I1*w1+I2*w2)/(I1+I2);
  S(R_LSTICK,"omega=(I1*omega1+I2*omega2)/(I1+I2)");
  S(R_CONT,"omega=("+fp(I1)+"*"+fp(w1)+"+"+fp(I2)+"*"+fp(w2)+")/("+fp(I1)+"+"+fp(I2)+")");
  S(R_CONT,"omega="+fu(w,"(rad/s)"));
  ans("omega",w,"(rad/s)");
  const double K1=I1*w1*w1/2+I2*w2*w2/2,K2=(I1+I2)*w*w/2;
  S(R_KBEFORE,"K1="+fu(K1,"J"));
  S(R_KAFTER,"K2=1/2*("+fp(I1)+"+"+fp(I2)+")*"+fp(w)+"^2="+fu(K2,"J"));
  S(R_KLOST,"K1-K2="+fu(K1-K2,"J"));
  ans("K_lost",K1-K2,"J");
  return true;
}
// the unit conversions (more > Unit conversions): convert(x, kind): kinds as UCONV
struct uconv { const char * from, * to, * fu, * tu; double f; const char * why; };
static const uconv UC[] = {
  {"rpm","rad/s","rpm","(rad/s)",PI_/30,"(2*pi)/60"},
  {"rad/s","rpm","(rad/s)","rpm",30/PI_,"60/(2*pi)"},
  {"rev","rad","rev","rad",2*PI_,"2*pi"},
  {"rad","rev","rad","rev",1/(2*PI_),"1/(2*pi)"},
  {"deg","rad","deg","rad",PI_/180,"pi/180"},
  {"rad","deg","rad","deg",180/PI_,"180/pi"},
  {"rev/s","rad/s","(rev/s)","(rad/s)",2*PI_,"2*pi"},
  {"cm","m","cm","m",0.01,"1/100"},
  {"g","kg","g","kg",0.001,"1/1000"},
  {"km/h","m/s","(km/h)","(m/s)",1/3.6,"1000/3600"},
  {"kg","N","kg","N",G,"g"},
};
enum { NUC = sizeof(UC)/sizeof(UC[0]) };
static NOINL bool p_conv(){
  double x,k;
  if (!need(get("x",x),"Fill in the value") || get("k",k)<=0) { if (!msg_) msg_="Pick a conversion"; return false; }
  const int i=(int)k;
  if (i<0 || i>=NUC){ msg_="No such conversion"; return false; }
  const uconv & c=UC[i];
  S(R_GIVEN,fu(x,c.fu));
  S(R_CONV,fp(x)+"*"+c.why+"="+fu(x*c.f,c.tu));
  if (res_) res_->push_back(P(fu(x*c.f,c.tu)));
  return true;
}

// the formula reference (more > Formula reference): a topic's formulas as steps, each its title
// (when to use it) over the formula
struct fref { const char * when, * f; };
static const fref REF[] = { // topics separated by {0,0}
  {"Angle, speed, acceleration","omega=(Dtheta)/t"},{"","alpha=(Domega)/t"},{"Revolutions and degrees","1_rev=2*pi*rad"},{"Period and frequency","omega=2*pi*f"},{"","T=1/f"},{0,0},
  {"Missing the angle","omega=omega0+alpha*t"},{"Missing the final speed","Dtheta=omega0*t+1/2*alpha*t^2"},{"Missing the time","omega^2=omega0^2+2*alpha*Dtheta"},{"Missing the acceleration","Dtheta=1/2*(omega0+omega)*t"},{0,0},
  {"A point at distance r","x=theta*r"},{"","v_t=omega*r"},{"","a_t=alpha*r"},{"Toward the center","a_c=v_t^2/r"},{"","a_c=omega^2*r"},{"Rolling without slipping","v=omega*R"},{0,0},
  {"Torque of one force (phi: between r and F)","tau=r*F*sin(phi)"},{"Several forces","tau_net=tau_ccw-tau_cw"},{"Gravity acts at the center of mass","tau_g=m*g*r*sin(phi)"},{0,0},
  {"Center of mass","x_com=(m1*x1+m2*x2)/(m1+m2)"},{"Static equilibrium","sum(F)=0"},{"Pick the pivot where an unknown force acts","sum(tau)=0"},{0,0},
  {"Point masses","I=sum(m*r^2)"},{"Parallel axis","I=I_com+M*h^2"},{"Hoop","I=M*R^2"},{"Solid disk or cylinder","I=1/2*M*R^2"},{"Solid sphere","I=2/5*M*R^2"},{"Rod, through its center","I=1/12*M*L^2"},{"Rod, through one end","I=1/3*M*L^2"},{0,0},
  {"Newton's second law","tau_net=I*alpha"},{"Rope on a pulley","a=alpha*R"},{"Hanging mass","m*g-T=m*a"},{0,0},
  {"Rotational kinetic energy","K_rot=1/2*I*omega^2"},{"Work by a torque","W=tau*Dtheta"},{"Energy","DeltaE_mech=W_nc"},{"Rolling: both kinds","K=1/2*m*v^2+1/2*I*omega^2"},{0,0},
  {"Angular momentum","L=I*omega"},{"No outside torque","I1*omega1=I2*omega2"},{"A torque over a time","tau*Dt=DeltaL"},{0,0},
};
int phys_ref; // focus_steps.cc: the page shows the formula reference (its parts: FORMULA, not STEP)
static const char * const REFT[] = {"Rotation", "Kinematics", "Linking", "Torque", "Equilibrium", "Inertia", "Dynamics", "Energy", "Momentum"};
static NOINL bool p_ref(){
  double k;
  if (get("k",k)<=0) return false;
  int t=(int)k,i=0;
  for (int n=0;n<t && i<(int)(sizeof(REF)/sizeof(REF[0]));++i) if (!REF[i].when) ++n;
  if (i>=(int)(sizeof(REF)/sizeof(REF[0]))) return false;
  if (out_) for (;REF[i].when;++i){
    step_t s; s.rule=R_REF; s.depth=0; s.rel=0; s.a=REF[i].f?P(REF[i].f):gen(0); s.b=0; s.c=0; s.v=0;
    s.b=gen(i); // (its text: REF[i].when)
    out_->push_back(s);
  }
  return true;
}

// ---------------------------------------------------------------- the problems by name
struct prob { const char * name, * args; bool (*f)(); };
static const prob PROBS[] = {
  {"kinematics","Dtheta omega0 omega alpha t",p_kin},
  {"spinpoint","r omega alpha",p_spin},
  {"period","T f omega",p_period},
  {"torque","r F phi phi_perp",p_torque},
  {"nettorque","F r phi",p_nettorque},
  {"com","m x y",p_com},
  {"beam","s M L m x",p_beam},
  {"tipping","h w",p_tip},
  {"inertia","m r",p_inertia},
  {"parallelaxis","I_com M h I",p_pax},
  {"taualpha","tau I alpha",p_taualpha},
  {"pulley","m M R I h",p_pulley},
  {"atwood","m1 m2 M R I",p_atwood},
  {"rotke","I omega tau Dtheta",p_rotke},
  {"rolling","h k R M",p_roll},
  {"shapechange","I1 omega1 I2 omega2",p_shape},
  {"disks","I1 omega1 I2 omega2",p_disks},
  {"unitconv","x k",p_conv},
  {"reference","k",p_ref},
};
enum { NPROBS = sizeof(PROBS)/sizeof(PROBS[0]) };
// the problem g calls (an index in PROBS), or -1; its arguments in w
static NOINL int prob_of(const gen & g,vecteur & w){
  if (!g.is_symb_of_sommet(at_of) || g._SYMBptr->feuille.type!=_VECT || g._SYMBptr->feuille._VECTptr->size()!=2) return -1;
  const gen & fn=g._SYMBptr->feuille._VECTptr->front(),& as=g._SYMBptr->feuille._VECTptr->back();
  if (fn.type!=_IDNT) return -1;
  for (int k=0;k<NPROBS;++k)
    if (!strcmp(fn._IDNTptr->id_name,PROBS[k].name)){
      w=as.type==_VECT && as.subtype==_SEQ__VECT?*as._VECTptr:vecteur(1,as);
      return k;
    }
  return -1;
}
// run problem k on arguments w: its answer (res), its steps (out); false and msg_ if it cannot
static NOINL bool run(int k,vecteur & w,std::vector<step_t> * out,vecteur * res,GIAC_CONTEXT){
  ctx_=contextptr; out_=out; res_=res; msg_=0; av_=&w; nm_=PROBS[k].args;
  const bool ok=PROBS[k].f();
  out_=0; res_=0;
  return ok;
}
int phys_problem(const char * s){ // focus_steps.cc: the text calls a problem (steps for it)
  for (int k=0;k<NPROBS;++k){
    const int n=strlen(PROBS[k].name);
    if (!strncmp(s,PROBS[k].name,n) && s[n]=='(') return 1;
  }
  return 0;
}
const char * phys_ref_text(const step_t & s){ // the formula reference's line: when to use it
  if (s.rule!=R_REF || s.b.type!=_INT_ || s.b.val<0 || s.b.val>=(int)(sizeof(REF)/sizeof(REF[0]))) return 0;
  return REF[s.b.val].when;
}

// ---------------------------------------------------------------- steps.h
int steps_method,steps_methods;
const char * steps_method_name="";
const char * steps_title(const step_t & s){ if (s.rule==R_REF){ const char * t=phys_ref_text(s); return t?t:""; } return s.rule<R_NRULES?PR[s.rule].t:""; }
const char * steps_formula(const step_t & s){ return s.rule<R_NRULES?PR[s.rule].f:0; }
gen steps_math(const step_t & s,GIAC_CONTEXT){ return s.a; }
bool steps_answer(const std::vector<step_t> & v,gen & r){ return false; }
std::string steps_table(const step_t & s,GIAC_CONTEXT){ return std::string(); }
int steps_release(){ return 0; }
gen steps_derivative(const gen & f,const gen & x,GIAC_CONTEXT){ return derive(f,x,contextptr); }
bool steps_of(const gen & input,std::vector<step_t> & out,GIAC_CONTEXT){
  vecteur w;
  const int k=prob_of(input,w);
  if (k<0) return false;
  vecteur res;
  return run(k,w,&out,&res,contextptr) && !out.empty();
}
// the problem as the steps' page heads it: what is given (the first step's line), or as typed
gen steps_problem(const gen & input,GIAC_CONTEXT){
  vecteur w;
  const int k=prob_of(input,w);
  phys_ref=k>=0 && PROBS[k].f==p_ref;
  if (phys_ref){ // its topic's name
    ctx_=contextptr; av_=&w; nm_=PROBS[k].args; msg_=0;
    double t;
    if (get("k",t)>0 && t>=0 && t<(int)(sizeof(REFT)/sizeof(REFT[0]))) return P(REFT[(int)t]);
  }
  return input;
}
// answer_before: a problem's answer (its values: one alone, else a list of name = value); 0 for
// anything else; *msg why it cannot be answered (the answer is then the message)
gen steps_app(const gen & g,const char ** msg,GIAC_CONTEXT){
  *msg=0;
  vecteur w;
  const int k=prob_of(g,w);
  if (k<0) return 0;
  vecteur res;
  if (!run(k,w,0,&res,contextptr)){ *msg=msg_?msg_:"Fill in the form"; return 0; }
  if (res.empty()) return 0;
  if (res.size()==1){ // its value alone (as area(...) gives its area)
    const gen & e=res.front();
    return e.is_symb_of_sommet(at_equal)?e._SYMBptr->feuille._VECTptr->back():e;
  }
  return gen(res,0);
}
#endif
