
/*
 * INCLUDE MATHJS LIBRARY
 */
 math = mathjs();

/*
 * METHODS TO COMPUTE CARLSON'S ELLIPTIC INTEGRALS 
 * ADAPTED FROM NUMERICAL RECIPES
 */

/*
 * Evaluate Carlson's elliptic integral (1st kind)
 */
function RF(x, y, z){

  // define constants
  ERRTOL = 0.08;
  TINY   = 1.5e-38;
  BIG    = 3.0e37;
  THIRD  = (1.0/3.0);
  C1     = (1.0/24.0);
  C2     = 0.1;
  C3     = (3.0/44.0);
  C4     = (1.0/14.0);
  
  // float alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

  if (math.min(math.min(x,y),z) < 0.0 || math.min(math.min(x+y,x+z),y+z) < TINY ||
    math.max(math.max(x,y),z) > BIG)
      console.log("invalid arguments in rf");
  xt=x;
  yt=y;
  zt=z;
  do {
    sqrtx=math.sqrt(xt);
    sqrty=math.sqrt(yt);
    sqrtz=math.sqrt(zt);
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    xt=0.25*(xt+alamb);
    yt=0.25*(yt+alamb);
    zt=0.25*(zt+alamb);
    ave=THIRD*(xt+yt+zt);
    delx=(ave-xt)/ave;
    dely=(ave-yt)/ave;
    delz=(ave-zt)/ave;
  } while (math.max(math.max(math.abs(delx),math.abs(dely)),math.abs(delz)) > ERRTOL);
  e2=delx*dely-delz*delz;
  e3=delx*dely*delz;
  return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/math.sqrt(ave);
}

/*****************************************************************************/

/*
 * Evaluate Carlson's elliptic integral (2nd kind)
 */
function RD(x, y, z){

  // define constants
  ERRTOL = 0.05;
  TINY   = 1.0e-25;
  BIG    = 4.5e21;
  C1     = (3.0/14.0);
  C2     = (1.0/6.0);
  C3     = (9.0/22.0);
  C4     = (3.0/26.0);
  C5     = (0.25*C3);
  C6     = (1.5*C4);

  // float alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
  //   sqrtz,sum,xt,yt,zt;

  if (math.min(x,y) < 0.0 || math.min(x+y,z) < TINY || math.max(math.max(x,y),z) > BIG)
    console.log("invalid arguments in rd");
  xt=x;
  yt=y;
  zt=z;
  sum=0.0;
  fac=1.0;
  do {
    sqrtx=math.sqrt(xt);
    sqrty=math.sqrt(yt);
    sqrtz=math.sqrt(zt);
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    sum += fac/(sqrtz*(zt+alamb));
    fac=0.25*fac;
    xt=0.25*(xt+alamb);
    yt=0.25*(yt+alamb);
    zt=0.25*(zt+alamb);
    ave=0.2*(xt+yt+3.0*zt);
    delx=(ave-xt)/ave;
    dely=(ave-yt)/ave;
    delz=(ave-zt)/ave;
  } while (math.max(math.max(math.abs(delx),math.abs(dely)),math.abs(delz)) > ERRTOL);
  ea=delx*dely;
  eb=delz*delz;
  ec=ea-eb;
  ed=ea-6.0*eb;
  ee=ed+ec+ec;
  return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
    +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*math.sqrt(ave));
}

/*****************************************************************************/

/*
 * Evaluate Carlson's elliptic integral (3rd kind)
 */
function RJ(x, y, z, p){

  // define constants
  ERRTOL = 0.05;
  TINY   = 2.5e-13;
  BIG    = 9.0e11;
  C1     = (3.0/14.0);
  C2     = (1.0/3.0);
  C3     = (3.0/22.0);
  C4     = (3.0/26.0);
  C5     = (0.75*C3);
  C6     = (1.5*C4);
  C7     = (0.5*C2);
  C8     = (C3+C3);

  // float rc(float x, float y);
  // float rf(float x, float y, float z);
  // float a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
  //   ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;

  if (math.min(math.min(x,y),z) < 0.0 || math.min(math.min(x+y,x+z),math.min(y+z,math.abs(p))) < TINY
    || math.max(math.max(x,y),math.max(z,math.abs(p))) > BIG)
      console.log("invalid arguments in rj");
  sum=0.0;
  fac=1.0;
  if (p > 0.0) {
    xt=x;
    yt=y;
    zt=z;
    pt=p;
  } else {
    xt=math.min(math.min(x,y),z);
    zt=math.max(math.max(x,y),z);
    yt=x+y+z-xt-zt;
    a=1.0/(yt-p);
    b=a*(zt-yt)*(yt-xt);
    pt=yt+b;
    rho=xt*zt/yt;
    tau=p*pt/yt;
    rcx=RC(rho,tau);
  }
  do {
    sqrtx=math.sqrt(xt);
    sqrty=math.sqrt(yt);
    sqrtz=math.sqrt(zt);
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    alpha=math.pow( (pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz) ,2);
    beta=pt*math.pow( (pt+alamb), 2);
    sum += fac*RC(alpha,beta);
    fac=0.25*fac;
    xt=0.25*(xt+alamb);
    yt=0.25*(yt+alamb);
    zt=0.25*(zt+alamb);
    pt=0.25*(pt+alamb);
    ave=0.2*(xt+yt+zt+pt+pt);
    delx=(ave-xt)/ave;
    dely=(ave-yt)/ave;
    delz=(ave-zt)/ave;
    delp=(ave-pt)/ave;
  } while (math.max(math.max(math.abs(delx),math.abs(dely)),
    math.max(math.abs(delz),math.abs(delp))) > ERRTOL);
  ea=delx*(dely+delz)+dely*delz;
  eb=delx*dely*delz;
  ec=delp*delp;
  ed=ea-3.0*ec;
  ee=eb+2.0*delp*(ea-ec);
  ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
    +delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*math.sqrt(ave));
  if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-RF(xt,yt,zt)));
  return ans;
}

/*****************************************************************************/


/*
 * Evaluate Carlson's elliptic integral (degenerate)
 */
function RC(x, y){

  // define constants
  ERRTOL = 0.04;
  TINY   = 1.69e-38;
  SQRTNY = 1.3e-19;
  BIG    = 3.e37;
  TNBG   = (TINY*BIG);
  COMP1  = (2.236/SQRTNY);
  COMP2  = (TNBG*TNBG/25.0);
  THIRD  = (1.0/3.0);
  C1     = 0.3;
  C2     = (1.0/7.0);
  C3     = 0.375;
  C4     = (9.0/22.0);

  // float alamb,ave,s,w,xt,yt;

  if (x < 0.0 || y == 0.0 || (x+math.abs(y)) < TINY || (x+math.abs(y)) > BIG ||
    (y<-COMP1 && x > 0.0 && x < COMP2))
      console.log("invalid arguments in rc");
  if (y > 0.0) {
    xt=x;
    yt=y;
    w=1.0;
  } else {
    xt=x-y;
    yt = -y;
    w=math.sqrt(x)/math.sqrt(xt);
  }
  do {
    alamb=2.0*math.sqrt(xt)*math.sqrt(yt)+yt;
    xt=0.25*(xt+alamb);
    yt=0.25*(yt+alamb);
    ave=THIRD*(xt+yt+yt);
    s=(yt-ave)/ave;
  } while (math.abs(s) > ERRTOL);
  return w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/math.sqrt(ave);
}
