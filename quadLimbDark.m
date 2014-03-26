% Based on Mandel & Agol's Model for Quadratic Limb Darkening
% Mandel K. & Agol E. 2002, ApJ, 580, L171; please cite this
% paper if you make use of this in your research.  Also, a
% thanks to Gil Nachmani who wrote this routine would be appreciated.
%
% [phi,F] are the observed phase and relative flux
% p is the planet's radius in units of the star's radius (Rs)
% ap is the planet's orbital radius in units of Rs, assuming zero eccentricity
% P is the planet's period in days
% i is the inclination of the orbit in degrees
% gamma1, gamma2 are the Quadratic Limb Darkening coefficients:
% I(r) = 1-gamma1*(1-mu)-gamma2*(1-mu)^2, where: mu=cos(th)=sqrt(1-r^2);
% E.g. gamma1=0.296, gamma2=0.34 for HD209458
% n is the number of phase points in the resulting lightcurve
% percentOfOrbit is the percentage of orbital phase to be used for the flux
% calculations, e.g. for full orbit use 100(%).
%
% Gil Nachmani, April 2011

function [phi,F, Z]=quadLimbDark(p,ap,P,i,gamma1,gamma2,n,percentOfOrbit)

    if nargin<8;
      percentOfOrbit = 150 * ( p + 1 ) / ap/2/pi;
    end

    endVal = P*percentOfOrbit/100;
    t = -endVal : (2*endVal)/n : endVal;
    phi = t/P;
    Z = ap*( sin( 2*pi/P*t ).^2 + ( cos(pi/180*i) * cos(2*pi/P*t) ).^2 ).^(1/2);
    
    c1 = 0;
    c2 = gamma1 + 2*gamma2;
    c3 = 0;
    c4 = -gamma2; % I(r) = 1-sum(cn*(1-mu^(n/2)))
    c0 = 1 - c1 - c2 - c3 - c4;
    ohmega = c0 / ( 0 + 4 ) + c1 / ( 1 + 4 ) + c2 / ( 2 + 4 ) + c3 / ( 3 + 4 ) + c4 / ( 4 + 4 );

    for j = 1 : n
        z = Z(j);
        a = ( z - p )^2;
        b = ( z + p )^2;
        k = sqrt( (1-a) / 4/z/p );
        q = p^2 - z^2;
        k1 = acos( ( 1 - p^2 + z^2 )/ 2/z );
        k0 = acos( ( p^2 + z^2 - 1 )/ 2/p/z );
        
%         disp( strcat('k:', num2str(k) ) )

        % Evaluating lam_e
        if 1+p < z || abs(phi(j)) > (p+1)/ap/2/pi;
            lam_e = 0;
        elseif abs(1-p) < z && z <= 1+p;
            lam_e = 1/pi * ( p^2 * k0 + k1 - 1/2 * sqrt( 4*z^2 - ( 1 + z^2 - p^2 )^2 ) );
        elseif z <= 1-p && z > p-1;
            lam_e = p^2;
        elseif z <= p-1;
            lam_e = 1;
        end
        
%         disp( strcat('lam_e:', num2str(lam_e) ) )


        % Evaluating lam_d and eta_d
        if z >= 1+p || p == 0 || abs(phi(j)) > (p+1)/ap/2/pi; %1
            lam_d = 0;
            eta_d = 0;
        elseif z >= 1/2 + abs(p-1/2) && z < 1+p; %2
            lam_d = lam1(p,z,a,b,k,q);
            eta_d = eta1(p,z,a,b,k1,k0);
        elseif p < 0.5 && z > p && z < 1-p; %3
            lam_d = lam2(p,z,a,b,k,q);
            eta_d = eta2(p,z);
        elseif p < 0.5 && z == 1-p; %4
            lam_d = lam5(p);
            eta_d = eta2(p,z);
        elseif p < 0.5 && z == p; %5
            lam_d = lam4(p);
            eta_d = eta2(p,z);
        elseif p == 0.5 && z == 0.5; %6
            lam_d = 1/3-4/pi/9;
            eta_d = 3/32;
        elseif p > 0.5 && z == p; %7
            lam_d = lam3(p);
            eta_d = eta1(p,z,a,b,k1,k0);
        elseif p > 0.5 && z >= abs(1-p) && z < p; %8
            lam_d = lam1(p,z,a,b,k,q);
            eta_d = eta1(p,z,a,b,k1,k0);
        elseif p < 1 && z > 0 && z <= 1/2 - abs(p-1/2); %9
            lam_d = lam2(p,z,a,b,k,q);
            eta_d = eta2(p,z);
        elseif p < 1 && z == 0; %10
            lam_d = lam6(p);
            eta_d=eta2(p,z);
        elseif p > 1 && z <= p-1; %11
            lam_d = 0;
            eta_d=1/2;
        end
        
%         disp( strcat('lam_d:', num2str(eta_d) ) )


      F(j) = 1 - 1/(4*ohmega) * ( (1-c2) * lam_e + c2 * ( lam_d+2/3 * heaviside(p-z) ) - c4*eta_d );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lam = lam1(p,z,a,b,k,q)
  lam = 1/9/pi/sqrt(p*z) * ( ( (1-b) * ( 2*b + a - 3 ) - 3*q * (b-2) ) * K(k) + 4*p*z*( z^2 + 7*p^2 - 4 ) * E(k) - 3 * q/a * PI( (a-1)/a,k ) );
end

function lam = lam2(p,z,a,b,k,q)
  lam = 2/9/pi/sqrt(1-a) * ( ( 1 - 5*z^2 + p^2 + q^2 ) * K(1/k) + (1 - a)*( z^2 + 7 * p^2 - 4 ) * E(1/k) - 3*q/a * PI( (a-b)/a,1/k) );
end

function lam = lam3(p)
  lam = 1/3 + 16*p/9/pi * ( 2*p^2 - 1 )*E(1/2/p) - ( 1 - 4*p^2 )*( 3 - 8*p^2 )/9/pi/p*K( 1/2/p );
end

function lam = lam4(p)
  lam = 1/3 + 2/9/pi * ( 4*( 2*p^2 - 1 )*E(2*p) + ( 1 - 4*p^2 )*K(2*p) );
end

function lam = lam5(p)
  lam = 2/3/pi * acos(1-2*p) - 4/9/pi*( 3 + 2*p - 8*p^2 ) * sqrt( p*(1 - p) ) - 2/3*heaviside( p - 1/2 );
end

function lam = lam6(p)
  lam = -2/3*(1-p^2)^(3/2);
end

function eta = eta1(p,z,a,b,k1,k0)
  eta = 1/2/pi * ( k1 + 2*eta2(p,z)*k0 - 1/4*( 1 + 5*p^2 + z^2 ) * sqrt( (1-a)*(b-1) ) );  
end

function eta = eta2(p,z)
  eta = p^2/2*(p^2+2*z^2);
end

%
% Functions to evaluate elliptic ingegrals (adapted from Numerical Recipes)
%

% first
function result = K(k)
    phi = pi/2;  % evaluate complete integral
    s = sin(phi); 
    result = s*RF( cos(phi)^2, (1.0-s*k)*(1.0+s*k), 1.0 );
%     disp( strcat( ' K(', num2str(k) ,') = ', num2str(result) ) )
%     disp( strcat('k = ', num2str(k)) )
end

% second
function result = E(k)
    phi = pi/2;
    s = sin(phi); 
    cc = cos(phi)^2;
    q = (1.0-s*k) * (1.0+s*k);
    result = s * ( RF(cc,q,1.0) - (s*k)^2 * RD(cc,q,1.0) / 3.0 );
%     disp( strcat( ' E(', num2str(k) ,') = ', num2str(result) ) )
%     disp( strcat('k = ', num2str(k)) )
end

% third
function result = PI(n, k)    
    n = -n; % sign convension opposite to Abramowitz and Stegun
    phi = pi/2;
    s=sin(phi);
    enss=n*s*s;
    cc=cos(phi)^2;
    q=(1.0-s*k)*(1.0+s*k);
    result = s*(RF(cc,q,1.0)-enss*RJ(cc,q,1.0,1.0+enss)/3.0);
%     disp( strcat( ' PI(', num2str(n) , ',', num2str(k) ,') = ', num2str(result) ) )
%     disp( strcat('k = ', num2str(k)) )
end

% %
% % Use built-in Matlab routines
% %
% % first
% function f = K(k)
%   f = mfun('EllipticK',k);
%   disp( strcat( ' K(k) = ', num2str(f) ) )
% end
% 
% % second
% function f = E(k)
%   f = mfun('EllipticE',k);
%   disp( strcat( ' E(k) = ', num2str(f) ) )
% end
% 
% % third
% function f = PI(n,k)
%   f = mfun('EllipticPi',n,k);
%   disp( strcat( ' PI(n,k) = ', num2str(f) ) )
% end

%
% Computes Carlson?s elliptic integral of the first kind
%
function result = RF(x,y,z) 
    ERRTOL = 0.08;
    TINY   = 1.5e-38;
    BIG    = 3.0e37;
    THIRD  = (1.0/3.0);
    C1     = (1.0/24.0);
    C2     = 0.1;
    C3     = (3.0/44.0);
    C4     = (1.0/14.0);
    
    if min(min(x,y),z) < 0.0 || min(min(x+y,x+z),y+z) < TINY || max(max(x,y),z) > BIG
        disp('invalid arguments in RF()')
        return
    end
    
    xt=x;
    yt=y;
    zt=z;
    
    sqrtx=sqrt(xt);
    sqrty=sqrt(yt);
    sqrtz=sqrt(zt);
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    xt=0.25*(xt+alamb);
    yt=0.25*(yt+alamb);
    zt=0.25*(zt+alamb);
    ave=THIRD*(xt+yt+zt);
    delx=(ave-xt)/ave;
    dely=(ave-yt)/ave;
    delz=(ave-zt)/ave;

    while max(max(abs(delx),abs(dely)),abs(delz)) > ERRTOL
        sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        ave=THIRD*(xt+yt+zt);
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
    end    
    
    e2=delx*dely-delz*delz;
    e3=delx*dely*delz;
        
    result = (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);   
    
%     disp( strcat('RF: ', num2str(result) ) )

end

%
% Computes Carlson?s elliptic integral of the second kind
%
function result = RD(x,y,z) 
    ERRTOL = 0.05;
    TINY   = 1.0e-25;
    BIG    = 4.5e21;
    C1     = 3/14;
    C2     = 1/6;
    C3     = 9/22;
    C4     = 3/26;
    C5     = 0.25*C3;
    C6     = 1.5*C4;
  
    
    if min(x,y) < 0.0 || min(x+y,z) < TINY || max( max(x,y), z ) > BIG
        disp('invalid arguments in RD()')
        return
    end
    
    xt = x;
    yt = y;
    zt = z;
    sum = 0.0;
    fac = 1.0;
    
    sqrtx=sqrt(xt);
    sqrty=sqrt(yt);
    sqrtz=sqrt(zt);
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    sum = sum + fac/(sqrtz*(zt+alamb));
    fac=0.25*fac;
    xt=0.25*(xt+alamb);
    yt=0.25*(yt+alamb);
    zt=0.25*(zt+alamb);
    ave=0.2*(xt+yt+3.0*zt);
    delx=(ave-xt)/ave;
    dely=(ave-yt)/ave;
    delz=(ave-zt)/ave;
    
    while max( max( abs(delx), abs(dely) ), abs(delz) ) > ERRTOL
        sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        sum = sum + fac/(sqrtz*(zt+alamb));
        fac=0.25*fac;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        ave=0.2*(xt+yt+3.0*zt);
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
    end
    
    ea=delx*dely;
    eb=delz*delz;
    ec=ea-eb;
    ed=ea-6.0*eb;
    ee=ed+ec+ec;
    result = 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee) + delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
end

%
% Computes Carlson?s elliptic integral of the third kind
%
function result = RJ(x,y,z,p) 
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
    
    if (min(min(x,y),z) < 0.0 || min(min(x+y,x+z),min(y+z,abs(p))) < TINY || max(max(x,y),max(z,abs(p))) > BIG)
        disp('invalid arguments in RJ()')
        return
    end
        
    sum=0.0;
    fac=1.0;
    if (p > 0.0)
        xt=x;
        yt=y;
        zt=z;
        pt=p;
    else
        xt=min(min(x,y),z);
        zt=max(max(x,y),z);
        yt=x+y+z-xt-zt;
        a=1.0/(yt-p);
        b=a*(zt-yt)*(yt-xt);
        pt=yt+b;
        rho=xt*zt/yt;
        tau=p*pt/yt;
        rcx=RC(rho,tau);
    end
    
    sqrtx=sqrt(xt);
    sqrty=sqrt(yt);
    sqrtz=sqrt(zt);
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)^2;
    beta=pt*(pt+alamb)^2;

    
%     disp( strcat('alpha: ', num2str(alpha) ) )
%     disp( strcat('x: ', num2str(x) ) )


    sum = sum + fac*RC(alpha,beta);
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

    while max(max(abs(delx),abs(dely)), max(abs(delz),abs(delp))) > ERRTOL
        sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)^2;
        beta=pt*(pt+alamb)^2;
        sum = sum + fac*RC(alpha,beta);
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
    end  
        
    ea=delx*(dely+delz)+dely*delz;
    eb=delx*dely*delz;
    ec=delp*delp;
    ed=ea-3.0*ec;
    ee=eb+2.0*delp*(ea-ec);
    result=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
    
    if (p <= 0.0) 
        result=a*(b*result+3.0*(rcx-RF(xt,yt,zt)));
    end
end


%
% Computes Carlson?s degenerate elliptic integral
%
function result = RC(x,y) 
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
    
    if x < 0.0 || y == 0.0 || (x+abs(y)) < TINY || (x+abs(y)) > BIG || (y<-COMP1 && x > 0.0 && x < COMP2)
        disp('invalid arguments in RC()')
        return
    end
    
    if (y > 0.0)
        xt=x;
        yt=y;
        w=1.0;
    else
        xt=x-y;
        yt = -y;
        w=sqrt(x)/sqrt(xt);
    end
    
    alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
    xt=0.25*(xt+alamb);
    yt=0.25*(yt+alamb);
    ave=THIRD*(xt+yt+yt);
    s=(yt-ave)/ave;
    
    while (abs(s) > ERRTOL)
        alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        ave=THIRD*(xt+yt+yt);
        s=(yt-ave)/ave;
    end
    
    result = w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
end
