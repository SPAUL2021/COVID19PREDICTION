function [pfinal, error] = invfirst(nval)
%
global T TD TDD
%
parameter = load('p_T_two.txt');
data = load('casescanada.csv');

pz = parameter;

%
m = 365;

for i = 1:m
    ta(i) = i;
end    

TT9 = data(:,1);

for i = 1:m
    %TT(i) = TT9(i)/STT;
     T9(i) = TT9(i);
end    


 del = 0.1;
 nn = (m/del) +1;

for i = 1:nn
    u(i) = (i-1)*del;
    xa(i) = u(i);
end

 eps = 0.01;

ydu=0;

mm=length(xa);

h9 = 20;
       
    for ii = 1:mm
        
    q = eps*(xa(ii) + 1);
    sums = 0;
    sumd = 0;
    sumdd = 0;
    for j = 1:h9
        e1 = 1;
        e2 = exp((-pz(h9+j)*q)+pz(2*h9+j));
        ef = e1 + e2;
        aa = 2*pz(j);
        bb = 2*pz(j)*pz(h9+j)*eps*e2;
        fac = 2*pz(j)*pz(h9+j)*pz(h9+j)*eps*eps;
        rr = aa*(1/ef); 
        tt = aa*(1/ef) + (xa(ii)*bb)/ef^2;
        ss = 2*(bb/ef^2) + xa(ii)*fac*((-e2/ef^2) + (2*e2^2/ef^3));
        sums = sums + rr;
        sumd = sumd + tt;
        sumdd = sumdd + ss;
    end
   yb(ii) = xa(ii)*sums;
   ydb(ii) = sumd;
   yddb(ii) = sumdd;

    end
    
    T = yb;
    TD = ydb;
    TDD = yddb;
    
 figure(3)
 plot(ta,T9,'r+',xa,T,'b.')
 
 figure(4)
 plot(xa,TD,'k.')
 
 figure(5)
 plot(xa,TDD,'r.')

%

H = 20;
v = 0.1*rand(H,1);
w = 1*rand(H,1);
b = 1*rand(H,1);
%
v1d = 0.1*rand(H,1);
w1d = 1*rand(H,1);
b1d = 1*rand(H,1);

v2d = 0.1*rand(H,1);
w2d = 1*rand(H,1);
b2d = 1*rand(H,1);

v3d = 0.1*rand(H,1);
w3d = 1*rand(H,1);
b3d = 1*rand(H,1);

%
pinit = [v; w; b; v1d; w1d; b1d; v2d; w2d; b2d; v3d; w3d; b3d];
pd = pinit;
  
%

Niter=10;

for i=1:Niter

    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');

    [pd] = fminunc(@(pd) Ed(pd,xa,H), pinit,options);
    
    iteration = i;
    iter(i) = i;
    vv1(i) = pd(1);
    vv2(i) = pd(6);
    vv3(i) = pd(11);
    vv4(i) = pd(16);
    vv5(i) = pd(21);
    vv6(i) = pd(26);
    pinit = pd;
    
    vde(i)= Ed(pinit,xa,H);
      err = Ed(pinit,xa,H);
   
 end;
 
 pfinal = pd;
  error = err;
  
  fileID = fopen ('p_inv_test.txt','w');
  
     fprintf(fileID,'%12s\n','error H = 20/abc');  
     fprintf(fileID,'%12.8f\n',error); 
     fprintf(fileID,'%12s\n','parameters');
     fprintf(fileID,'%12.8f\n',pfinal);
     fprintf(fileID,'%12s\n','----------');
     
 fclose(fileID);    
 
 figure(9)
 plot(iter, vv1,'.',iter, vv2, 'o', iter, vv3,'+',iter, vv4, '--',iter, vv5,'x', iter, vv6, '*')
 
 figure(8)
 loglog(iter,vde)
 
 for i = 1:12*H
     nn(i) = i;
 end
 

 figure(10)
 plot(nn,pfinal,'o')

 
 [x, a4, b4, c4, d4, r0] = abvalue(pd,xa,H);
 
 figure(11)
 plot(x,a4,'.')
 
 figure(12)
 plot(x,b4,'.')
 
 figure(13)
 plot(x,c4,'.')
 
 figure(14)
 plot(x,d4,'.')
 
 figure(15)
 plot(x,r0,'.')

end
%
 
 
function ydfinal = Ed(p,xb,h);

global T TD TDD
 
 eps = 0.01;

ydu9=0;

m=length(xb);

h9 = 20;
       
    for ii = 1:m
        
    q = eps*(xb(ii) + 1);
 
  a = 0;
  b1 = 0;
  bd = 0;
  cd = 0;
  dd = 0;
  
  for i=1:h
   
      a = a + 2*p(i)./ (exp(p(h+i)*q + p(2*h+i)) + exp(-p(h+i)*q - p(2*h+i)) );
	  b1 = b1 + 2*p(3*h+i)./ (exp(p(4*h+i)*q + p(5*h+i)) + exp(-p(4*h+i)*q - p(5*h+i)) );
        e11 = exp((p(4*h+i)*q)+p(5*h+i));
        e22 = exp((-p(4*h+i)*q)-p(5*h+i));
        eff = e11 + e22;
       
        b22 = 2*p(3*h+i)*p(4*h+i)*eps*(e11-e22);
        rrr = - b22/eff^2;
        bd = bd + rrr;

      cd = cd + 2*p(6*h+i)./ (exp(p(7*h+i)*q + p(8*h+i)) + exp(-p(7*h+i)*q - p(8*h+i)) );
      
      dd = dd + 2*p(9*h+i)./ (exp(p(10*h+i)*q + p(11*h+i)) + exp(-p(10*h+i)*q - p(11*h+i)) );
      

  end
  
  sc = 1/37894799;
  
  % beta = a, delta = b1
  aa1 = a/b1;
  bb1 = a;
  cc = b1 - a - (bd/b1);
  
  %if (TD(ii) > 0)
      
  TD(ii) = TD(ii)+cd;

  ydu9 = ydu9 + (TDD(ii)+aa1*sc*TD(ii)*TD(ii)+bb1*sc*T(ii)*TD(ii)+cc*TD(ii)+dd)^2;
  
  %ydu = ydu + (a*ydb(ii)+b*yb(ii)+c)^2;
  %end
 
end
ydfinal = sqrt(ydu9/m);
end 
 
 
function [x, a3, b3, c3, d3, r0] = abvalue(p,x,h);
 
 eps = 0.01;

for j=1:length(x)
    
    q = eps*(x(j)+1);
    
  
  aa2 = 0;
  bb2 = 0;
  cc2 = 0;
  dd2 = 0;
  for i=1:h
  
      aa2 = aa2 + 2*p(i)./ (exp(p(h+i)*q + p(2*h+i)) + exp(-p(h+i)*q - p(2*h+i)) );
	  bb2 = bb2 + 2*p(3*h+i)./ (exp(p(4*h+i)*q + p(5*h+i)) + exp(-p(4*h+i)*q - p(5*h+i)) );
      cc2 = cc2 + 2*p(6*h+i)./ (exp(p(7*h+i)*q + p(8*h+i)) + exp(-p(7*h+i)*q - p(8*h+i)) );
      dd2 = dd2 + 2*p(9*h+i)./ (exp(p(10*h+i)*q + p(11*h+i)) + exp(-p(10*h+i)*q - p(11*h+i)) );

  end

  a3(j) = aa2;
  b3(j) = bb2;
  c3(j) = cc2;
  d3(j) = dd2;
  r0(j) = aa2/bb2;
 
 
end

 end 


 function [xb, yb] = ytilda(p, h)
 
 eps = 0.01;

sumf = 0;
for i = 1:200
    u(i) = (i-1)*0.05;
    xb(i) = u(i);
    q = eps*(u(i) + 1);
    sums = 0;
    for j = 1:h
        e1 = exp((p(h+j)*q)+p(2*h+j));
        e2 = exp((-p(h+j)*q)-p(2*h+j));
        ef = e1 +e2;
        aa = 2*p(j);
       % bb = 2*u(i)*p(h+j)*p(j)*eps*(e1-e2);
        rr = aa*(1.0/ef); %- bb/ef^2;
        sums = sums + rr;
    end
   yb(i) = xb(i)*sums;
end 
 
 end 

