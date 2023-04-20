clear all
%

parameter = load('p_T_test.txt');

pvalue = load('p_inv_two.txt');

pzero = parameter;

%
m = 395;

 del = 0.1;
 nn = (m/del) +1;

for i = 1:nn
    u(i) = (i-1)*del;
    xa(i) = u(i);
end

%

H = 30;
vt = 1*rand(H,1);

sumv = 0;
for i = 1:H
    sumv = sumv + vt(i);
end

v = (0.5/sumv)*vt

w = 1*rand(H,1)
b = 1*rand(H,1)
%
v1dt = 1*rand(H,1);

sumv1 = 0;
for i = 1:H
    sumv1 = sumv1 + v1dt(i);
end

v1d = (0.2/sumv1)*v1dt

w1d = 1*rand(H,1)
b1d = 1*rand(H,1)

v2d = 1*rand(H,1)
w2d = 1*rand(H,1)
b2d = 1*rand(H,1)

v3d = 1*rand(H,1)
w3d = 1*rand(H,1)
b3d = 1*rand(H,1)


%
%pinit = [v; w; b; v1d; w1d; b1d; v2d; w2d; b2d; v3d; w3d; b3d]
pinit = pvalue;
pd = pinit
  
%

Niter=1;

for i=1:Niter

    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');

    [pd] = fminunc(@(pd) Ed(pd,xa,pzero,H), pinit,options);
    
    iteration = i;
    iter(i) = i;
    vv1(i) = pd(1);
    vv2(i) = pd(6);
    vv3(i) = pd(11);
    vv4(i) = pd(16);
    vv5(i) = pd(21);
    vv6(i) = pd(26);
    pinit = pd;
    
    vde(i)= Ed(pinit,xa,pzero,H);
      err = Ed(pinit,xa,pzero,H);
   
 end;
 
 pfinal = pd;
  error = err;
  
  fileID = fopen ('p_inv_test.txt','w');
  
     fprintf(fileID,'%12s\n','error H = 30/abcd');  
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

 
 [x, a, b, c, r0] = abvalue(pd,xa,H);
 
 figure(11)
 plot(x,a,'.')
 ylabel('Beta')
 xlabel('Time (days)')
 
 figure(12)
 plot(x,b,'.')
 ylabel('Delta')
 xlabel('Time (days)')
 
 figure(13)
 plot(x,c,'.')
 ylabel('c')
 xlabel('Time (days)')
 
 
 figure(15)
 plot(x,r0,'.')
 ylabel('Rzero')
 xlabel('Time (days)')


%
 
 
function ydfinal = Ed(p,xb,pz,h);
 
 eps = 0.01;

ydu=0;

m=length(xb);

h9 = 20;
       
    for ii = 1:m
        
    q = eps*(xb(ii) + 1);
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
        tt = aa*(1/ef) + (xb(ii)*bb)/ef^2;
        ss = 2*(bb/ef^2) + xb(ii)*fac*((-e2/ef^2) + (2*e2^2/ef^3));
        sums = sums + rr;
        sumd = sumd + tt;
        sumdd = sumdd + ss;
    end
   yb(ii) = xb(ii)*sums;
   ydb(ii) = sumd;
   yddb(ii) = sumdd;
  
  a = 0;
  b = 0;
  bd = 0;
  c2 = 0;
  cd = 0;
  dd = 0;
  for i=1:h
   
  a = a + 2*p(i)./(exp(p(h+i)*q + p(2*h+i)) + exp(-p(h+i)*q - p(2*h+i)));
  b = b + 2*p(3*h+i)./(exp(p(4*h+i)*q+p(5*h+i))+exp(-p(4*h+i)*q-p(5*h+i)));
  
        e1 = exp((p(4*h+i)*q)+p(5*h+i));
        e2 = exp((-p(4*h+i)*q)-p(5*h+i));
        ef = e1 + e2;
        b22 = 2*p(3*h+i)*p(4*h+i)*eps*(e1-e2);
        rr = - b22/ef^2;
        bd = bd + rr;

  c2 = c2 + 2*p(6*h+i)./(exp(p(7*h+i)*q+p(8*h+i))+exp(-p(7*h+i)*q-p(8*h+i)));
      
        e1d = exp((p(7*h+i)*q)+p(8*h+i));
        e2d = exp((-p(7*h+i)*q)-p(8*h+i));
        efd = e1d + e2d;
        b22d = 2*p(6*h+i)*p(7*h+i)*eps*(e1d-e2d);
        rrd = - b22d/efd^2;
        cd = cd + rrd;
        
dd = dd + 2*p(9*h+i)./(exp(p(10*h+i)*q+p(11*h+i))+exp(-p(10*h+i)*q-p(11*h+i)));        

  end
  
  sc = 1/37894799;
  
  % beta = a, delta = b
  aa = a/b;
  bb = a;
  cc = b - a - (bd/b);
  
  ydb(ii) = ydb(ii)+c2;
  
  yddb(ii) = yddb(ii)+cd;

 st = (yddb(ii)+aa*sc*ydb(ii)*ydb(ii)+bb*sc*yb(ii)*ydb(ii)+cc*ydb(ii)+dd)^2;
  
 ydu = ydu + st;
  %ydu = ydu + (a*ydb(ii)+b*yb(ii)+c)^2;
 
 
end
ydfinal = sqrt(ydu/m);
end 
 
 
function [x, a, b, c, r0] = abvalue(p,x,h);
 
 eps = 0.01;

for j=1:length(x)
    
    q = eps*(x(j)+1);
    
  
  aa = 0;
  bb = 0;
  cc = 0;

  for i=1:h
  
      aa = aa + 2*p(i)./ (exp(p(h+i)*q + p(2*h+i)) + exp(-p(h+i)*q - p(2*h+i)) );
	  bb = bb + 2*p(3*h+i)./ (exp(p(4*h+i)*q + p(5*h+i)) + exp(-p(4*h+i)*q - p(5*h+i)) );
      cc = cc + 2*p(6*h+i)./ (exp(p(7*h+i)*q + p(8*h+i)) + exp(-p(7*h+i)*q - p(8*h+i)) );

  end

  a(j) = aa;
  b(j) = bb;
  c(j) = cc;
  r0(j) = aa/bb;
 
 
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
