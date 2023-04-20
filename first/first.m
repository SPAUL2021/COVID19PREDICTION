function [pfinal, error] = first(nval)
%

% Loading data
data = load('casescanada.csv');
m1 = size(data);
m2 = m1(1)

m = nval

for i = 1:m
    ta(i) = i;
end    

TT9 = data(:,1);
RR9 = data(:,2);
DT9 = data(:,3);

RR(236) = 123592;

STT = 37894799;
fact = TT9(m2)/STT;
for i = 1:m
    %TT(i) = TT9(i)/STT;
     TT(i) = TT9(i);
end    

H = 20;

[pd, iter, vde, err] = opti(ta, TT, H);
 pfinal = pd
 
 figure(8)
 loglog(iter,abs(vde))
 
 n1 = 3*H;
 
 for i = 1:n1
     nn(i) = i;
 end
 
 figure(10)
 plot(nn,pfinal,'ro')

 [tb, TTb] = ytilda(pfinal, H, m);
 
 TTmin = min(TTb)
 
 if TTmin < 0
     eadjust = 6000;
 else
     eadjust = 0;
 end
 
 error = err + eadjust;
 
 figure(3)
 plot(ta,TT,'r+',tb,TTb,'b.')

%
function [pd, iter, vde, err] = opti(xa, ya, H)

%v = fact*rand(H,1);
v = rand(H,1);
w = rand(H,1);
b = rand(H,1);
%

%
pinit = [v; w; b];
pd = pinit;
  
%

Niter=25;

for i=1:Niter

    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');

    [pd] = fminunc(@(pd) Ed(pd,xa,ya,H), pinit,options);
    
    iteration = i;
    iter(i) = i;
  
    pinit = pd;
    
    vde(i)= Ed(pinit,xa,ya,H);
    err = Ed(pinit,xa,ya,H);
   
end
 
end 
%
 
 
function ydfinal = Ed(p,x,y,h)
 
 eps = 0.01;

ydu=0;

m=length(x);
for j=1:m    
    
    q = eps*(x(j)+1);
    
 z=0;
for i=1:h
  
	  z = z + 2*p(i)./ (1 + exp(-p(h+i)*q + p(2*h+i)));

end

  yu(j) = (x(j)*z - y(j))^2;    
  
  
     sums = 0;
    total = 0;
for i = 1:h
        e1 = exp((p(h+i)*q)+p(2*h+i));
        e2 = exp((-p(h+i)*q)-p(2*h+i));
        ef = e1 +e2;
        aa = 2*p(i);
        bb = 2*x(j)*p(h+i)*p(i)*eps*(e1-e2);
        rr = aa*(1.0/ef) - bb/ef^2;
        sums = sums + rr;
        cc = 2*p(h+i)*p(i)*eps*(e1-e2);
        fac = 2*p(i)*p(h+i)^2*eps^2;
        tt = -2*cc/ef^2 - x(j)*fac*((1.0/ef) - 2*(e1-e2)^2/ef^3);
        total = total +tt;
end
   yd(j) = sums;
   ydd(j) = total;  
    
  ydu = ydu + yu(j);
 
 
end
ydfinal = sqrt(ydu/m);
end 


 function [xb, yb] = ytilda(p, h, mm)
 
 eps = 0.01;
 
 del = 0.05;
 nn = (mm/del) +1;

sumf = 0;
for i = 1:nn
    u(i) = (i-1)*del;
    xb(i) = u(i);
    q = eps*(u(i) + 1);
    sums = 0;
    for j = 1:h
        %e1 = exp((p(h+j)*q)+p(2*h+j));
        e1 = 1;
        e2 = exp((-p(h+j)*q)+p(2*h+j));
        ef = e1 + e2;
        aa = 2*p(j);
       % bb = 2*u(i)*p(h+j)*p(j)*eps*(e1-e2);
        rr = aa*(1/ef); %- bb/ef^2;
        sums = sums + rr;
    end
   yb(i) = xb(i)*sums;
end 
 
 end 

end
