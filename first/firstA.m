%function [pfinal, error] = first(nval)
%

% Loading data
data = load('casescanada.csv');
parameter = load('p_T_two.txt');
m1 = size(data);
m2 = m1(1)

%m = nval
 m = 365

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

%[pd, iter, vde, err] = opti(ta, TT, H);
 pfinal = parameter
 
 for i = 1:H
     pv(i) = pfinal(i);
     pw(i) = pfinal(H+i);
     pb(i) = pfinal(2*H+i);
 end
 
 %figure(8)
 %loglog(iter,abs(vde))
 
 n1 = H;
 
 for i = 1:n1
     nn(i) = i;
 end
 



 [tb, TTb, TTdb, TTddb] = ytilda(pfinal, H, m);
 
 TTmin = min(TTb)
 
 %if TTmin < 0
 %    eadjust = 10;
 %else
 %    eadjust = 0;
 %end
 
 %error = err + eadjust;
 
 sp = 320;
 
 [tp, TTp] = ypred(pfinal, H, sp);
 
  figure(1)
  subplot(2,2,1)
 plot(nn,pv,'o','MarkerSize',8,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
hold on
 plot(nn,pw,'s','MarkerSize',8,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
hold on
plot(nn,pb,'d','MarkerSize',8,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k')
legend('p_i','w_i','b_i')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(a)','FontSize',16);
xlabel('i','FontSize',16);
%ylabel('Percentile','FontSize',16);
 
  subplot(2,2,2)
 plot(ta,TT,'ro','MarkerSize',8)
 hold on
 plot(tb,TTb,'b.','MarkerSize',8)
 legend('Data','T^a(t,P^{est}_T)')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(b)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Total cases','FontSize',16);
 
  subplot(2,2,3)
 plot(tb,TTdb,'b.','MarkerSize',8)
% legend('Data','T^a(t,P^{est}_T)')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(c)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('dT^a/dt','FontSize',16);
 
  subplot(2,2,4)
 plot(tb,TTddb,'r.','MarkerSize',8)
 %legend('Data','T^a(t,P^{est}_T)')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(d)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('d^2T^a/dt^2','FontSize',16);

%
function [pd, iter, vde, err] = opti(xa, ya, H)

v = fact*rand(H,1);
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


 function [xb, yb, ydb, yddb] = ytilda(p, h, mm)
 
 eps = 0.01;
 
 del = 0.01;
 nn = (mm/del) +1;

for i = 1:nn
    u(i) = (i-1)*del;
    xb(i) = u(i);
    q = eps*(u(i) + 1);
    sums = 0;
    sumd = 0;
    sumdd = 0;
    for j = 1:h
        e1 = 1;
        e2 = exp((-p(h+j)*q)+p(2*h+j));
        ef = e1 + e2;
        aa = 2*p(j);
        bb = 2*p(j)*p(h+j)*eps*e2;
        fac = 2*p(j)*p(h+j)*p(h+j)*eps*eps;
        rr = aa*(1/ef); 
        tt = aa*(1/ef) + (xb(i)*bb)/ef^2;
        ss = 2*(bb/ef^2) + xb(i)*fac*((-e2/ef^2) + (2*e2^2/ef^3));
        sums = sums + rr;
        sumd = sumd + tt;
        sumdd = sumdd + ss;
    end
   yb(i) = xb(i)*sums;
   ydb(i) = sumd;
   yddb(i) = sumdd;
end 
 
 end 
 
 
 
 function [xp, yp] = ypred(p, h, mm)
 
 eps = 0.01;
 
 del = 0.01;
 mm = mm/del;
 nn = mm + (14/del) +1;

for i = mm:nn
    u(i) = (i-1)*del;
    xb(i) = u(i);
    xp(i+1) = i*del;
    q = eps*(u(i) + 1);
    sums = 0;
    sumd = 0;
    sumdd = 0;
    for j = 1:h
        e1 = 1;
        e2 = exp((-p(h+j)*q)+p(2*h+j));
        ef = e1 + e2;
        aa = 2*p(j);
        bb = 2*p(j)*p(h+j)*eps*e2;
        fac = 2*p(j)*p(h+j)*p(h+j)*eps*eps;
        rr = aa*(1/ef); 
        tt = aa*(1/ef) + (xb(i)*bb)/ef^2;
        ss = 2*(bb/ef^2) + xb(i)*fac*((-e2/ef^2) + (2*e2^2/ef^3));
        sums = sums + rr;
        sumd = sumd + tt;
        sumdd = sumdd + ss;
    end
   yb(i) = xb(i)*sums;
   ydb(i) = sumd;
   yddb(i) = sumdd;
   
   yp(i+1) = yb(i) + del*ydb(i) + 0.5*del*del*yddb(i) ;
end 
 
 end 

