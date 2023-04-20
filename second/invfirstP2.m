clear all
%

parameter = load('p_T_test.txt');

pzero = parameter;

pvalue = load('p_inv_two.txt');
data = load('casescanada.csv');

TT9 = data(:,1);

%
m = 450;
m1 = 395;

for i = 1:m
    t(i) = i;
    TT(i) = TT9(i); 
end

TTD(1) = TT(1)

for i = 1:m-1
    TTD(i+1) = TT(i+1) - TT(i);
end

TTDM = movmean(TTD,7);

 del = 0.1;
 nn = (m/del) +1;

for i = 1:nn
    u(i) = (i-1)*del;
    xa(i) = u(i);
end

 nn1 = (m1/del) +1;

for i = 1:nn1
    u(i) = (i-1)*del;
    xd(i) = u(i);
end

%


pd = pvalue;
H = 30;
  
%

[yb, ydb, st] = test(pd,xa,pzero,H);

ybd(1) = yb(1);
for i = 1:m-1
    ybd(i+1) = yb(10*i+1)-yb(10*(i-1)+1);
end

m11 = 382;%381
mm = m11/del;
nn2 = mm + ((450-m11)/del) +1;

for i = mm:nn2
    xp(i-mm+1) = i*del;
end

xp1 = xp;


ypr = predic(pd,xp,pzero,H);

yprd(1) = TTD(m1);

for i = 1:56
    tp(i) = m1+(i-1)
end

for i = 1:55
    yprd(i+1) = ypr(10*i+1)-ypr(10*(i-1)+1);
end

t1 = length(t)
TTD1 = length(TTD)
TTDM1= length(TTDM)
xa1 = length(xa)
ybd1 = length(ybd)
tp11 = length(tp)
yprd1 = length(yprd)


 
 [x, a, ad, b, bd, c, d, r0] = abvalue(pd,xd,H);
 
 m33 = m11-154 %225 / 155;
 
 
 
 yprn = predicn(pd,xp,pzero,H,a,ad,b,bd,m11,m33);
 
 yprn1 = length(yprn);
 xpl = length(xp);
 xp1 = xp;
 
 for i = 1:69
    tpn(i) = m11+(i-1);
 end
tpn1 = tpn;

yprnd(1) = TTD(m11);

for i = 1:68
    yprnd(i+1) = yprn(10*i+1)-yprn(10*(i-1)+1);
end

yprnd1 = length(yprnd);


box1=[0 0 382 382];
box3=[382 382 395 395]
box2=[395 395 470 470];
boxy=[0 1200000 1200000 0];
boxy9=[0 14000 14000 0];


figure(52)
subplot(1,3,1)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,yprn,'k.','MarkerSize',8)
hold on
patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,yprn,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{SE}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(b)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Total cases','FontSize',16);

subplot(1,3,2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,yprn,'k.','MarkerSize',8)
hold on
patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,yprn,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{SE}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(c)','FontSize',16);
xlabel('time (days)','FontSize',16);
%ylabel('Total cases','FontSize',16);

subplot(1,3,3)
plot(t,TTDM,'ro','MarkerSize',8)
hold on
plot(t,ybd,'b.','MarkerSize',8)
hold on
plot(tpn,yprnd,'k.','MarkerSize',8)
hold on
patch(box1,boxy9,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy9,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy9,[0 0 1],'FaceAlpha',0.2)
plot(t,TTDM,'ro','MarkerSize',8)
hold on
plot(t,ybd,'b.','MarkerSize',8)
hold on
plot(tpn,yprnd,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{SE}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(d)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Daily cases','FontSize',16);


 
m22 = 360;
m2 = m22/del;
n2 = m2 + (60/del) +1;

tau = (360-118)/del;

for i = m2:n2
    xbeta(i-m2+1) = i*del;
end

n3 = length(xbeta);

at(1) = a(m2)
bt(1) = b(m2)

for i = 1:n3-1
    at(i+1) = at(i) + del*ad(m2+i-1-tau);
    bt(i+1) = bt(i) + del*bd(m2+i-1-tau);
end   

 
 N = 37894799;
 
 for i = 1:nn1
     inf(i) = (1/b(i))*(ydb(i)+c(i));
     ydbb(i) = ydb(i);
     sus(i) = N - inf(i) - yb(i);
 end
 
%
 
 
function [yb, ydb, st] = test(p,xb,pz,h);
 
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

 st(ii) = abs(yddb(ii)+aa*sc*ydb(ii)*ydb(ii)+bb*sc*yb(ii)*ydb(ii)+cc*ydb(ii)+dd);
  
 
end

end 

function ypr = predic(p,xb,pz,h);
 
 eps = 0.01;
 
 del = 0.1;

m=length(xb);

xb1 = xb(1)

h9 = 20;
       
    for ii = 1:1
        
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

 st(ii) = abs(yddb(ii)+aa*sc*ydb(ii)*ydb(ii)+bb*sc*yb(ii)*ydb(ii)+cc*ydb(ii)+dd);
  
    end

ypr(1) = yb(1)
yprd(1) = ydb(1)

for ii = 1:m-1
    ff = -(aa*sc*yprd(ii)*yprd(ii)+bb*sc*ypr(ii)*yprd(ii)+cc*yprd(ii)+dd);
  ypr(ii+1) = ypr(ii)+del*yprd(ii)+0.5*del*del*ff;
 yprd(ii+1) = yprd(ii)+del*ff;
end
    

end 
 
%----------------------------------------------------------
function ypr = predicn(p,xb,pz,h,a,ad,b,bd,m11,m33);
 
 eps = 0.01;
 
 del = 0.1;

m=length(xb);

xb1 = xb(1)

h9 = 20;
       
    for ii = 1:1
        
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
  
  c2 = 0;
  cd = 0;
  dd = 0;
  for i=1:h

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

  
  ydb(ii) = ydb(ii)+c2;
  
  yddb(ii) = yddb(ii)+cd;

 %st(ii) = abs(yddb(ii)+aa*sc*ydb(ii)*ydb(ii)+bb*sc*yb(ii)*ydb(ii)+cc*ydb(ii)+dd);
  
    end
    
tau = (m11-m33)/del;

m2 = m11/del;

at(1) = a(m2);
bt(1) = b(m2);
bdt(1) = bd(m2);

ypr(1) = yb(1)
yprd(1) = ydb(1)

for ii = 1:m-1
  aa = at(ii)/bt(ii);
  bb = at(ii);
  cc = bt(ii) - at(ii) - (bdt(ii)/bt(ii));
    ff = -(aa*sc*yprd(ii)*yprd(ii)+bb*sc*ypr(ii)*yprd(ii)+cc*yprd(ii)+dd);
  ypr(ii+1) = ypr(ii)+del*yprd(ii)+0.5*del*del*ff;
 yprd(ii+1) = yprd(ii)+del*ff;
%     at(ii+1) = at(ii) + del*ad(m2+ii-1-tau);
%    bt(ii+1) = bt(ii) + del*bd(m2+ii-1-tau);
    at(ii+1) = a(m2+ii-tau);
    bt(ii+1) = b(m2+ii-tau);
    bdt(ii+1) = bd(m2+ii-tau);
end
    

end 
%----------------------------------------------------------
 
function [x, a, ad, b, bd, c, d, r0] = abvalue(p,x,h);
 
 eps = 0.01;

for j=1:length(x)
    
    q = eps*(x(j)+1);
    
  
  aa = 0;
  bb = 0;
  cc = 0;
  dd = 0;
  aad = 0;
  bbd = 0;
  
  for i=1:h
  
      aa = aa + 2*p(i)./ (exp(p(h+i)*q + p(2*h+i)) + exp(-p(h+i)*q - p(2*h+i)) );
      
        e1a = exp((p(h+i)*q)+p(2*h+i));
        e2a = exp((-p(h+i)*q)-p(2*h+i));
        efa = e1a + e2a;
        b22a = 2*p(i)*p(h+i)*eps*(e1a-e2a);
        rra = - b22a/efa^2;
        aad = aad + rra;
      
	  bb = bb + 2*p(3*h+i)./ (exp(p(4*h+i)*q + p(5*h+i)) + exp(-p(4*h+i)*q - p(5*h+i)) );
      
        e1 = exp((p(4*h+i)*q)+p(5*h+i));
        e2 = exp((-p(4*h+i)*q)-p(5*h+i));
        ef = e1 + e2;
        b22 = 2*p(3*h+i)*p(4*h+i)*eps*(e1-e2);
        rr = - b22/ef^2;
        bbd = bbd + rr;
      
      cc = cc + 2*p(6*h+i)./ (exp(p(7*h+i)*q + p(8*h+i)) + exp(-p(7*h+i)*q - p(8*h+i)) );
      dd = dd + 2*p(9*h+i)./(exp(p(10*h+i)*q+p(11*h+i))+exp(-p(10*h+i)*q-p(11*h+i))); 
      
  end

  a(j) = aa;
  ad(j) = aad;
  b(j) = bb;
  bd(j) = bbd;
  c(j) = cc;
  d(j) = dd;
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
