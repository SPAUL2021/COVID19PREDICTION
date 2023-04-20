clear all
%

parameter = load('p_T_two.txt');

pzero = parameter;

pvalue = load('p_inv_two.txt');
data = load('casescanada.csv');

TT9 = data(:,1);

%
m = 450;
m1 = 365;

TT0(1) = TT9(1);
 t0(1) = 1;
for i = 1:m1-1
    TT0(i+1) = TT9(i+1)-TT9(i);
     t0(i+1) = i+1;
end    


for i = 1:m
    t(i) = i;
    TT(i) = TT9(i); 
end

TTD(1) = TT(1)

for i = 1:m-1
    TTD(i+1) = TT(i+1) - TT(i);
end

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


m11 = 360;
mm = m11/del;
nn2 = mm + ((450-m11)/del) +1;

for i = mm:nn2
    xp(i-mm+1) = i*del;
end

xp1 = xp


ypr = predic(pd,xp,pzero,H);

yprd(1) = TTD(m1);

for i = 1:56
    tp(i) = m1+(i-1);
end

for i = 1:55
    yprd(i+1) = ypr(10*i+1)-ypr(10*(i-1)+1);
end

t1 = length(t)
TTD1 = length(TTD)
xa1 = length(xa)
ybd1 = length(ybd)
tp11 = length(tp)
yprd1 = length(yprd)

box1=[0 0 m11 m11];
box3=[m11 m11 365 365];
box2=[365 365 470 470];
boxy=[0 1200000 1200000 0];

boxy9=[0 14000 14000 0];

figure(41)
subplot(1,3,1)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,ypr,'k.','MarkerSize',8)
hold on
patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,ypr,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{360}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(a)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Total cases','FontSize',16);

subplot(1,3,2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,ypr,'k.','MarkerSize',8)
hold on
patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,ypr,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{360}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(b)','FontSize',16);
xlabel('time (days)','FontSize',16);
%ylabel('Total cases','FontSize',16);

subplot(1,3,3)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,ypr,'k.','MarkerSize',8)
hold on
patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xp,ypr,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{360}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(c)','FontSize',16);
xlabel('time (days)','FontSize',16);
%ylabel('Total cases','FontSize',16);

figure(9)
plot(t,TTD,'ro','MarkerSize',8)
hold on
plot(t,ybd,'b.','MarkerSize',8)
hold on
plot(tp,yprd,'k.','MarkerSize',8)
hold on
patch(box1,boxy9,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy9,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy9,[0 0 1],'FaceAlpha',0.2)
plot(t,TTD,'ro','MarkerSize',8)
hold on
plot(t,ybd,'b.','MarkerSize',8)
hold on
plot(tp,yprd,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{360}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(c)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Daily cases','FontSize',16);

 
 [x, a, ad, b, bd, c, d, r0] = abvalue(pd,xd,H);
 
m22 = 360;
m2 = m22/del;
n2 = m2 + (55/del) +1;

tau = (360-118)/del;

for i = m2:n2
    xbeta(i-m2+1) = i*del;
end

n3 = length(xbeta);

at(1) = a(m2);
bt(1) = b(m2);

for i = 1:n3-1
    at(i+1) = at(i) + del*ad(m2+i-1-tau);
    bt(i+1) = bt(i) + del*bd(m2+i-1-tau);
end   
 
for i = 1:21
    xs(i) = 20*(i-1);
    ys(i) = 1;
end    

 figure(2)
 
 subplot(2,2,1)
 plot(x,a,'b.','MarkerSize',8)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(a)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('\beta^a(t,P^{est}_{\beta})','FontSize',16);

 
 subplot(2,2,2)
 plot(x,b,'b.','MarkerSize',8)
  set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(b)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('\delta^a(t,P^{est}_{\delta})','FontSize',16);

 subplot(2,2,3)
  plot(x,c,'b.','MarkerSize',8)
  hold on
  plot(x,d,'r.','MarkerSize',8)
  legend('E_1^a(t,P^{est}_{E_1})','E_2^a(t,P^{est}_{E_2})')
  set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(c)','FontSize',16);
xlabel('time (days)','FontSize',16);
%ylabel('\delta^a(t,P^{est}_{\delta})','FontSize',16);
 
 subplot(2,2,4)
 plot(x,r0,'r.','MarkerSize',8)
 hold on
 plot(xs,ys,'b.','MarkerSize',8)
  set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(d)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('R_0(t)','FontSize',16);

 
 for i=1:m1
     beta(i) = a(10*(i-1)+1);
     delta(i) = b(10*(i-1)+1);
 end    
 
 beta1 = length(beta);
 T1 = length(TT0);
 t01= length(t0);
 
 N = 37894799;
 
 for i = 1:nn1
     inf(i) = (1/b(i))*(ydb(i)+c(i));
     %ydbb(i) = ydb(i);
     sus(i) = N - inf(i) - yb(i);
 end
 
 figure(3)
 subplot(2,2,1)
 plot(xd,sus,'b.','MarkerSize',8)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(a)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Susceptible','FontSize',16);

 
 subplot(2,2,2)
 plot(xd,inf,'r.','MarkerSize',8)
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(b)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Infected','FontSize',16);

 
 subplot(2,2,3)
 plot3(t0,beta,TT0,'b-','LineWidth',2)
 grid on
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(c)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('\beta^a(t,P^{est}_{\beta})','FontSize',16);
zlabel('Daily cases','FontSize',16)
 
 subplot(2,2,4)
 plot3(t0,delta,TT0,'r-','LineWidth',2)
 grid on
 set(gca,'LineWidth',2,'FontSize',16,'Box','on');
 title('(d)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('\delta^a(t,P^{est}_{\delta})','FontSize',16);
zlabel('Daily cases','FontSize',16)
 
%
figure(42)

%%
subplot(1,3,1)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(0,0,'.','MarkerSize',8,'Color',[0.4940, 0.1840, 0.5560])
hold on
plot(0,0,'k.','MarkerSize',8)
hold on


n2 = (450-365)/del+1;

sumy = zeros(n2,1);
npath = 10;

box1=[0 0 316 316];
box3=[316 316 365 365];
box2=[365 365 470 470];
boxy=[0 1800000 1800000 0];

patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)

for jj = 1:npath
m11 = (365 - npath + 1 ) + 1*(jj-1); % 3xx + npath - 1 = 365
mm = m11/del;
nn2 = mm + ((450-m11)/del);

len = nn2 - mm +1;
xp = zeros(len,1);
for i = mm:nn2
    xp(i-mm+1) = i*del;
end

xp1 = length(xp);


ypr = predic(pd,xp,pzero,H);

nn3 = length(ypr);
m44 = (365-m11)/del+1;
mn = nn3 - m44; 

for k = 1:mn
         k1 = k+m44-1;
    sumy(k) = sumy(k)+ypr(k1);
end   
plot(xp,ypr,'.','MarkerSize',2,'Color',[0.4940, 0.1840, 0.5560])

end

xpp = xp;
xpp1 = length(xpp);
yprav = (1/npath)*sumy;
yprav1 = length(yprav);

plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xpp,yprav,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{paths}','T^{est}_{av10}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(d)','FontSize',16);
xlabel('time (days)','FontSize',16);
ylabel('Total cases','FontSize',16);

%%
subplot(1,3,2)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(0,0,'.','MarkerSize',8,'Color',[0.4940, 0.1840, 0.5560])
hold on
plot(0,0,'k.','MarkerSize',8)
hold on

n2 = (450-365)/del+1;

sumy = zeros(n2,1);
npath = 30;

box1=[0 0 316 316];
box3=[316 316 365 365];
box2=[365 365 470 470];
boxy=[0 1800000 1800000 0];

patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)

for jj = 1:npath
m11 = (365 - npath + 1 ) + 1*(jj-1); % 3xx + npath - 1 = 365
mm = m11/del;
nn2 = mm + ((450-m11)/del);

len = nn2 - mm +1;
xp = zeros(len,1);
for i = mm:nn2
    xp(i-mm+1) = i*del;
end

xp1 = length(xp);


ypr = predic(pd,xp,pzero,H);

nn3 = length(ypr);
m44 = (365-m11)/del+1;
mn = nn3 - m44 ;

for k = 1:mn
         k1 = k+m44-1;
    sumy(k) = sumy(k)+ypr(k1);
end   
plot(xp,ypr,'.','MarkerSize',2,'Color',[0.4940, 0.1840, 0.5560])

end

xpp = xp;
xpp1 = length(xpp);
yprav = (1/npath)*sumy;
yprav1 = length(yprav);

plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xpp,yprav,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{paths}','T^{est}_{av30}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(e)','FontSize',16);
xlabel('time (days)','FontSize',16);
%ylabel('Total cases','FontSize',16);

%%
subplot(1,3,3)
plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(0,0,'.','MarkerSize',8,'Color',[0.4940, 0.1840, 0.5560])
hold on
plot(0,0,'k.','MarkerSize',8)
hold on

n2 = (450-365)/del+1;

sumy = zeros(n2,1);
npath = 50;

box1=[0 0 316 316];
box3=[316 316 365 365];
box2=[365 365 470 470];
boxy=[0 1800000 1800000 0];

patch(box1,boxy,[0.4660 0.8740 0.1880],'FaceAlpha',0.7)
patch(box3,boxy,[1 0 0],'FaceAlpha',0.2)
patch(box2,boxy,[0 0 1],'FaceAlpha',0.2)

for jj = 1:npath
m11 = (365 - npath + 1 ) + 1*(jj-1); % 3xx + npath - 1 = 365
mm = m11/del;
nn2 = mm + ((450-m11)/del);

len = nn2 - mm +1;
xp = zeros(len,1);
for i = mm:nn2
    xp(i-mm+1) = i*del;
end

xp1 = length(xp);


ypr = predic(pd,xp,pzero,H);

nn3 = length(ypr);
m44 = (365-m11)/del+1;
mn = nn3 - m44 ;

for k = 1:mn
         k1 = k+m44-1;
    sumy(k) = sumy(k)+ypr(k1);
end   
plot(xp,ypr,'.','MarkerSize',2,'Color',[0.4940, 0.1840, 0.5560])

end

xpp = xp;
xpp1 = length(xpp);
yprav = (1/npath)*sumy;
yprav1 = length(yprav);

plot(t,TT,'ro','MarkerSize',8)
hold on
plot(xa,yb,'b.','MarkerSize',8)
hold on
plot(xpp,yprav,'k.','MarkerSize',8)
legend('Data', 'T^a','T^{est}_{paths}','T^{est}_{av50}')
set(gca,'LineWidth',2,'FontSize',16,'Box','on');
title('(f)','FontSize',16);
xlabel('time (days)','FontSize',16);
%ylabel('Total cases','FontSize',16);

%%
 
 
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

xb1 = xb(1);

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

ypr(1) = yb(1);
yprd(1) = ydb(1);

for ii = 1:m-1
    ff = -(aa*sc*yprd(ii)*yprd(ii)+bb*sc*ypr(ii)*yprd(ii)+cc*yprd(ii)+dd);
  ypr(ii+1) = ypr(ii)+del*yprd(ii)+0.5*del*del*ff;
 yprd(ii+1) = yprd(ii)+del*ff;
end
    

end 
 
 
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
