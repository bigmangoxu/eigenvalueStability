load force.dat


N=size(force,1);

for n=1:N
    x(n)=force(n,8);
end

t=1:N;
t=(t-1)/1000;

qe=0.5*1.225*34.03*34.03;

c0=force(1,8);
c1=0.00000000001;
lam_r=1.64;
lam_i=12.39;
phi=-1.45;
y=c1*exp(lam_r*t).*sin(lam_i*t+phi);

plotopt=1;
if(plotopt==1) 
    plot(t,log10(abs(x/qe-x(1)/qe)),'r.')
    hold
    plot(t,log10(abs(y/qe)),'k-')
    axis([0 25 -20 5])
    legend('unsteady CFD','ROM (eigen)')
    xlabel('time(sec)')
    ylabel('log(|cl|)')
    axis square
else
    plot(t,x/qe,'r.')
    hold
    plot(t,y/qe,'k-')
    %axis([0 25 0 2.8])
    axis([0 15 -0.0008 0.0008])
    legend('unsteady CFD','ROM (eigen)')
    xlabel('time(sec)')
    ylabel('cl')
        axis square

end


figure;
count=0;
for i=2:24999
if((x(i)-x(i-1))>0 && (x(i)-x(i+1))>0)
count=count+1;
ind(count)=i;
if(count>1)
T(count)=0.001*(i-ind(count-1));
end
end
end
plot(ind*0.001,1./T*2*pi,'r-o');
axis([10 25 12 18])
xlabel('Time(sec)')
ylabel('Angular frequency [rad/s]')
axis square
