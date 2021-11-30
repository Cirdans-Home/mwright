%% Fractional heat conduction in nonhomogeneous media under perfect thermal contact

p0 = 1;
rho = 0.5;
a1 = 3;
a2 = 1;
k1 = 2;
k2 = 6;
eps = (k1*sqrt(a2))/(k2*sqrt(a1));
alpha = 0.5;

T1 = @(x,t) (p0./(2*sqrt(a1))).*...
    (mwright((x - rho)./sqrt(a1),t,-alpha/2,1-alpha/2) + ...
    ((eps-1)/(eps+1)).*mwright((x + rho)./sqrt(a1),t,-alpha/2,1-alpha/2));
T2 = @(x,t) (eps*p0./((eps+1)*sqrt(a1))).*mwright(abs(x)./sqrt(a2) + rho,t,-alpha/2,1-alpha/2);

xplus = linspace(0,10,100);

t = [0,logspace(-4,0,10)];
maxyval = max(T1(xplus,t(1)));
yval = linspace(0,maxyval,10);

ls = {'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
col = {'r','b','r','b','r','b','r','b','r','b','r','b'};
widths = linspace(2,1,length(t));
figure(1)
clf
for i=1:length(t)
    hold on
    plot(-xplus,T2(xplus,t(i)),sprintf("%s%s",col{i},ls{i}),...
        'LineWidth',widths(i),'DisplayName','removeme');
    plot(xplus,T1(xplus,t(i)),sprintf("%s%s",col{i},ls{i}),...
        'LineWidth',widths(i),'DisplayName',sprintf('t = %1.2e',t(i)));
    axis([-10 10 0 maxyval])
    hold off
end
vline(0,'k-');
xlabel('x')
legend('Interpreter','LaTeX','FontSize',14)