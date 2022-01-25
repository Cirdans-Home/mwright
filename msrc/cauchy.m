%% Cauchy Problem
% Solution of the Cauchy problem using the Wright function representation
% of the solution

clear; clc; close all;

% Azimuth and elevation to produce visually pleasing figure
caz =  134.8992;
cel =  21.8656;
% Do you want to plot the solution for every time step? (SLOW!)
dotimesteplot = false;

xlr=5; np=256; dx=2*xlr/np;
x=linspace(-xlr,xlr-dx,np);

subpindex = 1;
for nu = [0.1,0.3,0.5]
        
    Gc = @(x,t) mwright(x,t,-nu,1-nu)/2;
    g = @(x) real((abs(x) <= 1));
    
    if nu == 0.1
        N = 10000;
    else
        N = 5000;
    end
    T = linspace(0,1,N);
    xlr=5; np=256; dx=2*xlr/np;
    x=linspace(-xlr,xlr-dx,np);
    u(1,:) = g(x);
    fu1 = fft(u(1,:));
    
    if (dotimesteplot)
        figure(1)
        plot(x,u(1,:),'r','LineWidth',2);
        axis([min(x) max(x) 0 1])
    end
    for i=2:length(T)
        Gcx = Gc(x,T(i));
        fGcx  = fft(Gcx);
        fgt = fGcx.*fu1;
        u(i,:) = dx*fftshift( ifft( fgt ) );
        if (dotimesteplot)
            figure(1)
            plot(x,u(i,:),'r','LineWidth',2);
            axis([min(x) max(x) 0 5])
            pause(1/25);
        end
    end
    
    [xt1,xt2] = meshgrid(x,T);
    figure(subpindex)
    mesh(xt1,xt2,u)
    xlabel('x');
    ylabel('t');
    title(sprintf('$\\nu = %1.2f$',nu),'Interpreter','latex');
    axis tight
    view(caz,cel);
    set(gcf,'units','centimeters','position',[0,0,5,5]); 
    set(gca,'FontSize',10);
    set(gcf, 'PaperSize', [5 5]);
    print(sprintf("cauchyproblem-%d",subpindex),'-dpdf','-opengl','-r600');
    subpindex = subpindex + 1;
    clear u
end