%% HeatMap for :math:`mu \in \mathbb{C}`
% This code produces the heatmaps with the convergence results for the
% computation of the Wright function with complex values of :math:`mu`.

lambda = linspace(-0.6,-0.1,3);
mure = linspace(-3,3,20);
muim = linspace(-3,3,20);

n = 1000;
prec = 128;

%% Set of points for the evaluation
m = 100;
x = linspace(-5,0,m).';
fid = fopen('input.inp','w+');
fprintf(fid,"%d",m);
for i=1:m
    fprintf(fid,"%1.16f\n",x(i));
end

for l = 1:length(lambda)
    for i=1:length(mure)
        for j=1:length(muim)
            stringa = sprintf("../build/zwrighttest %1.16f %1.16f %1.16f %d %d input.inp",...
                lambda(l), mure(i), muim(j), n, prec);
            system(stringa);
            
            system('sed -i "s/- /-/g" wrighttest.out');
            system('sed -i "s/+ /+/g" wrighttest.out');
            
            %% Read the result to MATLAB for plotting
            wseries = zeros(m,1);
            fid = fopen('wrighttest.out','r');
            for k=1:m
                line = fgetl(fid);
                res = sscanf(line,"(%e %e)");
                wseries(k) = res(1) + res(2)*1i;
            end
            fclose(fid);
            
            t = 1;
            [wmatlab,N(i,j,l)] = mwright(x,t,lambda(l),mure(i)+muim(j)*1i);
            
            ERR(i,j,l) = norm(wmatlab-wseries,'inf')./norm(wseries,'inf');
        end
    end
    figure(1)
    subplot(1,3,l)
    hmap = heatmap(round(mure,2),round(muim,2),log10(ERR(:,:,l)));
    xlabel('Re(μ)')
    ylabel('Im(μ)');
    title(sprintf("λ = %1.2f",lambda(l)));
    hmap.CellLabelFormat = '%1.1f';
    figure(2)
    subplot(1,3,l)
    heatmap(mure,muim,N(:,:,l))
    xlabel('Re(μ)')
    ylabel('Im(μ)');
    title(sprintf("λ = %1.2f",lambda(l)));
end


