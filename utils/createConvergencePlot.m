function createConvergencePlot(type,layer,options,noRuns,fileName)

plotFarField = layer{1}.calc_p_0;
warning('off', 'e3Dss:N_max_reached')
options.Display = 'none';
switch type
    case '2D'
        count = 1;
        Error = zeros(noRuns,2);
        N_arr = 0:noRuns-1;
        [layer,~,~,relTermMaxArr] = e3Dss(layer, options);
        if plotFarField
            p = layer{1}.p_0;
        else
            p = layer{1}.p;
        end
        for N = N_arr
            options.N_max = N;
            layer = e3Dss(layer, options);
            if plotFarField
                p_N = layer{1}.p_0;
            else
                p_N = layer{1}.p;
            end
            Error(count,:) = norm2((p - p_N).')./norm2(p.');
            count = count + 1;
        end
        relTermMaxArr = [relTermMaxArr, zeros(size(relTermMaxArr,1),noRuns-size(relTermMaxArr,2),1)];
        semilogy(N_arr, Error(:,1), N_arr, Error(:,2), N_arr, relTermMaxArr(1,:), N_arr, relTermMaxArr(2,:))
        set(0,'defaulttextinterpreter','latex')
        xlabel('$$N$$')
        ylabel('$$\frac{\|p_1-p_1^{(N)}\|_2}{\|p_1\|_2}$$')
        legend({'$$k_1 = 15\mathrm{m}^{-1}$$', '$$k_1 = 20\mathrm{m}^{-1}$$', 'Relative term max $$k_1 = 15\mathrm{m}^{-1}$$', 'Relative term max $$k_1 = 20\mathrm{m}^{-1}$$'},'interpreter','latex')
        yLim = ylim;
        ylim([yLim(1),1e1])
        if ~isempty(fileName)
            for i = 1:size(Error,2)
                printResultsToFile([fileName '_Errors_' num2str(i)], {'x',N_arr.', 'y', Error(:,i)})
            end
        end
    case '3D'
        %% Create convergence plot
        R_i = layer{1}.R_i;
        nFreqs = numel(options.omega);
        Error = zeros(noRuns,nFreqs);
        N_arr = 0:noRuns-1;
        layer = e3Dss(layer, options);
        if plotFarField
            p = layer{1}.p_0;
        else
            p = layer{1}.p;
        end
        for N = N_arr
            options.N_max = N;
            layer = e3Dss(layer, options);
            if plotFarField
                p_N = layer{1}.p_0;
            else
                p_N = layer{1}.p;
            end
            Error(N+1,:) = norm2((p - p_N).')./norm2(p.');
        end
        N_max = min(find(max(Error,[],2) == 0));
        N_arr = 0:N_max;
        Error = Error(1:N_max+1,:);
        
        k = options.omega/layer{1}.c_f;
        [NN,kk] = meshgrid(N_arr,k);
        B = log10(Error);
        B(B == -Inf) = Inf;
        B(B == Inf) = min(min(B))+min(min(B))/1000;
        surf(kk*R_i,NN,B.','EdgeColor','none')
        
        set(0,'defaulttextinterpreter','latex')
        xlabel('$$k_1R_{0,1}$$')
        ylabel('$$N$$')
        zlabel('$$\|p-p_N\|_2/\|p\|_2$$')
        set(gca, 'Color', 'none')
        set(gca, 'Layer', 'top')
        box on
        
        grid off
        xlim([0 round(max(k*R_i))])
        ylim([min(N_arr) max(N_arr)])
        view(0,90)
        
        h = colorbar('XTickLabel',{'10^{-18}','10^{-16}','10^{-14}','10^{-12}','10^{-10}','10^{-8}','10^{-6}','10^{-4}','10^{-2}','1'}, ...
               'XTick', -18:2:0);
        ylabel(h, '$\frac{\left\|p_1-p_1^{(N)}\right\|_2}{\|p_1\|_2}$','interpreter','latex')
        extraAxisOptions = {...
            'axis on top=true', ...
            'at={(0,0)}', ...
            'colorbar style={ylabel={$\frac{\left\|p_1-p_1^{(N)}\right\|_2}{\|p_1\|_2}$}, ytick={-18,-16,...,0}, yticklabels={$10^{-18}$, $10^{-16}$, $10^{-14}$, $10^{-12}$, $10^{-10}$, $10^{-8}$, $10^{-6}$, $10^{-4}$, $10^{-2}$, $10^0$}}'};
        
        colormap default
        map = colormap;
        map = [1,1,1; map];
        colormap(map)
        if ~isempty(fileName)
            matlab2tikz([fileName '.tex'], 'height', '3.2094in', 'width', '3.2094in', ...
                'extraAxisOptions', extraAxisOptions)
        end
end
warning('on', 'e3Dss:N_max_reached') 
