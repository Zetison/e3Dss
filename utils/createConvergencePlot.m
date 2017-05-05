function createConvergencePlot(type,options,v,noRuns,fileName)
 
switch type
    case '2D'
        count = 1;
        Error = zeros(noRuns,2);
        N_arr = 0:noRuns-1;
        data = e3Dss(v, options);
        p = data(1).p;
        for N = N_arr
            N
            options.N_max = N;
            data = e3Dss(v, options);
            p_N = data(1).p;
            Error(count,:) = norm2((p - p_N).', 1)./norm2(p.', 1);
            count = count + 1;
        end
        semilogy(N_arr, Error(:,1), N_arr, Error(:,2))
        set(0,'defaulttextinterpreter','latex')
        title('Error plots')
        xlabel('$$N$$')
        ylabel('$$\|p-p_N\|_2/\|p\|_2$$')
        legend({'$$k_1 = 15$$', '$$k_1 = 20$$'},'interpreter','latex')
        
        for i = 1:size(Error,2)
            printResultsToFile([fileName '_Errors_' num2str(i)], N_arr.', Error(:,i), NaN, 0, 1)
        end
            
    case '3D'
        %% Create convergence plot
        R_o = options.R_o;
        count = 1;
        nFreqs = numel(options.omega);
        Error = zeros(noRuns,nFreqs);
        N_arr = 0:noRuns-1;
        data = e3Dss(v, options);
        p = data(1).p;
        for N = N_arr
            N
            options.N_max = N;
            data = e3Dss(v, options);
            p_N = data(1).p;
            Error(N+1,:) = norm2((p - p_N).', 1)./norm2(p.', 1);
        end
%         keyboard
        N_max = min(find(max(Error,[],2) == 0));
        N_arr = 0:N_max;
        Error = Error(1:N_max+1,:);
        
        k = options.omega/options.c_f(1);
        [NN,kk] = meshgrid(N_arr,k);
        B = log10(Error);
        B(B == -Inf) = Inf;
        B(B == Inf) = min(min(B))+min(min(B))/1000;
        surf(kk*R_o,NN,B.','EdgeColor','none')
%         axis image
        
        set(0,'defaulttextinterpreter','latex')
        xlabel('$$k_1R_{0,1}$$')
        ylabel('$$N$$')
        zlabel('$$\|p-p_N\|_2/\|p\|_2$$')
        set(gca, 'Color', 'none')
        set(gca, 'Layer', 'top')
        box on
        
        grid off
        xlim([0 round(max(k*R_o))])
        ylim([min(N_arr) max(N_arr)])
        view(0,90)
        colorbar

        extraAxisOptions = {...
            'axis on top=true', ...
            'at={(0,0)}', ...
            'colorbar style={ylabel={$\frac{\left\|p_1-p_1^{(N)}\right\|_2}{\|p_1\|_2}$}, ytick={-18,-16,...,0}, yticklabels={$10^{-18}$, $10^{-16}$, $10^{-14}$, $10^{-12}$, $10^{-10}$, $10^{-8}$, $10^{-6}$, $10^{-4}$, $10^{-2}$, $10^0$}}'};
        
        colormap gray
        map = colormap;
        map = flipud(map);
        colormap(map)
        matlab2tikz([fileName '_bw.tex'], 'height', '3.2094in', ...
            'extraAxisOptions', extraAxisOptions)
        
        colormap default
        map = colormap;
        map = [1,1,1; map];
        colormap(map)
        matlab2tikz([fileName '.tex'], 'height', '3.2094in', 'width', '3.2094in', ...
            'extraAxisOptions', extraAxisOptions)
        
        title('Error plots')
end
        
