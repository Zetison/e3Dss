% close all
clear all %#ok
% close all
clf
% 
startup
resultsFolder = [folderName '/NearFieldPlots_1D'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end
%% Test benchmark models
npts = 1000;
% npts = 10;
plotType = 'angular';
% plotType = 'radial';
polarPlot = false;
generateVideo = 0;
plotTimeOscillations = 0;
noSteps = 30;
P_inc = 1;
model = 'Skelton1997tao';
% model = 'S15';
% model = 'S1';
applyLoad = 'planeWave';
% applyLoad = 'radialPulsation';
switch model
    case 'S1'
        layer = setS1Parameters();
        f_arr = 1e3; % frequencies in Hertz
    case 'S15'
        layer = setS15Parameters();
        f_arr = 1e3; % frequencies in Hertz
    case 'Skelton1997tao'
        layer = setSkelton1997taoParameters();
        withCoating = false;
        if ~withCoating
            layerSSBC = layer([1,3]);
            layerSSBC{1}.R = layer{2}.R;
            layer = layerSSBC;
        end
        BC = 'SSBC';
        layer{1}.calc_p_0 = true; % Calculate the far field pattern
        layer{1}.calc_p = true;
        % npts_f = 2;
        load('miscellaneous/Skelton_extremas')
        npts_f = 2000;
%         npts_f = 10;
        f_max = 25e3;
%         f_max = 25e2;
        f_arr = linspace(f_max/npts_f,f_max,npts_f);
        f_arr = sort(unique([f_arr, specialValues(specialValues <= f_max).']));
        
        f_max = 25e4;
        f_arr = [f_arr, linspace(f_arr(end)+(f_max-f_arr(end))/npts_f,f_max,npts_f)];
        f_arr = sort(unique([f_arr, specialValues(specialValues <= f_max).']));
        
        layer{2}.calc_u_r = true;
        layer{2}.calc_u_t = true;
        layer{2}.calc_sigma_s = [1,0,0,0,0,0];
end
if strcmp(model,'S1') || strcmp(model,'S15')
    layer{2}.calc_sigma_s = [1,0,0,0,0,0];
    layer{1}.calc_p = true;
    layer{1}.calc_p_inc = true;
    layer{3}.calc_p = true;
    BC = 'SSBC';
    layer = defineBCstring(layer,BC);
    for i = 1:numel(layer)
        layer{i}.X = X{i};
    end
end
omega = 2*pi*f_arr;
k = omega./layer{1}.c_f;
options = struct('BC', BC, ...
                 'd_vec', [0,0,1].', ...
                 'omega', omega, ...
                 'Display','iter', ...
                 'P_inc', P_inc, ...
                 'applyLoad',applyLoad);


switch plotType
    case 'radial'
        beta_f = pi/2;
        alpha_f = 0;
        r_arr{1} = linspace(layer{1}.R,2*layer{1}.R,npts).';
        r_arr{2} = linspace(layer{2}.R,layer{1}.R,npts).';
        r_arr{3} = linspace(layer{3}.R,layer{2}.R,npts).';
        X = cell(1,3);
        for i = 1:numel(layer)
            X{i} = r_arr{i}*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];
        end
    case 'angular'
        beta_f = linspace(pi/2,-pi/2,npts).';
        beta_f = linspace(-pi/4,-pi/2,npts).';
%             beta_f = beta_f(1:end-10);
        theta = pi/2-beta_f;
        alpha_f = 0;
        X{1} = layer{1}.R*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];
        if polarPlot
            X{2} = X{1};
            axpolar(1) = subplot(2,4,5,polaraxes);
            axpolar(2) = subplot(2,4,6,polaraxes);
            axpolar(3) = subplot(2,4,7,polaraxes);
            axpolar(4) = subplot(2,4,8,polaraxes);
        else
            X{2} = [X{1}; layer{2}.R*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))]];
        end
        X{3} = [];
end
for i = 1:numel(layer)
    layer{i}.X = X{i};
end
if 0
    layer = e3Dss(layer, options);
    save([resultsFolder, '/layer'])
else
    load([resultsFolder, '/layer'],'layer')
end
    
if polarPlot
    rmax(1) = max(20*log10(abs(layer{2}.u_r(:))));
    rmax(2) = max(20*log10(abs(layer{2}.u_t(:))));
    rmax(3) = max(20*log10(abs(layer{2}.sigma_s{1}(:))));
    rmax(4) = max(20*log10(abs(layer{1}.p(:))));
    rmin(1) = min(20*log10(abs(layer{2}.u_r(:))));
    rmin(2) = min(20*log10(abs(layer{2}.u_t(:))));
    rmin(3) = min(20*log10(abs(layer{2}.sigma_s{1}(:))));
    rmin(4) = min(20*log10(abs(layer{1}.p(:))));
else
    rmax(1) = max(real(layer{2}.u_r(:)));
    rmin(1) = min(real(layer{2}.u_r(:)));
end

TS = 20*log10(abs(layer{1}.p_0));
h = subplot(2,4,[1,4]);
ax2 = subplot(2,4,[5,8]);
xlabel(h,'$f$','interpreter','latex')
ylabel(h,'TS','interpreter','latex')

if plotTimeOscillations
    load('miscellaneous/Skelton_extremas','specialValues')
    specialValues = sort(specialValues);
    indices = [13:20,63,200,201,400,431,530];
%     indices = [431,530];
    analyze_f = specialValues(indices).';
    if generateVideo
        vObj = VideoWriter([resultsFolder '/' model '_oscillations.avi']);
        vObj.FrameRate = 15;
        open(vObj);
    end
    for f = analyze_f
        i_f = find(f==f_arr);
        omega = 2*pi*f;
        i = 0;

        u_r_inner = layer{2}.u_r(1:npts,i_f);
        u_r_outer = layer{2}.u_r(npts+1:end,i_f);
        max_u_r = max(max(abs(u_r_inner)),max(abs(u_r_outer)));
        while i < 3*noSteps
            plot(h,f_arr,TS(end,:))
            hold(h,'on')
            plot(h,f,TS(end,i_f),'o','color','red')
            hold(h,'off')
            
            t = i/noSteps*2*pi/omega;
            plot(ax2,theta*180/pi,real(u_r_inner*exp(-1i*omega*t)),'DisplayName','Outer surface')
            hold(ax2,'on')
            plot(ax2,theta*180/pi,real(u_r_outer*exp(-1i*omega*t)),'DisplayName','Inner surface')
            h.YLim = [-30,10];
            hold(ax2,'off')
            ax2.YLim = [-max_u_r,max_u_r];
            ax2.XLim = [min(theta*180/pi),max(theta*180/pi)];
            xlabel(ax2,'$\theta$','interpreter','latex')
            ylabel(ax2,'$\mathrm{Re}(u_r e^{i\omega t})$','interpreter','latex')
            t_string = sprintf('f = %0.3fkHz, $t = %.2f$ms',f/1e3,1e3*t);
            title(ax2,t_string,'interpret','latex');
            legend(ax2,'show')
            pause(0.001)
            if generateVideo
                frame = getframe(gcf);
                writeVideo(vObj,frame);
            end
            i = i + 1;
        end
        pause(0.5)
    end
    if generateVideo
        close(vObj);
    end
    return
end

if generateVideo
    vObj = VideoWriter([resultsFolder '/' model '_animation.avi']);
    vObj.FrameRate = 15; 
    open(vObj);
end
for i_f = 1:numel(f_arr)
    f = f_arr(i_f);
    plot(h,f_arr,TS(end,:))
    hold(h,'on')
    plot(h,f,TS(end,i_f),'o','color','red')
    hold(h,'off')
    h.YLim = [-30,10];
    switch plotType
        case 'radial'
            switch numel(layer)
                case 1
                    plot(r_arr{1}, real(layer{1}.p_inc + layer{1}.p))
                case 2
                    plot([r_arr{2}; r_arr{1}], [-real(layer{2}.sigma_s{1}); real(layer{1}.p_inc + layer{1}.p)])
                otherwise
                    plot([r_arr{3}; r_arr{2}; r_arr{1}], [real(layer{3}.p); -real(layer{2}.sigma_s{1}); real(layer{1}.p_inc + layer{1}.p)])
            end
        case 'angular'
            if polarPlot
                polarplot(axpolar(1), theta, 20*log10(abs(layer{2}.u_r(:,i_f))))
                polarplot(axpolar(2), theta, 20*log10(abs(layer{2}.u_t(:,i_f))))
                polarplot(axpolar(3), theta, 20*log10(abs(layer{2}.sigma_s{1}(:,i_f))))
                polarplot(axpolar(4), theta, 20*log10(abs(layer{1}.p(:,i_f))))
            else
                u_r_inner = real(layer{2}.u_r(1:npts,i_f));
                u_r_outer = real(layer{2}.u_r(npts+1:end,i_f));
                max_u_r = max(max(abs(u_r_inner)),max(abs(u_r_outer)));
                plot(ax2,theta*180/pi,u_r_inner,'DisplayName','Outer surface')
                hold(ax2,'on')
                plot(ax2,theta*180/pi,u_r_outer,'DisplayName','Inner surface')
                ax2.YLim = [-max_u_r,max_u_r];
                
                hold(ax2,'off')
            end
    end
%     legend show
    switch plotType
        case 'radial'
            xlim([-90,90])
        case 'angular'
            if polarPlot
                for i = 1:numel(axpolar)
                    axpolar(i).ThetaZeroLocation = 'top';
                    axpolar(i).ThetaDir = 'clockwise'; % 90 degrees at the right
                    axpolar(i).ThetaLim = [0,180];
                    axpolar(i).RLim = [rmin(i),rmax(i)];
                end
                title(axpolar(1),'$20\log_{10}(|u_r|)$','interpret','latex');
                title(axpolar(2),'$20\log_{10}(|u_t|)$','interpret','latex');
                title(axpolar(3),'$20\log_{10}(|\sigma_{rr}|)$','interpret','latex');
                title(axpolar(4),'$20\log_{10}(|p|)$','interpret','latex');
            else
%                 ax2.YLim = [rmin(1),rmax(1)];
%                 ax2.YLim = [-1e-11,1e-11];
                ax2.XLim = [min(theta*180/pi),max(theta*180/pi)];
                xlabel(ax2,'$\theta$','interpreter','latex')
                ylabel(ax2,'$\mathrm{Re}(u_r)$','interpreter','latex')
                legend(ax2,'show')
            end
    end
    pause(0.001)
%     return
    if generateVideo
        frame = getframe(gcf);
        writeVideo(vObj,frame);
    end
end
if generateVideo
    close(vObj);
end
