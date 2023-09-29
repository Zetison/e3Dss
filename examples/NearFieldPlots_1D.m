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
plotInTimeDomain = 0;
noSteps = 30;
P_inc = 1;
type = 1;
% model = 'Skelton1997tao';
model = 'S15';
% model = 'S1';
% applyLoad = 'planeWave';
applyLoad = 'surfExcitation';
% applyLoad = 'radialPulsation';
switch model
    case 'S1'
        layer = setS1Parameters();
        f_arr = 1e3; % frequencies in Hertz
        BC = 'SSBC';
    case 'S15'
        layer = setS15Parameters();
        f_arr = 1e3; % frequencies in Hertz
        BC = 'SHBC';
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
        
end

R = layer{1}.R;
R_a = 1.5*R;
P_inc = 1;
theta_s = NaN(1,2);
r_s = 2*R; 
if strcmp(applyLoad,'pointCharge')
    if intermediatePointCharge
        r_s = layer{2}.R*1/3 + layer{3}.R*2/3;
    end
    P_inc = P_inc*r_s;
elseif strcmp(applyLoad,'surfExcitation')
    r_s = layer{1}.R;
    theta_s = [40,60]*pi/180;
elseif strcmp(applyLoad,'mechExcitation')
    r_s = layer{1}.R; 
end
if plotInTimeDomain
    f_c = 1500;
    omega_c = 2*pi*f_c;
    T = 120/f_c;
    % T = 1;
    N = 2^10;
%     N = 2^5;
    B = N/T; % bandwidth
    dt = T/N;
    f_L = -B/2;
    f_R = B/2;
    df = 1/T;
    f_arr = linspace(0,f_R-df,N/2);
end
    
layer{1}.calc_p = true;
layer{1}.calc_p_0 = true;
layer{1}.calc_p_inc = true;
layer{2}.calc_u_r = true;
layer{2}.calc_u_t = true;
layer{2}.calc_sigma_s = [1,0,0,0,0,0];
layer{3}.calc_p = true;
layer{3}.calc_p_inc = true;
if strcmp(model,'S1') || strcmp(model,'S15')
    layer = defineBCstring(layer,BC);
end
theta_s = [40,60]*pi/180;
omega = 2*pi*f_arr;
k = omega./layer{1}.c_f;
options = struct('BC', BC, ...
                 'd_vec', [0,0,1].', ...
                 'r_s', layer{1}.R, ...
                 'theta_s', theta_s, ...
                 'omega', omega, ...
                 'Display','iter', ...
                 'P_inc', P_inc, ...
                 'applyLoad',applyLoad);
if plotInTimeDomain
    options.omega = omega(2:end);
end


switch plotType
    case 'radial'
        beta_f = pi/2;
        alpha_f = 0;
        r_arr{1} = linspace(-layer{1}.R,-2*layer{1}.R,npts).';
        r_arr{2} = linspace(-layer{2}.R,-layer{1}.R,npts).';
        r_arr{3} = linspace(-layer{3}.R,-layer{2}.R,npts).';
        X = cell(1,3);
        for i = 1:numel(layer)
            X{i} = r_arr{i}*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];
        end
    case 'angular'
        beta_f = linspace(pi/2,-pi/2,npts).';
%         beta_f = linspace(-pi/4,-pi/2,npts).';
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
if plotInTimeDomain
    options.P_inc = @(omega) P_inc_(omega,omega_c,P_inc,type);
end
if 1
    layer = e3Dss(layer, options);
    save([resultsFolder, '/layer_' model '_' BC])
else
    load([resultsFolder, '/layer' model '_' BC],'layer')
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
if plotInTimeDomain
    totField = cell(1,3);
    totField{1} = zeros(npts,N/2);
    totField{2} = zeros(npts,N/2);
    totField{3} = zeros(npts,N/2);
    PincField = zeros(npts,N/2);

    for n = 0:N-1
        f = f_L + (f_R-f_L)/N*n;
        omega = 2*pi*f;
        if n >= N/2+1
            for m = 1:numel(layer)
                if strcmp(layer{m}.media,'fluid')
                    totField{m}(:,n-N/2+1) = layer{m}.p(:,n-N/2) + layer{m}.p_inc(:,n-N/2);
                    if any(layer{m}.p_inc ~= 0)
                        PincField(:,n-N/2+1) = layer{m}.p_inc(:,n-N/2);
                    end
                else
                    totField{m}(:,n-N/2+1) = layer{m}.sigma_s{1}(:,n-N/2);
                end
            end
        end
    end
    totFieldTime = cell(1,3);
    for m = 1:numel(layer)
        totFieldTime{m} = 2/T*fft(totField{m},N,2);
    end
    PincFieldTime = 2/T*fft(PincField,N,2);

    startIdx = 2*round(1900/2);
    startIdx = round(N*0.9);
    if N < 100
        startIdx = 1;
    end
    for i = 1:3
        temp = totFieldTime{i};
        totFieldTime{i}(:,1:N-startIdx+1) = temp(:,startIdx:end);
        totFieldTime{i}(:,N-startIdx+2:end) = temp(:,1:startIdx-1);
    end
    temp = PincFieldTime;
    PincFieldTime(:,1:N-startIdx+1) = temp(:,startIdx:end);
    PincFieldTime(:,N-startIdx+2:end) = temp(:,1:startIdx-1);

    figure(3)
    m_arr = 0:N-1;
    for m = m_arr
        t = dt*m;
        switch plotType
            case 'radial'
                plot(r_arr{3},real(totFieldTime{3}(:,m+1)),'blue','DisplayName','p_{tot}');
                hold on
                plot(r_arr{2},-real(totFieldTime{2}(:,m+1)),'blue','DisplayName','-\sigma_{rr}');
                plot(r_arr{1},real(totFieldTime{1}(:,m+1)),'blue','DisplayName','p_{tot}');
                if strcmp(applyLoad,'pointCharge')
                    plot(r_arr{3},real(PincFieldTime(:,m+1)),'red','DisplayName','p_{inc}');
                elseif strcmp(applyLoad,'radialPulsation') || strcmp(applyLoad,'planeWave')
                    plot(r_arr{1},real(PincFieldTime(:,m+1)),'red','DisplayName','p_{inc}');
                end
                ylim([-2,2])
                xlim([r_arr{1}(end),r_arr{3}(1)])
                plot(-layer{1}.R*[1,1],ylim,'--','color','black')
                plot(-layer{2}.R*[1,1],ylim,'--','color','black')
                plot(-layer{3}.R*[1,1],ylim,'--','color','black')
                hold off
            %     legend show % This is really slow...
                drawnow
            case 'angular'
                plot(theta*180/pi,real(totFieldTime{1}(:,m+1)),'blue','DisplayName','p_{tot}');
                ylim([-2,2])
                xlim([0,180])
        end
        titleStr = sprintf('Time step %4d, t = %5.3fs. T = %5.3fs.',m,t,T);
        title(titleStr)
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
            xlim([r_arr{1}(end),r_arr{end}(1)])
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
