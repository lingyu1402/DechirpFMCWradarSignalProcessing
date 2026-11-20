function proj2_rangedoppler()
%PROJECT2_RANGEDOPPLER 项目2：2D-FFT（Range-Doppler 距离-多普勒图）
%
%   快时间 FFT → 距离；慢时间 FFT → 多普勒。二者级联得到 R-D 图。
%
%   三个场景：
%     1) 静止人体 + 呼吸/心跳微动  → 躯干在 0 多普勒，微动边带被淹没
%     2) 匀速走动（不模糊）        → 明显多普勒峰 + 距离迁移
%     3) 超速走动（多普勒模糊）    → 速度折叠到负轴
%
%   【本脚本要回答的核心问题】
%   Q：生命体征为什么通常不走 R-D，而是在固定距离 bin 取慢时间相位？
%   A：呼吸 0.25 Hz / 心跳 1.2 Hz 对应的多普勒只有
%       v = λ·f/2 ≈ 0.0016 / 0.0075 m/s，远小于速度分辨率
%       Δv = λ/(2·T_CPI) ≈ 0.0012 m/s，几乎落在 0 频 bin 内；
%       同时躯干静态反射在 0 Hz 有强峰，其旁瓣把微动边带盖住。
%       所以正确做法是：CFAR 定 bin → 该 bin 慢时间相位 → unwrap/detrend → 谱分析。
%       R-D 用来「粗定位运动目标」，不是用来测呼吸心跳。
%
%   【公式辨析】
%       Δv   = λ/(2·T_CPI)      速度「分辨率」（不同速度能否分开）
%       vmax = λ·PRF/4          最大「不模糊」速度（超过就折叠）
%       盲速（MTI 语境）= n·λ·PRF/2（被对消掉的速度）
%       三个概念不同，别把 λ/(2·T_CPI) 叫盲速。
%
%   运行：addpath(文件所在目录); proj2_rangedoppler

    close all; clc;
    p = radar_param();

    Nf = p.Nfast;   % 500  快时间
    Ns = p.Nslow;   % 1024 慢时间

    fprintf('=== 项目2：2D-FFT Range-Doppler ===\n');
    fprintf('ΔR = %.3f m, Δv(分辨率) = %.5f m/s, vmax(不模糊) = %.3f m/s\n', ...
            p.dR, p.vres, p.vmax);
    fprintf('呼吸 %.2f Hz → v=%.5f m/s; 心跳 %.2f Hz → v=%.5f m/s\n\n', ...
            p.f_resp, p.lambda*p.f_resp/2, p.f_heart, p.lambda*p.f_heart/2);

    % 轴：距离 0 ~ Nf/2 个 bin；速度（fftshift 后从负到正）
    Raxis = (0:Nf/2-1)*p.dR;
    vaxis = (-Ns/2:Ns/2-1)*(p.Fslow/Ns)*p.lambda/2;

    %% ============ 场景1：静止人体 + 呼吸心跳 ============
    cube1 = gen_cube(p,Nf,Ns,3.0,0,4e-4,1e-4);
    [RF1,RD1] = rd_process(cube1,Nf,Ns);

    fh = figure('Name','场景1：R-D 静止+微动','NumberTitle','off');
    plot_rd(fh,vaxis,Raxis,RD1,'场景1：静止人体 + 呼吸/心跳微动（0 多普勒附近）');

    % 固定距离 bin 取慢时间相位 → 这才是生命体征的正确链路
    bin_h = round(3.0/p.dR);
    z = RF1(bin_h,:);
    phi = unwrap(angle(z));
    phi = phi - polyval(polyfit((0:Ns-1),phi,1),(0:Ns-1));   % 线性去趋势
    faxis = (0:Ns/2-1)*p.Fslow/Ns;
    figure('Name','场景1：慢时间相位谱','NumberTitle','off');
    subplot(2,1,1);
        plot((0:Ns-1)/p.Fslow,phi); grid on;
        xlabel('慢时间 (s)'); ylabel('相位 (rad)'); title('固定距离 bin 的去趋势相位');
    subplot(2,1,2);
        nshow = round(3*Ns/p.Fslow);
        Xphi  = abs(fft(phi))/Ns*2;
        plot(faxis(1:nshow),Xphi(1:nshow),'b','LineWidth',1); hold on;
        plot([p.f_resp p.f_resp],ylim,'r--');
        plot([p.f_heart p.f_heart],ylim,'g--');
        grid on; xlabel('频率 (Hz)'); ylabel('幅度');
        legend('相位谱','呼吸 0.25 Hz','心跳 1.2 Hz');
        title('生命体征在相位谱里可见，在 R-D 图上不可见');

    %% ============ 场景2：匀速走动 v = 0.4 m/s（不模糊）============
    v2 = 0.4;
    cube2 = gen_cube(p,Nf,Ns,3.0,v2,4e-4,1e-4);
    [RF2,RD2] = rd_process(cube2,Nf,Ns);

    fh = figure('Name','场景2：R-D 走动','NumberTitle','off');
    plot_rd(fh,vaxis,Raxis,RD2,sprintf('场景2：走动 v=%.2f m/s（< vmax，不模糊）',v2));
    fh = figure('Name','场景2：距离迁移','NumberTitle','off');
    plot_slowtime(fh,RF2,p,'场景2：走动目标的距离-慢时间轨迹（距离迁移）');

    %% ============ 场景3：超速走动 v = 0.9 m/s（多普勒模糊）============
    v3 = 0.9;
    cube3 = gen_cube(p,Nf,Ns,3.0,v3,4e-4,1e-4);
    [RF3,RD3] = rd_process(cube3,Nf,Ns);

    v_period = p.lambda*p.Fslow/2;             % 速度混叠周期 = 2·vmax = λ·PRF/2
    v_fold = v3 - round(v3/v_period)*v_period; % 折叠回 [-vmax, vmax]
    fh = figure('Name','场景3：R-D 模糊','NumberTitle','off');
    plot_rd(fh,vaxis,Raxis,RD3,...
        sprintf('场景3：真实 v=%.2f m/s，折叠显示为 %.2f m/s（模糊）',v3,v_fold));

    fprintf('场景2：真实 0.40 m/s → R-D 峰值位置在 +0.40 m/s 附近 ✓\n');
    fprintf('场景3：真实 %.2f m/s → R-D 峰值位置在 %.2f m/s（错！这就是多普勒模糊）\n',v3,v_fold);
end

%% ==================== 本地函数 ====================

function cube = gen_cube(p,Nf,Ns,R0,v,a_resp,a_heart)
% 生成 FMCW 数据立方体 cube(Nf,Ns)：行=快时间(距离)，列=慢时间(chirp 序号)
% 第 m 个 chirp 的目标距离：
%   R(m) = R0 + v·t_slow(m) + A_resp·sin(2πf_resp t) + A_heart·sin(2πf_heart t)
% 中频相位包含两项：
%   快时间项 exp(j2π·fbeat·t_fast)  → 决定「在哪个距离 bin」
%   慢时间项 exp(j4πR/λ)            → 决定「多普勒频率 = 2v/λ」
    t_fast = (0:Nf-1)/p.Fs;
    t_slow = (0:Ns-1)/p.Fslow;
    cube = zeros(Nf,Ns);
    for m = 1:Ns
        Rm = R0 + v*t_slow(m) ...
           + a_resp*sin(2*pi*p.f_resp*t_slow(m)) ...
           + a_heart*sin(2*pi*p.f_heart*t_slow(m));
        fbeat = 2*p.slope*Rm/p.c;
        phase = 4*pi*Rm/p.lambda;               % 往返路程 → 4πR/λ
        cube(:,m) = exp(1j*2*pi*fbeat*t_fast) * exp(1j*phase);
    end
    sigma = 0.2;
    cube = cube + sigma*(randn(Nf,Ns)+1j*randn(Nf,Ns));
end

function [RF,RD] = rd_process(cube,Nf,Ns)
% 2D-FFT：快时间 FFT（加汉宁窗降旁瓣）→ 慢时间 FFT（加窗 + fftshift）
    RF = fft(bsxfun(@times,cube,hann_win(Nf)),Nf,1);
    RF = RF(1:Nf/2,:);                          % 只取正距离半谱
    RD = fft(bsxfun(@times,RF,hann_win(Ns).'),Ns,2);
    RD = fftshift(RD,2);                        % 零多普勒移到中间
end

function plot_rd(fh,vaxis,Raxis,RD,tag)
    figure(fh);                                 % 在调用方指定的窗口里画，不再新建
    clf(fh);
    img = 20*log10(abs(RD)+eps);
    img = img - max(img(:));            % 归一化到峰值 0 dB：FFT 未按 N 缩放，绝对量级会把色标全顶成红色
    imagesc(vaxis,Raxis,img); axis xy;
    xlabel('速度 (m/s)'); ylabel('距离 (m)');
    title(tag,'FontSize',11);
    c = colorbar; c.Label.String = '幅度 (dB)';
    caxis([-55 0]); colormap jet;
end

function plot_slowtime(fh,RF,p,tag)
    figure(fh);                                 % 同上：不新建窗口
    clf(fh);
    img = 20*log10(abs(RF)+eps);
    img = img - max(img(:));            % 同上：归一化到峰值 0 dB
    imagesc((0:size(RF,2)-1)/p.Fslow,(0:size(RF,1)-1)*p.dR,img);
    axis xy; xlabel('慢时间 (s)'); ylabel('距离 (m)');
    title(tag,'FontSize',11);
    c = colorbar; c.Label.String = '幅度 (dB)'; caxis([-45 0]); colormap jet;
end

function w = hann_win(N)
% 本地汉宁窗，避免依赖 Signal Processing Toolbox
    if N == 1, w = 1; return; end
    w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
end
