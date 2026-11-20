function proj1_cfar()
%PROJECT1_CFAR 项目1：CFAR 恒虚警检测仿真
%
%   目标：在距离谱上自动「找到人在哪个 bin」，而不是对所有 bin 盲算相位。
%   这是生命体征链路和跌倒检测的前置步骤。
%
%   实现两种最经典的检测器：
%     CA-CFAR  （单元平均）  —— 均匀噪声里最优，杂波边缘/多目标会失效
%     OS-CFAR  （有序统计）  —— 用第 k 大代替均值，抗边缘、抗多目标掩蔽
%
%   本脚本分三部分：
%     场景A  均匀噪声 + 单目标       → 看阈值怎么形成、虚警怎么控制
%     场景B  杂波边缘 + 双目标       → CA / OS 的失效与稳健对比
%     校验   蒙特卡洛                → 实测虚警率是否等于设计 Pfa
%
%   运行：addpath(文件所在目录); proj1_cfar

    close all; clc;
    p = radar_param();   %#ok<NASGU>  % 统一参数（本脚本主要用 dR 标距离轴）

    % ---------- 检测器参数 ----------
    Nr   = 100;             % 距离 bin 数 → 0~15 m
    Nref = 16;              % 参考单元总数（左右各 8）
    Ng   = 2;               % 单侧保护单元数（防目标能量泄漏进参考窗）
    Pfa  = 1e-4;            % 设计虚警概率
    k_os = round(0.75*Nref);% OS 取第 k 小（= 较大的那个，偏保守）

    Raxis = (0:Nr-1)*p.dR;

    fprintf('=== 项目1：CFAR 恒虚警检测 ===\n');
    fprintf('Nref=%d, 单侧保护=%d, Pfa=%.0e, OS 的 k=%d\n', Nref, Ng, Pfa, k_os);
    fprintf('CA 缩放因子 alpha = %.4f\n', Nref*(Pfa^(-1/Nref)-1));
    fprintf('OS 缩放因子 alpha = %.4f\n\n', os_alpha(Nref,k_os,Pfa));

    %% ============ 场景A：均匀噪声 + 单目标 ============
    rng(7);
    mu = 1;                                  % 噪声平均功率
    Pa = exp_pow(Nr,mu);                     % 距离维功率：复高斯噪声 → 指数分布
    bin_t = round(3/p.dR);                   % 3 m 处的人 → bin 20
    Pa(bin_t) = Pa(bin_t) + 15;              % 注入目标，SNR = 15 (约 11.8 dB)

    [det_ca_a, thr_ca_a] = ca_cfar(Pa,Nref,Ng,Pfa);
    [det_os_a, thr_os_a] = os_cfar(Pa,Nref,Ng,Pfa,k_os);

    fh = figure('Name','场景A：单目标','NumberTitle','off');
    plot_profile(fh, Raxis, Pa, thr_ca_a, det_ca_a, 'CA-CFAR：均匀噪声 + 单目标');
    fh = figure('Name','场景A：OS','NumberTitle','off');
    plot_profile(fh, Raxis, Pa, thr_os_a, det_os_a, 'OS-CFAR：均匀噪声 + 单目标');

    fprintf('场景A：目标在 %.2f m (bin %d)\n', Raxis(bin_t), bin_t);
    fprintf('  CA 检出：%s\n', binstr(find(det_ca_a)));
    fprintf('  OS 检出：%s\n\n', binstr(find(det_os_a)));

    %% ============ 场景B：杂波边缘 + 双目标 ============
    rng(11);
    Pb = exp_pow(Nr,mu);
    Pb(50:90) = exp_pow(41,mu*8);            % bin 50~90 强杂波平台（墙/家具）
    Pb(30) = Pb(30) + 30;                    % 低杂波区目标
    Pb(70) = Pb(70) + 120;                   % 强杂波区目标
    Pb(76) = Pb(76) + 150;                   % 距大目标 6 bin：落在参考窗内但避开保护单元 → 演示「掩蔽效应」

    [det_ca_b, thr_ca_b] = ca_cfar(Pb,Nref,Ng,Pfa);
    [det_os_b, thr_os_b] = os_cfar(Pb,Nref,Ng,Pfa,k_os);

    fh = figure('Name','场景B：CA 边缘失效','NumberTitle','off');
    plot_profile(fh, Raxis, Pb, thr_ca_b, det_ca_b, 'CA-CFAR：杂波边缘效应 + 多目标掩蔽');
    draw_edge(fh, Raxis([50 90]));

    fh = figure('Name','场景B：OS 稳健','NumberTitle','off');
    plot_profile(fh, Raxis, Pb, thr_os_b, det_os_b, 'OS-CFAR：杂波边缘更稳健');
    draw_edge(fh, Raxis([50 90]));

    fprintf('场景B：低杂波目标 bin 30，强杂波目标 bin 70，掩蔽目标 bin 76，杂波边缘 bin 50/90\n');
    fprintf('  CA 检出：%s\n', binstr(find(det_ca_b)));
    fprintf('  OS 检出：%s\n\n', binstr(find(det_os_b)));

    %% ============ 蒙特卡洛：虚警率校验 ============
    M = 400;
    fa_ca = 0; fa_os = 0; ncut = 0;
    for m = 1:M
        Pn = exp_pow(Nr,mu);                             % 纯噪声，无目标
        [d1,t1] = ca_cfar(Pn,Nref,Ng,Pfa);
        [d2,t2] = os_cfar(Pn,Nref,Ng,Pfa,k_os);
        valid = ~isnan(t1);                              % 边缘 CUT 不参与统计
        fa_ca = fa_ca + sum(d1 & valid);
        fa_os = fa_os + sum(d2 & valid);
        ncut  = ncut  + sum(valid);
    end
    fprintf('蒙特卡洛 (%d 次)：CA 实测 Pfa = %.2e，OS 实测 Pfa = %.2e，设计值 = %.0e\n', ...
            M, fa_ca/ncut, fa_os/ncut, Pfa);
end

%% ==================== 本地函数 ====================

function P = exp_pow(n,mu)
% 生成 n 个独立指数分布功率样本，均值 mu
% 物理来源：复高斯噪声 n = a + jb，|n|^2 服从指数分布
    P = -log(rand(n,1))*mu;
end

function [det,thr] = ca_cfar(P,Nref,Ng,Pfa)
% CA-CFAR：参考单元求均值 → 乘缩放因子得阈值
% 推导：ref 单元功率 iid 指数(均值 mu)，Z=mean(ref)
%       Pfa = E[exp(-T/mu)] = (1+alpha)^(-Nref)
%       其中 T = alpha*Z  ⇒  alpha = Nref*(Pfa^(-1/Nref)-1)
    N = numel(P); half = floor(Nref/2);
    alpha = Nref*(Pfa^(-1/Nref)-1);
    thr = nan(N,1);
    for i = 1:N
        idx = [(i-Ng-half):(i-Ng-1), (i+Ng+1):(i+Ng+half)];
        idx = idx(idx>=1 & idx<=N);
        if numel(idx) < Nref, continue; end    % 边缘参考窗不足，跳过
        thr(i) = alpha*mean(P(idx));
    end
    det = P > thr;                             % NaN 比较自然为 false
end

function [det,thr] = os_cfar(P,Nref,Ng,Pfa,k)
% OS-CFAR：参考单元升序排序，取第 k 个作为杂波功率估计
% k 偏大 → 阈值高、虚警低、抗多目标掩蔽好，但检测损失略增
    N = numel(P); half = floor(Nref/2);
    alpha = os_alpha(Nref,k,Pfa);
    thr = nan(N,1);
    for i = 1:N
        idx = [(i-Ng-half):(i-Ng-1), (i+Ng+1):(i+Ng+half)];
        idx = idx(idx>=1 & idx<=N);
        if numel(idx) < Nref, continue; end
        Z = sort(P(idx),'ascend');
        thr(i) = alpha*Z(k);
    end
    det = P > thr;
end

function a = os_alpha(N,k,Pfa)
% 用精确 OS 虚警公式 + 二分法求缩放因子 a
% 模型：N 个 iid Exp(均值1)，X_(k) 为第 k 小，T = a*X_(k)
%   Pfa = P(x > a X_(k)) = N!/(N-k)! * Gamma(a+N-k+1)/Gamma(a+N+1)
% 用 gammaln 避免大数溢出。a 越大 → 阈值越高 → Pfa 越小。
    lo = 0; hi = 1e4;
    for it = 1:80
        a = (lo+hi)/2;
        lg = gammaln(N+1)-gammaln(N-k+1) ...
           + gammaln(a+N-k+1)-gammaln(a+N+1);
        if exp(lg) > Pfa
            lo = a;
        else
            hi = a;
        end
    end
    a = (lo+hi)/2;
end

function plot_profile(fh,R,P,thr,det,tag)
% 画距离功率谱 + CFAR 阈值 + 检测标记（统一 dB 坐标）
    figure(fh);
    dBp = 10*log10(P+eps);
    plot(R,dBp,'b','LineWidth',1.1); hold on;
    plot(R,10*log10(thr+eps),'r--','LineWidth',1.3);
    ii = find(det);
    if ~isempty(ii)
        plot(R(ii),dBp(ii),'ro','MarkerFaceColor','r','MarkerSize',7);
    end
    xlabel('距离 (m)'); ylabel('功率 / 阈值 (dB)');
    title(tag,'FontSize',11); legend('距离谱','CFAR 阈值','检测点','Location','best');
    grid on; ylim([min(dBp)-5 max(dBp)+5]);
end

function draw_edge(fh,xedge)
% 在已有图上叠加杂波边缘竖线（不依赖 xline，兼容旧版 MATLAB）
    figure(fh);
    yl = ylim;
    for x = xedge
        plot([x x],yl,'k--','LineWidth',1.2);
    end
    legend('距离谱','CFAR 阈值','检测点','杂波边缘','Location','best');
end

function s = binstr(ii)
% 把检测 bin 列表格式化为字符串，便于日志阅读
    if isempty(ii)
        s = '无';
    else
        s = sprintf('%d ',ii);
    end
end
