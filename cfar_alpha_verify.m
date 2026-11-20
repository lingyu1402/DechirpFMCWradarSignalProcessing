function cfar_alpha_verify()
%CFAR_ALPHA_VERIFY CFAR 系数 alpha 公式验证（CA 闭式解 / OS 数值解）
%
%   配合 CFAR_推导文档.html 使用：文档里所有理论值和蒙特卡洛实测值
%   都由本脚本生成，保证可复现。
%
%   验证内容：
%     1) CA：alpha = N*(Pfa^(-1/N)-1)，回代 Pfa=(1+alpha/N)^(-N)
%     2) OS：alpha 无闭式解，用精确 Gamma 公式 + 二分法求，回代验证
%     3) 蒙特卡洛：两种算法在纯噪声下的实测虚警率
%     4) 近似公式 alpha≈k*(Pfa^(-1/k)-1) 的误差评估
%     5) 检测概率 Pd 公式验证（SNR=10 dB）
%
%   中文注释。运行：cfar_alpha_verify

    clc;
    rng(2026);

    N = 16;               % 参考单元总数
    k = round(0.75*N);    % OS 取第 k 小 = 12
    mu = 1;               % 噪声平均功率（指数分布均值）
    M = 4e5;              % 蒙特卡洛次数
    snr_db = 10;          % 检测概率验证用的 SNR

    fprintf('=== CFAR alpha 验证 (N=%d, k=%d, M=%d, mu=%g) ===\n\n', N, k, M, mu);

    % 一次生成、多次复用，保证不同 Pfa 的对比公平
    ref  = -log(rand(M,N))*mu;         % M 次实验的 N 个参考单元功率，指数分布
    cut0 = -log(rand(M,1))*mu;         % H0：纯噪声 CUT
    snr  = 10^(snr_db/10);
    cut1 = -log(rand(M,1))*(mu*(1+snr)); % H1：SNR=10dB 的 CUT

    Zca = mean(ref,2);                  % CA 的杂波功率估计
    sref = sort(ref,2);                 % OS 需要排序
    Zos = sref(:,k);                    % OS 的杂波功率估计 = 第 k 小

    designs = [1e-2, 1e-3, 1e-4];

    fprintf('%-8s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
            'Pfa', 'CA_alpha', 'CA_thy', 'CA_mc', 'OS_alpha', 'OS_thy', ...
            'OS_mc', 'OSapprox');
    fprintf('%s\n', repmat('-',1,108));

    for p = designs
        % ---- CA：闭式解 ----
        a_ca = N*(p^(-1/N) - 1);
        p_ca_thy = (1 + a_ca/N)^(-N);              % 回代理论值
        p_ca_mc  = mean(cut0 > a_ca*Zca);          % 蒙特卡洛

        % ---- OS：精确 Gamma 公式 + 二分 ----
        a_os = os_alpha_exact(N,k,p);
        p_os_thy = os_pfa(N,k,a_os);               % 回代理论值
        p_os_mc  = mean(cut0 > a_os*Zos);          % 蒙特卡洛

        % ---- OS 近似公式 ----
        a_os_ap = k*(p^(-1/k) - 1);
        p_os_ap_mc = mean(cut0 > a_os_ap*Zos);

        fprintf('%-8.0e %-12.4f %-12.3e %-12.3e %-12.4f %-12.3e %-12.3e %-12.3e\n', ...
                p, a_ca, p_ca_thy, p_ca_mc, a_os, p_os_thy, p_os_mc, p_os_ap_mc);
    end

    % ---- 检测概率验证 ----
    fprintf('\n=== 检测概率 Pd (SNR=%d dB) ===\n', snr_db);
    for p = [1e-3, 1e-4]
        a_ca = N*(p^(-1/N) - 1);
        pd_ca_thy = (1 + a_ca/(N*(1+snr)))^(-N);
        pd_ca_mc  = mean(cut1 > a_ca*Zca);

        a_os = os_alpha_exact(N,k,p);
        beta = a_os/(1+snr);
        pd_os_thy = os_pfa(N,k,beta);              % 同构：alpha -> alpha/(1+SNR)
        pd_os_mc  = mean(cut1 > a_os*Zos);

        fprintf('Pfa=%.0e: CA Pd thy=%.4f mc=%.4f | OS Pd thy=%.4f mc=%.4f\n', ...
                p, pd_ca_thy, pd_ca_mc, pd_os_thy, pd_os_mc);
    end

    % ---- 数值表：给文档用的 alpha 速查 ----
    fprintf('\n=== alpha 速查表 (N=%d, k=%d) ===\n', N, k);
    fprintf('%-10s %-12s %-12s %-12s\n', 'Pfa', 'CA_alpha', 'OS_alpha', 'OS近似');
    for p = [1e-2,1e-3,1e-4,1e-5,1e-6]
        a_ca = N*(p^(-1/N) - 1);
        a_os = os_alpha_exact(N,k,p);
        a_ap = k*(p^(-1/k) - 1);
        fprintf('%-10.0e %-12.4f %-12.4f %-12.4f\n', p, a_ca, a_os, a_ap);
    end
end

%% ==================== 本地函数 ====================

function pfa = os_pfa(N,k,a)
% OS-CFAR 精确虚警概率：
%   Pfa = [N!/(N-k)!] * Gamma(a+N-k+1) / Gamma(a+N+1)
% 用 gammaln 避免阶乘/Gamma 溢出
    lg = gammaln(N+1) - gammaln(N-k+1) ...
       + gammaln(a+N-k+1) - gammaln(a+N+1);
    pfa = exp(lg);
end

function a = os_alpha_exact(N,k,Pfa)
% 二分法求 OS 缩放因子：Pfa(a) 随 a 单调递减
    lo = 0; hi = 1e4;
    for it = 1:100
        a = (lo+hi)/2;
        if os_pfa(N,k,a) > Pfa
            lo = a;
        else
            hi = a;
        end
    end
    a = (lo+hi)/2;
end
