function p = radar_param()
%RADAR_PARAM 人体毫米波雷达仿真统一参数
%
%   第一阶段全部仿真脚本共用此文件。集中管理，改一处全局生效。
%   面向「人体生命体征 / 跌倒检测」场景（非车载），24 GHz FMCW。
%
%   用法（在每个主脚本开头）：
%       p = radar_param();
%
%   注意：所有脚本请把本文件放在同一目录，或 addpath 该目录。

    % ---------- 射频与波长 ----------
    p.fc     = 24e9;             % 载频 (Hz)，24 GHz ISM 频段，人体雷达常用
    p.c      = 3e8;              % 光速 (m/s)
    p.lambda = p.c / p.fc;       % 波长 (m) = 1.25 cm

    % ---------- 快时间（距离维）----------
    p.B      = 1e9;              % 扫频带宽 (Hz)，1 GHz
    p.Tc     = 1e-3;             % 单个 chirp 扫频时长 (s)
    p.slope  = p.B / p.Tc;       % 调频斜率 (Hz/s) = 1e12
    p.Fs     = 500e3;            % ADC 采样率 (Hz)
    p.Nfast  = round(p.Fs*p.Tc); % 每 chirp 快时间采样点数 = 500

    % ---------- 慢时间（多普勒 / 生命体征维）----------
    p.Tchirp = 5e-3;             % chirp 重复间隔 (s)，含发射空闲 → PRF = 200 Hz
    p.Fslow  = 1/p.Tchirp;       % 慢时间采样率 (Hz) = 200
    p.Nslow  = 1024;             % 一个 CPI 内的 chirp 数
    p.Tcpi   = p.Nslow*p.Tchirp; % 相干处理间隔 (s) = 5.12 s

    % ---------- 派生量（务必分清概念）----------
    p.dR     = p.c/(2*p.B);              % 距离分辨率 (m) = 0.15 m
    p.Rmax   = p.Fs*p.c/(2*p.slope);     % 最大无模糊距离 (m) = 75 m
    p.vres   = p.lambda/(2*p.Tcpi);      % 速度「分辨率」Δv (m/s) ≈ 0.00122
                                         %   ← 注意：这是分辨率，不是盲速！
    p.vmax   = p.lambda*p.Fslow/4;       % 最大「不模糊」速度 (m/s) = 0.625
                                         %   ← 多普勒模糊 / 盲速看这个

    % ---------- 生命体征典型值 ----------
    p.f_resp  = 0.25;            % 呼吸率 15 次/分 (Hz)
    p.f_heart = 1.20;            % 心率 72 次/分 (Hz)
end
