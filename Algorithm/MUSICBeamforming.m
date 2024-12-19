function [y_music_db, theta] = MUSICBeamforming(rec_signal, frequency, sound_speed, sensor_distance, sensor_number, theta_start, theta_end, d_theta, dynamic_range)
    % MUSIC算法波束形成成像函数
    % 基本参数设置
    wavelength = sound_speed / frequency;
    d = (0:sensor_number-1) * sensor_distance / wavelength;
    theta = theta_start:d_theta:theta_end;
    [~, n_samples] = size(rec_signal);
    n_sources = 1; % 声源数量
    % 初始化MUSIC谱
    P_music = zeros(length(theta), n_samples);
    % 对每个时间点（距离点）进行MUSIC处理
    for sample_idx = 1:n_samples
        % 计算协方差矩阵 - 直接使用矩阵乘法
        X = rec_signal(:, sample_idx);
        Rxx = X * X';
        % 特征值分解 - 完全按照参考代码处理
        [EV, D] = eig(Rxx);
        EVA = diag(D)';
        [EVA, I] = sort(EVA);
        EV = fliplr(EV(:, I));
        % 获取噪声子空间
        En = EV(:, n_sources+1:end);
        Un = En * En';  % 预计算投影矩阵
        % 计算所有角度的MUSIC谱
        for theta_idx = 1:length(theta)
            a = exp(-1j * 2 * pi * d * sin(theta(theta_idx) * pi/180)).';
            P_music(theta_idx, sample_idx) = 1 / real(a' * Un * a);
        end
    end
    
    % 转换为dB
    y_music_db = trans2db(abs(P_music), dynamic_range);
end

function output = trans2db(input, threshold)
    % 本函数将input进行归一化，并转化为dB值。略去小于threshold (dB)的值。
    % 输出output为input对应的dB值。
    if ~isreal(input)
        warning('输入序列应为实数');
    end
    if nargin == 1  %输入变量的个数
        threshold = -40;
    end
    [rows, columns] = size(input);
    if rows == 1 || columns == 1
        output = 10 * log10(input./max(input));
        output = max(output,threshold);
    else
        output = 10 * log10(input./max(max(input)));
        output = max(output,threshold);
    end
end
