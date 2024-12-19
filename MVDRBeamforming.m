function [y_mvdr_db, theta] = MVDRBeamforming(received_signal, frequency, sound_speed, ...
                                          sensor_distance, sensor_number, ...
                                          angle_start, angle_end, d_theta, db_floor)
    % MVDRBeamforming - 最小方差无失真响应(MVDR)波束形成算法实现
    %
    % 输入参数:
    %   received_signal - 接收信号矩阵 [sensor_number × samples]
    %   frequency - 信号频率 (Hz)
    %   sound_speed - 声速 (m/s)
    %   sensor_distance - 阵元间距 (m)
    %   sensor_number - 阵元数量
    %   angle_start - 起始角度 (度)
    %   angle_end - 终止角度 (度)
    %   d_theta - 角度分辨率 (度)
    %   db_floor - 分贝下限值
    %
    % 输出参数:
    %   y_mvdr_db - 波束形成后的分贝值输出
    %   theta - 对应的角度向量
    % 生成角度向量
    theta = angle_start:d_theta:angle_end;
    sine_theta = sind(theta);
    % 生成阵元序号向量
    n = 0:(sensor_number-1);
    % 计算波数
    k = 2 * pi * frequency / sound_speed;
    % 计算导向矢量矩阵
    steering_matrix = (1 / sqrt(sensor_number)) * ...
        exp(-1i * k * n.' * sensor_distance * sine_theta);
    % 计算协方差矩阵
    R = received_signal * received_signal' / size(received_signal, 2);
    % 对协方差矩阵进行对角加载以提高稳定性
    diagonal_loading = 1e-3 * trace(R) / sensor_number;
    R = R + diagonal_loading * eye(sensor_number);
    % MVDR波束形成
    y_mvdr = zeros(length(theta), size(received_signal, 2));
    for i = 1:length(theta)
        a = steering_matrix(:, i);
        % 计算MVDR权值向量
        % 计算 R 的逆矩阵
        R_inv = inv(R);
        % 计算分子：R^(-1)a(θ)
        numerator = R_inv * a;
        % 计算分母：a(θ)^H * R^(-1)^H * a(θ)
        denominator = (a)' * R_inv' * a;
        % 计算最终结果
        w = numerator / denominator;
        % 波束形成
        y_mvdr(i, :) = w' * received_signal;
    end
    % 转换为分贝值
    y_mvdr_db = trans2db(abs(y_mvdr).^2, db_floor);
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
