function [y_cbf_db, theta] = RLCBFBeamforming(received_signal, frequency, sound_speed, ...
    sensor_distance, sensor_number, angle_start, angle_end, d_theta, db_floor)
    
    % 计算波数k
    k = 2 * pi * frequency / sound_speed;
    
    % 生成传感器位置向量n
    n = (0:sensor_number-1);
    
    % 生成角度向量theta和对应的sine值
    theta = angle_start:d_theta:angle_end;
    sine_theta = sind(theta);
    
    % 计算导向矢量
    steering_vector_1 = (1/sensor_number) * exp(-1i*k*n.'*sensor_distance*sine_theta);
    
    % 计算常规波束形成结果
    y_cbf = steering_vector_1' * received_signal;
    
    % 计算点扩散函数矩阵
    N = length(sine_theta);
    Bp_list = zeros(N);
    Yp1 = zeros(N);
    
    for loop = 1:N
        Yp = steering_vector_1' * exp(-1i*k*n.'*sensor_distance*sine_theta(loop));
        Yp1(loop, :) = Yp.';
        Bp_list(loop, :) = (abs(Yp).^2).';
    end
    
    % 准备RL解卷积输入
    y_cbf_11 = (abs(y_cbf).^2).';
    
    % 对每个距离进行RL解卷积
    decbf = zeros(size(y_cbf_11));
    for i_RL = 1:size(y_cbf_11, 1)
        decbf(i_RL,:) = mydeconv(y_cbf_11(i_RL,:), Bp_list, 100);
    end
    
    % 转换为dB尺度
    y_cbf_db = trans2db(decbf, db_floor);
end

function output = mydeconv(B, Bp_list, time)
    % RL解卷积函数
    % B: 输入空间谱
    % Bp_list: 点扩散函数矩阵
    % time: 迭代次数
    
    % 计算点扩散矩阵每列和
    M = sum(Bp_list, 1);
    Bp_list_M = Bp_list ./ M;
    
    % 初始化迭代值为CBF空间谱
    S = B;
    
    % 迭代求解
    for Loop = 1:time
        P_estimate = (Bp_list * S.').';
        R = B ./ P_estimate;
        temp = R * Bp_list_M;
        S = S .* temp;
    end
    
    output = S;
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