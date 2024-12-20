function [y_damas_db, theta] = DAMASBeamforming(received_signal, frequency, sound_speed, ...
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
    for loop = 1:N
        Yp = steering_vector_1' * exp(-1i*k*n.'*sensor_distance*sine_theta(loop));
        Bp_list(loop, :) = (abs(Yp).^2).';
    end
    
    % 准备DAMAS输入
    y_cbf_power = (abs(y_cbf).^2).';
    
    % 对每个距离进行DAMAS处理
    y_damas = zeros(size(y_cbf_power));
    for i_dist = 1:size(y_cbf_power, 1)
        y_damas(i_dist,:) = damas_deconv(y_cbf_power(i_dist,:), Bp_list, 5);
    end
    
    % 转换为dB尺度
    y_damas_db = trans2db(y_damas, db_floor);
end

function output = damas_deconv(Y, A, maxIter)
    % DAMAS解卷积函数
    % Y: 输入CBF结果
    % A: 点扩散函数矩阵
    % maxIter: 最大迭代次数
    
    % 初始化
    [N_angle, ~] = size(A);
    Q = zeros(1, N_angle);
    Q0 = real(Y);
    deps = 0.1;  % 收敛阈值
    
    % Gauss-Seidel迭代
    for iter = 1:maxIter
        for n = 1:N_angle
            sum_prev = sum(A(n, 1:n-1) .* Q(1:n-1));
            sum_next = sum(A(n, n+1:end) .* Q0(n+1:end));
            Q(n) = max(0, Y(n) - sum_prev - sum_next);
        end
        
        % 检查收敛性
        dX = (Q - Q0);
        maxd = max(abs(dX(:)))/mean(Q0(:));
        
        if maxd < deps
            break;
        end
        
        Q0 = Q;
    end
    
    output = Q;
end

function output = trans2db(input, threshold)
    % 本函数将input进行归一化，并转化为dB值。略去小于threshold (dB)的值。
    % 输出output为input对应的dB值。
    if ~isreal(input)
        warning('输入序列应为实数');
    end
    if nargin == 1
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