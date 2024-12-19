function [y_cbf_db, theta] = CBFBeamforming(received_signal, frequency, sound_speed, ...
                                          sensor_distance, sensor_number, ...
                                          angle_start, angle_end, d_theta, db_floor)
    % CBFBeamforming - 常规波束形成(CBF)算法实现
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
    %   y_cbf_db - 波束形成后的分贝值输出
    %   theta - 对应的角度向量
    
    % 生成角度向量
    theta = angle_start:d_theta:angle_end;
    sine_theta = sind(theta);
    
    % 生成阵元序号向量
    n = 0:(sensor_number-1);
    
    % 计算波数
    k = 2 * pi * frequency / sound_speed;
    
    % 计算导向矢量
    steering_vector = (1 / sensor_number) * ...
        exp(-1i * k * n.' * sensor_distance * sine_theta);
    
    % 波束形成
    y_cbf = steering_vector' * received_signal;
    
    % 转换为分贝值
    y_cbf_db = trans2db(abs(y_cbf).^2, db_floor);
end

