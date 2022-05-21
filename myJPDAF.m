clc;clear all;close all;

%%��������
T = 20;                                                            %��������
dt = 1;                                                            %TΪ�������
MC_number = 1;                                                     %Monte Carlo�������
target_num = 2;                                                    %Ŀ�����
target_initposition = [1500 500; 
                       300 400; 
                       500 1500; 
                       400 300];                              %Ŀ�����ʼλ�ú��ٶ�(x,vx,y,vy)         
Pd = 1;                                                            %������
Pg = 0.99;                                                         %��ȷ��������������ڵø���
g_sigma = 9.21;                                                    %����,����ֵ����R��Q�йأ������������
lambda = 2;
gamma = lambda*10^(-6);                                            %ÿһ����λ���(km^2)�ڲ���lambda���Ӳ� 

P = zeros(4, 4, target_num, T);                                       %Э�������
P_init = diag([10000 0 10000 0]);        %������ʼЭ������� ��x,vx,y,vy�ĸ������Э���������ʵֵ�͹���ֵ�Ĺ�������Э�������
for i = 1 : target_num
    P(:, :, i, 1) = P_init;
end
F = [1 dt 0 0;                                                           
     0 1 0 0;
     0 0 1 dt;
     0 0 0 1];                                             %״̬ת�ƾ���
Q = diag([4 0 4 0]);                                       %ϵͳ��������Э����
H = [1 0 0 0;
     0 0 1 0];                                             %�۲����
R = diag([100 100]);                                       %�۲�Э�������

X_est_out = zeros(4, target_num, T);                                   %�洢Ŀ��ĸ�ʱ�̵��˲�ֵ
X_est_out_MC = zeros(4, target_num, T, MC_number);                    %MC_number��Montle Carlo��������ȫ������洢����
Z_out = zeros(2, target_num, T);                           %�۲�ֵx,y
X_out = zeros(4, target_num, T);                                   %ʵ��״ֵ̬x,vx,y,vy   

%%������ʵ״̬��ģ������
X_out(:, :, 1) = target_initposition;                                          %ʵ��λ�þ����ʼ�� 
for i = 1 : target_num
    for t = 2 : T                                                                     %ʵ��λ�� 
        X_out(:, i, t) = (F * X_out(:, i, t-1)) + (sqrt(Q) * (randn(4, 1)));        
    end
end
%%��ͼ
X_target_1 = zeros(1, T);         
Y_target_1 = zeros(1, T);         
for t = 1 : T
    X_target_1(t) = X_out(1, 1, t);
    Y_target_1(t) = X_out(3, 1, t);
end
figure;
plot(X_target_1, Y_target_1, 'b-');
hold on;
X_target_2 = zeros(1, T);         
Y_target_2 = zeros(1, T);         
for t = 1:T
    X_target_2(t) = X_out(1, 2, t);
    Y_target_2(t) = X_out(3, 2, t);
end
plot(X_target_2, Y_target_2, 'r-');
xlabel('x(m)'),ylabel('y(m)');
legend('Ŀ��a��ʵ��λ��','Ŀ��b��ʵ��λ��');
grid;

%%��ʼ����
for m = 1 : MC_number    
    %%1.�����۲�ֵ
    for t = 1 : T %��������
        for i = 1 : target_num %Ŀ�������  %���������۲��λ��
            Z_out(:, i, t) = H * X_out(:, i, t) + sqrt(R) * randn(2, 1);
        end
    end
    %%2.�����Ӳ�,��ȷ����Ч�۲� 
    S = zeros(2, 2, target_num);
    Z_prediction = zeros(2, target_num);                                                               %�洢����Ŀ��Ĺ۲�Ԥ��ֵ,��ֻ����x,y����
    X_prediction = zeros(4, target_num);                                                               %�洢����Ŀ���״̬Ԥ��ֵ,������x,y�����x,y�����ٶ�
    ellipse_Volume = zeros(1, target_num);
    NOISE_target_1 = [];                                                                    %�洢Ŀ��1���Ӳ�
    NOISE_target_2 = [];                                                                    %�洢Ŀ��2���Ӳ�

    for t = 1 : T
        %�����Ӳ�����
        noise = [];               %�����Ӳ�����
        NOISE = [];               %����target���Ӳ�����
        for i = 1 : target_num      %����
            if t ~= 1
                X_prediction(:, i) = F * X_est_out(:, i, t-1);                                       %��ǰһʱ�̵��˲�ֵ��Ԥ�⵱ǰ��ֵ(kalman�˲��ĵ�һ�����ʽ)
                P_prediction = F * P(:, :, i, t-1) * F' + Q;
            else
                X_prediction(:, i) = target_initposition(:, i);                                       %��һ�β�����������ʵλ�õ�Ԥ��ֵ 
                P_prediction = P(:, :, i, 1);
            end
            Z_prediction(:, i) = H * X_prediction(:, i);                                                 %��ȡԤ��ֵ��x,y���꣬����x,y�ٶ�
            S(:, :, i) = H * P_prediction * H' + R;       
            ellipse_Volume(i) = pi * g_sigma * sqrt(det(S(:,:,i)));                              %ÿ��Ŀ�������Բ�����ŵ����   
            number_returns = floor(ellipse_Volume(i) * gamma + 1);                      %��Բ�������ڵĴ���ز���
            side = sqrt((ellipse_Volume(i) * gamma + 1) / gamma) / 2;                       %����Բ�����ŵ�ЧΪ�����Σ�����������α߳��Ķ���֮һ
            noise_x = X_prediction(1, i) + side - 2 * rand(1, number_returns) * side;         %��Ԥ��ֵ��Χ��������ز���
            noise_y = X_prediction(3, i) + side - 2 * rand(1, number_returns) * side;         %ע�⣺��ĳһ��number_returnsС�ڵ���0ʱ�����������һ�μ��ɡ�
            noise = [noise_x; noise_y];
            NOISE = [NOISE noise];
            if i == 1
                NOISE_target_1 = [NOISE_target_1 noise];
            else
                NOISE_target_2 = [NOISE_target_2 noise];
            end
        end
        current_measurements = [NOISE Z_out(:, 1, t) Z_out(:, 2, t)];  %��¼�˿̵Ĺ۲�ֵ�����Ӳ���Ŀ��1��Ŀ��2�Ĺ۲�ֵ��

        %%���ڲ�������ʼjpda�ļ��㣬ǰ��ֻ�ǲ���ʵ������Ҫ�����ݡ� 
        %%�����۲�ȷ�Ͼ���measure_confirm_matrix 
        valid_measurements = [];
        num_valid_measurements = 0;                                                                          %��¼��Ч�۲����
        [dim, num] = size(current_measurements); %dim��ʾ�۲�ά�ȣ�num��ʾ�۲�ֵ�ĸ���
        measure_confirm_matrix = zeros(1000, target_num + 1);
        for j = 1 : num %Ŀ����Ӳ�
            flag = 0;   %�۲�ֵ�Ƿ���Ч
            for i = 1 : target_num %Ŀ�꣬Ҳ������Ϊ�Ǵ��ڵĺ�����
                delta_measurement = current_measurements(:, j) - Z_prediction(:, i);  %����ֵ��Ԥ��ֵ�����ֵ��
                delta_measurement_cov = delta_measurement' / S(:, :, i) * delta_measurement; %ͨ��ÿ��Ŀ�����������Э������󣬵õ�y1(�Ӳ�����Ŀ��)��Z_predic(Ŀ��)���Ƴ̶ȡ�                      
                if delta_measurement_cov <= g_sigma                                                    
                    flag = 1;
                    measure_confirm_matrix(num_valid_measurements + 1, 1) = 1; %ȷ�Ͼ���ÿһ�е�һ��Ԫ�رض�Ϊ1��
                    measure_confirm_matrix(num_valid_measurements + 1, i+1) = 1; %��������޷�Χ������Ϊ1
                end
                %�Ƿ�Ӧ����else�����ڲ��ڲ��ŷ�Χ�ڵ�Ŀ�꣬�����µĺ��������������ʡ������º����ķ���������ʷ�켣�ʹ�����֣��Ƚϸ��ӣ��Ⱥ���
            end  %% Ҳ����������ѭ�������ȷ�Ͼ���ı�ը����Ϊ����һ֡ѭ����ʱ�򣬸�Ŀ��c��Ϊ��ǰ֡�õ������к���
            if flag == 1   %�ü��Ŀ�겻���Ǻ��ĸ�Ŀ��۲����ˣ�����¼������
                valid_measurements = [valid_measurements current_measurements(:, j)];                                           %������������е����лز�����y��
                num_valid_measurements = num_valid_measurements + 1;                                                            %��¼��Ч�۲����
            end
        end
        measure_confirm_matrix = measure_confirm_matrix(1 : num_valid_measurements, 1 : target_num + 1);

        %%4.������������interconnect_matrix,����num��ʾ���������¼�����
        interconnect_matrix = zeros(num_valid_measurements, 3, 10000);
        interconnect_matrix(:, 1, 1:10000) = 1;  %��һ�ж�����Ϊ1������ÿ��Ŀ�궼�п��������Ӳ���
        if num_valid_measurements ~= 0                                 %num_valid_measurements=0��ʾ����Ŀ�궼û�й۲⵽������������Ŀ����Ӳ���û�й۲⵽��
            event_num = 1; %���������Ҳ�����¼�����
            for i = 1 : num_valid_measurements
                if measure_confirm_matrix(i, 2) == 1 %��һ����������
                    interconnect_matrix(i, 2, event_num) = 1;
                    interconnect_matrix(i, 1, event_num) = 0;
                    event_num = event_num + 1;
                    for j = 1 : num_valid_measurements
                        if (i ~= j) & (measure_confirm_matrix(j, 3) == 1)
                            interconnect_matrix(i, 2, event_num) = 1;
                            interconnect_matrix(i, 1, event_num) = 0;
                            interconnect_matrix(j, 3, event_num) = 1;
                            interconnect_matrix(j, 1, event_num) = 0;
                            event_num = event_num + 1;
                        end
                    end
                end
            end                                   

            for i = 1 : num_valid_measurements
                if measure_confirm_matrix(i, 3) == 1
                    interconnect_matrix(i, 3, event_num) = 1;
                    interconnect_matrix(i, 1, event_num) = 0;
                    event_num = event_num + 1;
                end
            end
        else
            flag = 1;
        end
        interconnect_matrix = interconnect_matrix(:, :, 1:event_num);           %��ٷ���ֵĽ������interconnect_matrix��
        
        %%  5.����������Pr,����
        %%  num_false_measurements��ʾ������,
        %%  mea_indicator��ʾ�۲�ָʾ��,�Ƿ�ÿ������ֵ����Ŀ�����
        %%  target_indicator��ʾĿ��ָʾ�� ���Ƿ�ÿ��Ŀ�궼����⵽��
        %%  ����ÿ���¼������ĸ���
        Pr = zeros(1, event_num); %numΪ�¼�����
        for i = 1 : event_num 
            num_false_measurements = num_valid_measurements;
            N = 1;  %��⹫ʽ����֮һ
            for j = 1 : num_valid_measurements    %num_valid_measurementsΪ����Ч�����ڵ��Ӳ���Ŀ������
                mea_indicator = sum(interconnect_matrix(j, 2:3, i));                                      %�ο�������ʽ4-48
                if mea_indicator == 1 %��ʾ���������У��й����ϵ�Ŀ��
                    num_false_measurements = num_false_measurements - 1;
                    if interconnect_matrix(j, 2, i) == 1                                                  %����۲���Ŀ��1����
                        b = (valid_measurements(:, j) - Z_prediction(:, 1))' * inv(S(:, :, 1)) * (valid_measurements(:, j) - Z_prediction(:, 1));
                        N = N / sqrt(det(2 * pi * S(:,:,1))) * exp(-1/2 * b);                          %������̬�ֲ�����                         
                    else                                                                   %����۲���Ŀ��2����
                        b = (valid_measurements(:, j) - Z_prediction(:, 2))' * inv(S(:, :, 2)) * (valid_measurements(:, j) - Z_prediction(:, 2));
                        N = N / sqrt(det(2 * pi *S(:,:,2))) * exp(-1/2 * b);                          %������̬�ֲ�����                         
                    end                                                                        
                end
            end
            if Pd == 1
                a = 1;
            else
                a = 1;
                for j = 1 : target_num
                    target_indicator = sum(interconnect_matrix(:, j+1, i));                               %�ο�������ʽ4-49
                    a = a * Pd^target_indicator * (1 - Pd)^(1 - target_indicator);                   %���������
                end
            end                                                                            
            V = ellipse_Volume(1) + ellipse_Volume(2);                                         %��ʾ������������
            a1 = 1;
            for j = 1 : num_false_measurements
                a1 = a1*j;
            end
            Pr(i) = N * a * a1 / (V^num_false_measurements);   %mk��Ϊ�ǳ������ڹ�һ�������лᱻԼ���������ȥ�����������þ��ȷֲ���Ҳ��ȥ
        end
        Pr = Pr / sum(Pr);
        %%6.�����������U
        %%���Ϲ������ʼ�Ϊ���յ�ÿ��������Ŀ��֮������Ŀ�����ϵ����
        U = zeros(num_valid_measurements + 1, target_num);
        for i = 1 : target_num
            for j = 1 : num_valid_measurements
                for k = 1 : event_num
                    U(j, i) = U(j, i) + Pr(k) * interconnect_matrix(j, i+1, k);     %���ú������pr����ÿ����������ĸ��ʡ�
                end
            end
        end
        U(num_valid_measurements + 1, :) = 1 - sum(U(1:num_valid_measurements, 1:target_num));                %��������Ŀ��T�����Ĺ������ʴ���U��num_valid_measurements+1,:),��һ��
        %%7.Kalman�˲���ʼ
        for i = 1 : target_num    
            if t ~= 1                        
                P_prediction = F * P(:, :, i, t-1) * F' + Q;
            else
                P_prediction = P(:, :, i, 1);
            end                                                                      %����Э�������
            K(:, :, i) = P_prediction * H' * inv(S(:, :, i));
            P(:, :, i, t) = P_prediction - (1 - U(num_valid_measurements + 1, i)) * K(:, :, i) * S(:, :, i) * K(:, :, i)';
        end
        for i = 1 : target_num
            a = 0;         
            b = 0;
            x_est = 0;                                                                   %������õ��м����
            for j = 1 : num_valid_measurements
                x_est = x_est + U(j, i)*(X_prediction(:, i) + K(:, :, i)*(valid_measurements(:, j) - Z_prediction(:, i)));
            end
            x_est = U(j+1, i) * X_prediction(:, i) + x_est; 
            X_est_out(:, i, t) = x_est;
            for j = 1 : num_valid_measurements + 1
                if j == num_valid_measurements + 1
                    a = X_prediction(:, i);
                else
                    a = X_prediction(:, i) + K(:, :, i)*(valid_measurements(:, j) - Z_prediction(:, i));
                end
                b = b + U(j,i) * (a*a' - x_est*x_est');
            end
            P(:, :, i, t) = P(: ,: ,i ,t) + b; 
            X_est_out_MC(:, i, t, m) = X_est_out(:, i, t);   
        end
    end
end


%%��ͼ                                            %�˲�ֵ��ƽ��
%%1.�˲����
figure;
%Ŀ��a,b�Ĺ۲�λ��
for i = 1 : target_num
    measure_x = zeros(1, T);
    measure_y = zeros(1, T);
    for j = 1 : T
        measure_x(j) = Z_out(1, i, j);
        measure_y(j) = Z_out(2, i, j);
    end
    if i == 1
       plot(measure_x(:), measure_y(:), 'r:')
    else 
       plot(measure_x(:), measure_y(:), 'b:')
    end
    hold on;
end
%Ŀ��a,b���Ӳ�λ��
for i = 1 : target_num
    if i == 1
       plot(NOISE_target_1(1, :), NOISE_target_1(2, :), 'r.');
    else
       plot(NOISE_target_2(1, :), NOISE_target_2(2, :), 'b.');
   end
   hold on;
end
%Ŀ��a,b�Ĺ���λ��
for i = 1 : target_num
    estimate_x = zeros(1, T);
    estimate_y = zeros(1, T);
    for j = 1 : T
        estimate_x(j) = X_est_out(1, i, j);
        estimate_y(j) = X_est_out(3, i, j);
    end
    if i == 1
        plot(estimate_x(:), estimate_y(:), 'r-');
    else 
        plot(estimate_x(:), estimate_y(:), 'g-');
    end
    hold on;
end
xlabel('x/m'), ylabel('y/m');
legend('Ŀ��a�Ĺ۲�λ��', 'Ŀ��b�Ĺ۲�λ��', 'Ŀ��a���Ӳ�', 'Ŀ��b���Ӳ�', 'Ŀ��a�Ĺ���λ��', 'Ŀ��b�Ĺ���λ��');grid;

%%2.���
figure;
mse = zeros(target_num, T);
for i = 1 : target_num
    for t = 1 : T
        for M = 1 : MC_number                                                              %��С�������
            mse_t = (X_est_out_MC(1, i, t, M) - X_out(1, i, t))^2 + (X_est_out_MC(3, i, t, M) - X_out(3, i, t))^2;
            mse(i, t) = mse(i, t) + mse_t;
        end
        mse(i, j) = sqrt(mse(i, j) / MC_number);
    end
    %�����ֵ
    [max_num, max_index] = max(mse(i, :));
    str = strcat('\itT=', num2str(max_index), ' \itError=', num2str(max_num), '(m^2)');
    text(max_index, 0.8*max_num, str);
    hold on;
    if i == 1
        plot([max_index max_index],[0 max_num],'r');
        hold on;
        plot(1:T, mse(i, :), 'r:') 
        hold on;
    else 
        plot([max_index max_index],[0 max_num],'b');
        hold on;
        plot(1:T, mse(i, :), 'b:') 
        hold on;
    end
end
xlabel('times'),ylabel('����ֵ�����ֵ������/m^2');
legend('Ŀ��a��������ֵ','Ŀ��a�����','Ŀ��b��������ֵ','Ŀ��b�����');grid;

%rmse
for i = 1 : target_num    
    figure;
    rmse_est = zeros(1,T);
    rmse_measure = zeros(1,T);
    for t = 1 : T
        rmse_est(t) = sqrt((X_est_out(1, i, t) - X_out(1, i, t))^2 + (X_est_out(3, i, t) - X_out(3, i, t))^2);
        rmse_measure(t) = sqrt((Z_out(1, i, t) - X_out(1, i, t))^2 + (Z_out(2, i, t) - X_out(3, i, t))^2);
    end
    plot(1:T, rmse_est,'k:'); 
    if i == 1
        xlabel('times'), ylabel('RMSE of a');    
    else 
        xlabel('times'), ylabel('RMSE of b');
    end
    grid;
end
