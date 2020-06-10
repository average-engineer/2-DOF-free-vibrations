clear all
close all
clc

%system parameters
%damping of front wheel
c_f = 100;
%stiffness of front wheel
k_f = 60000;
%damping of rear wheel
c_r = 100;
%stiffness of rear wheel
k_r = 60000;
%mass
m = 1000;
%mass moment of intertia about COM
j = 1000;
%front wheel offset from COM [m]
l_f = 2.5;
%rear wheel offset from COM [m] 
l_r = 2.5;

%initial conditions
%bounce
x_0 = 0.1;
x_dot_0 = 0;
%pitch
p_0 = 1.0;%radians
p_dot_0 = 0;%radians

%sampling rate
fs = 100;

%time span
time_span = [0:1/fs:50];

%type of solver
solver = 'anal';

%number of test runs for each computation
n = 10;

switch solver  
    case 'anal'
        for kk = 1:n
            tic
        %system matrices
        %mass matrix
        M = [m,0;0,j];
        %stiffness matrix
        K = [(k_r + k_f),(k_f*l_f - k_r*l_r);(k_f*l_f - k_r*l_r),(k_f*l_f^2 + k_r*l_r^2)];
        %damping matrix
        C = [(c_r + c_f),(c_f*l_f - c_r*l_r);(c_f*l_f - c_r*l_r),(c_f*l_f^2 + c_r*l_r^2)];
        
        %using polyeigen function to find out the eigenvalues and
        %eigenvectoes of the system
        %in this system, the mass matrix is associated with lambda^2,
        %stiffness matrix is associated with lambda^0 and damping matrix is
        %associated with lambda
        [eigen_vector,lambda] = polyeig(K,C,M);
        
        %the eigen values
        %for this particular case, eigen values are coming out as conjugate
        %complex pairs with negative real parts and thus the system is
        %comparitively weakly damped
        l_1 = lambda(1);
        l_2 = lambda(2);
        l_3 = lambda(3);
        l_4 = lambda(4);
        
        %eigen vector elements
        %first eigen vector
        c_11 = eigen_vector(1,1);
        c_12 = eigen_vector(2,1);
        %second eigenvector
        c_21 = eigen_vector(1,2);
        c_22 = eigen_vector(2,2);
        %3rd eigen vector
        c_31 = eigen_vector(1,3);
        c_32 = eigen_vector(2,3);
        %4th eigenvector
        c_41 = eigen_vector(1,4);
        c_42 = eigen_vector(2,4);
        
        %initial condition matrix
        w_0 = [x_0;p_0;x_dot_0;p_dot_0];
        %intermediate eigen vector matrix
        A = [c_11,c_21,c_31,c_41;c_12,c_22,c_32,c_42;l_1*c_11,l_2*c_21,l_3*c_31,l_4*c_41;l_1*c_12,l_2*c_22,l_3*c_32,l_4*c_42];
        %integration constants matrix
        int_const = A\w_0;
        
        %EOM for bounce
        x_t_anal = real(c_11*int_const(1)*exp(l_1*time_span) + c_21*int_const(2)*exp(l_2*time_span) + c_31*int_const(3)*exp(l_3*time_span) + c_41*int_const(4)*exp(l_3*time_span));
        x_t_anal = x_t_anal';
        v_t_anal = real(l_1*c_11*int_const(1)*exp(l_1*time_span) + l_2*c_21*int_const(2)*exp(l_2*time_span) + l_3*c_31*int_const(3)*exp(l_3*time_span) + l_4*c_41*int_const(4)*exp(l_4*time_span));
        %EOM for pitch
        p_t_anal = real(c_12*int_const(1)*exp(l_1*time_span) + c_22*int_const(2)*exp(l_2*time_span) + c_32*int_const(3)*exp(l_3*time_span) + c_42*int_const(4)*exp(l_4*time_span));
        p_t_anal = p_t_anal';
        vp_t_anal = real(l_1*c_12*int_const(1)*exp(l_1*time_span) + l_2*c_22*int_const(2)*exp(l_2*time_span) + l_3*c_32*int_const(3)*exp(l_3*time_span) + l_4*c_42*int_const(4)*exp(l_4*time_span));
        %computation time for each test run/loop
        computational_time_anal(kk) = toc;
        %plotting
        figure(1)
subplot(1,2,1)
hold on
aa = plot(time_span,x_t_anal,'-*','color','k')
bb = plot(time_span,v_t_anal,'-*','color','r')
legend([aa,bb],'Bounce Motion','Bounce Velocity')
title('Bounce Motion')
subplot(1,2,2)
hold on
cc = plot(time_span,p_t_anal,'-*','color','k')
dd = plot(time_span,vp_t_anal,'-*','color','r')
legend([cc,dd],'Pitch Angular Motion','Pitch Angular Velocity') 
title('Pitch Motion')
        

end
%average computation time for analytical solution
anal_time = sum(computational_time_anal)/n;

end

%angular eigen frequencies
%only works for cases of damping = 0 or damping << spring stiffness 
eigenf_angular = abs(lambda);
%eigen frequencies
eigenf = eigenf_angular/(2*pi);

% %discrete fourier transform
% %number of points for fourier transform
% % N = length(x_t_anal);
% % 
% % %sampling rate
% % fs = 100;
% % 
% % %we can only analyse frequencies which are half of the sampling rate
% % freq = [0:fs/N:fs/2];
% % 
% % tic
% % %the frequency/fourier coefficient matrix
% % y = fft(x_t_anal);
% % DFT_clocker = toc
% % 
% % %half of the elements of the frequency coefficient matirx are taken as the
% % %frequency is limited to half of the sampling rate
% % y_half = y(1:floor(N/2) + 1);
% % 
% % %for the amplitude, we have to divide each element of the frequency
% % %coefficient matrix by the number of points taken
% % y_half = y_half/length(x_t_anal);
% % 
% % %amplitudes of the signal
% % amp = abs(y_half);
% % %amp(2:end) = 2*amp(2:end);
% 
% %dft for bounce motion
% %number of data points
% N1 = length(x_t_anal);
% 
% %fundamental frequency
% omega1 = exp(-i*(2*pi)/N1);
% for iii = 1:N1
%     for jjj = 1:N1
%         DFT1(iii,jjj) = omega1^((iii-1)*(jjj-1));
%     end
% end
% 
% freq_coeffs = DFT1*x_t_anal;
% freq_coeffs = freq_coeffs(1:floor(length(freq_coeffs)/2) + 1);
% %frequency vector
% freq = [0:fs/N1:fs/2];
% 
% %amplitude
% amp = abs(freq_coeffs/length(x_t_anal));
% amp(2:end) = 2*amp(2:end);
% 
% figure(2)
% plot(freq,amp)
% xlabel('Frequency,[Hz]')
% ylabel('Amplitude,[m]')
% title('1st Eigen Frequency')
% 
% %dft for pitch motion
% %number of data points
% N2 = length(p_t_anal);
% 
% %fundamental frequency
% omega2 = exp(-i*(2*pi)/N2);
% for iii = 1:N2
%     for jjj = 1:N2
%         DFT2(iii,jjj) = omega2^((iii-1)*(jjj-1));
%     end
% end
% 
% freq_coeffs1 = DFT2*p_t_anal;
% freq_coeffs1 = freq_coeffs1(1:floor(length(freq_coeffs1)/2) + 1);
% %frequency vector
% freq1 = 0:fs/N2:fs/2;
% 
% %amplitude
% amp1 = abs(freq_coeffs1/length(p_t_anal));
% amp1(2:end) = 2*amp1(2:end);
% 
% figure(3)
% plot(freq1,amp1)
% xlabel('Frequency,[Hz]')
% ylabel('Amplitude,[m]')
% title('2nd Eigen Frequency')


%Fast Fourier Transform
%bounce motion FFT
%here the number of points (or the number of elements of the data vector
%should be greater than 2^u (or length(x_t))
u = log2(length(x_t_anal));
N = 2^ceil(u);

%frequency vector
freq_fft = [0:fs/N:fs/2];

%applying fft algorithm
tic
freq_coeffs_fft = fft(x_t_anal,N);%matlab extends the datapoint vector (x_t) by adding zero elements
FFT_timer = toc

%the frequency coefficient vector is cut in half because the frequency
%range has exactly half the values 
freq_coeffs_fft = freq_coeffs_fft(1:floor(length(freq_coeffs_fft)/2) + 1);

%amplitude
amp_fft = abs(freq_coeffs_fft/length(x_t_anal));
amp_fft(2:end) = 2*amp_fft(2:end);

figure(2)
plot(freq_fft,amp_fft)
title('Bounce Motion FFT')
xlabel('Frequency,[Hz]')
ylabel('Amplitude,[m]')

%pitch motion FFT
%here the number of points (or the number of elements of the data vector
%should be greater than 2^u (or length(x_t))
u1 = log2(length(x_t_anal));
N1 = 2^ceil(u1);

%frequency vector
freq_fft1 = [0:fs/N1:fs/2];

%applying fft algorithm
tic
freq_coeffs_fft1 = fft(p_t_anal,N1);
FFT_timer1 = toc

%the frequency coefficient vector is cut in half because the frequency
%range has exactly half the values 
freq_coeffs_fft1 = freq_coeffs_fft1(1:floor(length(freq_coeffs_fft1)/2) + 1);

%amplitude
amp_fft1 = abs(freq_coeffs_fft1/length(p_t_anal));
amp_fft1(2:end) = 2*amp_fft1(2:end);

figure(3)
plot(freq_fft1,amp_fft1)
title('Pitch Motion FFT')
xlabel('Frequency,[Hz]')
ylabel('Amplitude,[rad]')
