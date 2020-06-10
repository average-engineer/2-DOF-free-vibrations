function [cost] = cost_function(damping,solver)
%system parameters
%damping of front wheel
c_f = damping;
%stiffness of front wheel
k_f = 60000;
%damping of rear wheel
c_r = damping;
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


%number of test runs for each computation
n = 10;

switch solver
   
    case 'ode45'
        for jj = 1:n
            tic
        %initial conditions vector
        w_0 = [x_0;p_0;x_dot_0;p_dot_0];
        [t,results] = ode45(@(t,w)statefunction(w,t,m,j,k_f,k_r,c_f,c_r,l_f,l_r),time_span,w_0);
        computational_time_num(jj) = toc;
        %bounce motion
        x_t_num = results(:,1);
        %pitch motion
        p_t_num = results(:,2);
        %bounce velocity
        v_t_num = results(:,3);
        %pitch velocity
        vp_t_num = results(:,4);
        
        %plotting the motions
        end
 
%average computation time for numerical solution
num_time = sum(computational_time_num)/n;

%cost function for numerical solution
%cost function is basically the absolute bounce and pitch motion of the
%mass
cost = sum(abs(x_t_num))*0.4 + sum(abs(p_t_num))*0.6;
        
        
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
        

end
%average computation time for analytical solution
anal_time = sum(computational_time_anal)/n;

%cost function for analytical solution
cost = sum(abs(x_t_anal))*0.4 + sum(abs(p_t_anal))*0.6;

end
end