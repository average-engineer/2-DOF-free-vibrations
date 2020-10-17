# 2-DOF-free-vibrations
Code for calculating vibration characteristics of 2 degree of freedom mechanical systems.
The 2 DOF system is assumed to be a simple car model with its mass concentrated in a rectangular mass which can translate vertically (bounce) and rotate about its central axis normal to the plane (pitch). The front and rear suspension of the car is modelled using a spring and damper system connected to the rectangle at a particular distance from the centre of mass.
There is no source of excitation for the model and the car has some initial bounce and pitch motion, thus the sytem will experience free vibrations which can be mathematically modelled by a 2nd order coupled (or uncoupled depending on the system parameters) homogeneous ordinary differential equation.
3 approaches have been taken in the code to solve for the system:

1). Analytical approach where the code is built according to the same exact steps which are taken while analytically solving the equation of motion (generating the characteristic polynomial and then calculating the eigen values from the polynomial). This approach has the fastest computing time, but is limited as it can be used only for simpler systems having shallow modelling depth.

2). Numerical approach where the equation of motion is represented in the state space vector form (`w_dot = A*w` where `w_dot` is the first derivative of the state vector, `w` is the state vector and `A` is the system matrix) and the `ode45` solver in MATLAB is used with the objective function whose independent variables are the time vector and the state vector while the dependent vector is the first derivate of the state vector. This approach is a little slower than the analytical approach but can be used with ease in systems having complex modelling depths.

3). Symbolic approach where the equations of motion (2nd order homogeneous ODE) is solved symbollically. This approach is much slower than the other two, and while can be used in complex problems, the excess time taken will not be worth it.

Explanation of the scripts included in the repository:
`two_dof_system_free_oscillations.m` is the main script which computes the bounce and pitch amplitudes and velocities using all the three approaches, analytical, numerical and symbolic. For the numerical solver, `statefunction.m` is the function script which describes the state space representation of the system. 

`fourier_transform_two_dof_free.m` is a script which applies the fast-fourier transform to the system to convert the time domain (amplitude vs time) to the frequency domain (amplitude vs frequency) which provides us with the information about the eigen-modes of the system. In the amplitude vs frequency plot, the eigen modes are represented by peaks in the plot. Since, this is a 2 DOF system, there will be 2 eigen modes, both of which are present in form of peaks for both the bounce and pitch motions vs frequency plots if the system is coupled. For this script, the analytical approach of solving the system is used since we can find the eigen frequencies of the system using the analytical method and can validate the eigenfrequencies calculated by the FFT method.

`two_dof_parameter_fit.m` is a script which varies the amount of damping present in the system (in rear and front suspensions) beginning with an initial value and constrained between an upper-bound value and a lower-bound value. A global optimization algorithm is applied for calculating the damping value for which the aggregate motion of the system is minimal. The aggregate motion of the system is represented by a cost function (`cost_function.m`) in which both bounce and pitch amplitudes have a contribution.

Apart from the MATLAB scripts, a simulink model `two_dof_system.slx` also has been included which is modelled on the basis of the equations of motions describing the system.

