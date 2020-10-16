# 2-DOF-free-vibrations
Code for calculating vibration characteristics of 2 degree of freedom mechanical systems.
The 2 DOF system is assumed to be a simple car model with its mass concentrated in a rectangular mass which can translate vertically (bounce) and rotate about its central axis normal to the plane (pitch). The front and rear suspension of the car is modelled using a spring and damper system connected to the rectangle at a particular distance from the centre of mass.
There is no source of excitation for the model and the car has some initial bounce and pitch motion, thus the sytem will experience free vibrations which can be mathematically modelled by a 2nd order coupled (or uncoupled depending on the system parameters) homogeneous ordinary differential equation.
3 approaches have been taken in the code to solve for the system:
1). Analytical approach where the code is built according to the same exact steps which are taken while analytically solving the equation of motion (generating the characteristic polynomial and then calculating the eigen values from the polynomial). This approach has the fastest computing time, but is limited as it can be used only for simpler systems having shallow modelling depth.
2). Numerical approach where the equation of motion is represented in the state space vector form (derivative of state_matrix = system_matrix*(state_matrix)) and the ode45 solver in MATLAB is used with the objective function whose independent variables are the time vector and the state vector while the dependent vector is the first derivate of the state vector. This approach is a little slower than the analytical approach but can be used with ease in systems having complex modelling depths.
3). Symbolic approach where the equations of motion (2nd order homogeneous ODE) is solved symbollically. This approach is much slower than the other two, and while can be used in complex problems, the excess time taken will not be worth it.

Apart from just solving the system, fast fourier transform is also applied to the system in order to convert from the time domain (x vs t) to the frequency domain (x vs f) which is much more useful as we are able to observe the natural frequencies of the system and thus design the suspension of the car keeping in mind the natural frequencies which must be avoided.


