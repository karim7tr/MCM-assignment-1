%% Exercises Modelling Part 1
% Rotation matrices, Equivalent angle-axis representations, Quaternions
addpath('include') %%DO NOT CHANGE STUFF INSIDE THIS PATH

%% Exercise 1
% Write a function ca
% lled ComputeAngleAxis() implementing the Rodrigues Formula, 
% taking in input the geometric unit vector v and the rotation angle theta
% and returning the orientation matrix.
% and test it for the following cases:

% 1.1.
%we had create the rotation matrix by giving the angle & the vector v

% vector v : The Geometric Unit Vector

% 1.2.
v=[1 0 0];
theta= pi/6; % 30°= pi/6 rad

aRb = ComputeAngleAxis(theta, v);
disp('aRb ex 1.2:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.3:');disp(theta);
disp('v ex 1.3:');disp(v);

% 1.3.
v=[0 1 0];
theta= pi/4;

aRb = ComputeAngleAxis(theta, v);
disp('aRb ex 1.3:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.4:');disp(theta);
disp('v ex 1.4:');disp(v);
% 1.4.
v=[0 0 1];
theta= pi/2;

aRb = ComputeAngleAxis(theta, v);

disp('aRb ex 1.4:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.5:');disp(theta);
disp('v ex 1.5:');disp(v);

% 1.5.

%from question 1.6 to 1.8 we have the rotation vector so we need to divide it
%to angel axis & vector v
% 1.6.

 p = [0, pi/2, 0];
 theta = norm(p);
 v = p/theta;

aRb = ComputeAngleAxis(theta, v);
disp('aRb ex 1.6:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.7:');disp(theta);
disp('v ex 1.7:');disp(v);

% 1.7.
 p = [0.4, -0.3, -0.3];

 theta = norm(p);
 v = p/theta;

aRb = ComputeAngleAxis(theta, v);
disp('aRb ex 1.7:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.8:');disp(theta);
disp('v ex 1.8:');disp(v);

% 1.8.
 p = [-pi/4, -pi/3, pi/8];

 theta = norm(p);
 v = p/theta;

aRb = ComputeAngleAxis(theta, v);
disp('aRb ex 1.8:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.8:');disp(theta);
disp('v ex 1.8:');disp(v);
%% Exercise 2
% 2.1. Write the relative rotation matrix aRb

% 2.2. Solve the Inverse Equivalent Angle-Axis Problem for the orientation matrix aRb. 

% 2.3. Repeat the exercises using the wRc instead of wRa (more general example)
% NB: check the notation used !

    % 2.1 
    % Initialize the rotation matrices, using the suggested notation
        %rotation matrix from <w> to frame <a> (scriviamo alla lavagna)
        %rotation matrix from <w> to <b> (represent 90° around z)
    % Compute the rotation matrix between frame <a> and <b>

    v=[0 0 1]; %rotation arround z-axis
    theta= pi/2;
    wRb=ComputeAngleAxis(theta,v);
    aRw=eye(3);% <a> & <w>  are identical and parallel
    aRb =aRw*wRb;
    disp('wRb ex 2.1:');disp(wRb);
    disp('aRb ex 2.1:');disp(aRb);
    % 2.2
    % Compute the inverse equivalent angle-axis repr. of aRb s
    [theta, v]=ComputeInverseAngleAxis(aRb);
    % Plot Results
    disp('aRb ex 2.1:');disp(aRb);
    plotRotation(theta,v,aRb);
    disp('theta ex 2.2:');disp(theta);
    disp('v ex 2.2:');disp(v); 

    % 2.3
    wRc = [0.835959, -0.283542, -0.469869; 
        0.271321, 0.957764, -0.0952472;
        0.47703, -0.0478627, 0.877583];

    % Compute the rotation matrix between frame <c> and <b>
     cRb = (wRc')*wRb ;
    % Compute inverse equivalent angle-axis repr. of cRb
    [theta, v] = ComputeInverseAngleAxis(cRb);
    % Plot Results
    plotRotation(theta,v,cRb);
    disp('theta ex 2.3:');disp(theta);
    disp('v ex 2.3:');disp(v); 

%% Exercise 3
% 3.1 Given two generic frames < w > and < b >, define the elementary 
% orientation matrices for frame < b > with respect to frame < w >, knowing
% that:
    % b. < b > is rotated of 45◦ around the y-axis of < w >
    % c. < b > is rotated of 15◦ around the x-axis of < w >
% 3.2 Compute the equivalent angle-axis representation for each elementary rotation
% 3.3 Compute the z-y-x (yaw,pitch,roll) representation using the already
% computed matrices and solve the Inverse Equivalent Angle-Axis Problem
% 3.4 Compute the z-x-z representation using the already computed matrices 
% and solve the Inverse Equivalent Angle-Axis Problem 

% 3.1
% hint: define angle of rotation in the initialization

    % a
        %rotation matrix from <w> to frame <b> by rotating around z-axes
        wRb_z = [cos(pi/6) -sin(pi/6) 0 ; sin(pi/6) cos(pi/6) 0; 0 0 1];

    % b
        %rotation matrix from <w> to frame <b> by rotating around y-axes
         wRb_y = [cos(pi/4) 0 sin(pi/4) ; 0 1 0 ; -sin(pi/4) 0 cos(pi/4)];

    % c
        %rotation matrix from <w> to frame <b> by rotating around x-axes
        wRb_x = [1 0 0 ; 0 cos(pi/12) -sin(pi/12) ; 0 sin(pi/12) cos(pi/12)];


    disp('es 3.1:');disp(wRb_z);disp(wRb_y);disp(wRb_x);
    
% 3.2
    % a
        [theta, v] = ComputeInverseAngleAxis(wRb_z);
        % Plot Results
        plotRotation(theta,v,wRb_z);
        disp('theta ex 3.2.a:');disp(theta);
        disp('v ex 3.2.a:');disp(v); 
    % b
        [theta, v] = ComputeInverseAngleAxis(wRb_y);
        % Plot Results
        plotRotation(theta,v,wRb_y);
        disp('theta ex 3.2.b:');disp(theta);
        disp('v ex 3.2.b:');disp(v); 
    % c
        [theta, v] = ComputeInverseAngleAxis(wRb_x);
        % Plot Results
        plotRotation(theta,v,wRb_x);
        disp('theta ex 3.2.c:');disp(theta);
        disp('v ex 3.2.c:');disp(v); 
% 3.3 
    % Compute the rotation matrix corresponding to the z-y-x representation;
    Rzyx= wRb_z*wRb_y*wRb_x;
    disp('Rzyx 3.3:');disp(Rzyx);
    % Compute equivalent angle-axis repr.
    [theta, v] = ComputeInverseAngleAxis(Rzyx);
    % Plot Results
    plotRotation(theta,v,Rzyx);
    disp('theta ex 3.3:');disp(theta);
    disp('v ex 3.3:');disp(v);
% 3.4
    % Compute the rotation matrix corresponding to the z-x-z representation;
    Rzxz= wRb_z*wRb_x*wRb_z;
    disp('Rzxz 3.4:');disp(Rzxz);
	% Compute equivalent angle-axis repr.
    [theta, v] = ComputeInverseAngleAxis(Rzxz);
    % Plot Results
    plotRotation(theta,v,Rzxz);
    disp('theta ex 3.4:');disp(theta);
    disp('v ex 3.4:');disp(v);
   
%% Exercise 4
% 4.1 Represent the following quaternion with the equivalent angle-axis
 %q = 0.8924 +  0.23912i +  0.36964j + 0.099046k;

% 4.2 Repeat the exercise using the built-in matlab functions see:
% https://it.mathworks.com/help/nav/referencelist.html?type=function&category=coordinate-system-transformations&s_tid=CRUX_topnav
% CHECK IF THE RESULT IS THE SAME 
    %%%%%%%%%%%%
    a = 0.8924; % real part % check notazione usata in classe
    b = 0.23912;
    c = 0.36964;
    d = 0.099046;
    q=[a b c d];
    %%%%%%%%%%%%%
    % Compute the rotation matrix associated with the given quaternion
    % using quatToRot fonction
    rotMatrix = quatToRot(a,b,c,d);
    disp('rot matrix es 4.1');disp(rotMatrix)
    % solve using matlab functions quat2rotm
     rotMatrix2= quat2rotm(q);

    % Evaluate angle-axis representation and display rotations - check if the
    % same results as before
    [theta, v] = ComputeInverseAngleAxis(rotMatrix);
    % Plot Results
    plotRotation(theta,v,rotMatrix);
    disp('rot matrix2 es 4.2');disp(rotMatrix2)