function [rot_matrix] = quatToRot(q0,q1,q2,q3)
% quatToRot convert a quaternion into a rotation matrix
    %Covert a quaternion into a full three-dimensional rotation matrix.
 
    %Input
    %:param Q: A 4 element array representing the quaternion (q0,q1,q2,q3) 
    q =[q0,q1,q2,q3];
    %Output
    %return: A 3x3 element matrix representing the full 3D rotation matrix. 
    %First row of the rotation matrix
    %...
    R00 = 2 * (q0 * q0 + q1 * q1) - 1;
    R01 = 2 * (q1 * q2 - q0 * q3);
    R02 = 2 * (q1 * q3 + q0 * q2);
    %Second row of the rotation matrix
    %...
    R10 = 2 * (q1 * q2 + q0 * q3);
    R11 = 2 * (q0 * q0 + q2 * q2) - 1;
    R12 = 2 * (q2 * q3 - q0 * q1);
     
    %Third row of the rotation matrix
    %...
    R20 = 2 * (q1 * q3 - q0 * q2);
    R21 = 2 * (q2 * q3 + q0 * q1);
    R22 = 2 * (q0 * q0 + q3 * q3) - 1;
    %3x3 rotation matrix
    [rot_matrix] =[R00 R01 R02 ; R10 R11 R12; R20 R21 R22];

end