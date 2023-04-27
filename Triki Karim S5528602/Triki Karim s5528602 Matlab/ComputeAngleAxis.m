function R = ComputeAngleAxis(theta,v)
%Implement here the Rodrigues formula
I = eye(3);
q=[0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R= I+ q*sin(theta)+ (1-cos(theta))*(q*q);
end
