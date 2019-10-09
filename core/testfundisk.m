function f = testfundisk(r,theta,kk,K)

%   INPUT    
%   r,theta     polar coordinates to evaluate the functions
%   kk          index for one of the test functions
%   K           parameter for the test function
%
%   OUTPUT     
%   f           resulting function values

% In case we need cartesian coordinates
x = r.*cos(theta); y = r.*sin(theta);

switch kk
               
    case 1  % Basis function I (cosine)
     
    f = cos(K(1)*acos(r)).*cos(K(2)*theta);

        
    case 2  % Basis function II (sine)
    
    f = cos(K(1)*acos(r)).*sin(K(2)*theta);
    
    
    case 3  % Two Gaussians
        
    f = 0.75*exp(-5*((x-0.6).^2 + (y-0.12).^2))+ 0.65*exp(-6*((x+0.5).^2 + (y-0.3).^2));
        
    case 4  % The cross (test function in form of a cross)
    
    f = max(exp(-40*x.^2),exp(-40*y.^2));
    
    case 5  % Rosenbrock function (very simple polynomial test function)
    
    a = 1;
    b = 100;
    lambda = 1;
        
    f = (a-lambda*x).^2 + b*(lambda*y-lambda^2*x.^2).^2;
    
    case 6  % Characteristic function of L (to show some Gibbs artifacts)
        
    f = (abs(x)<= 0.5).*(abs(y)<= 0.5)-(abs(x-0.25)<= 0.25).*(abs(y-0.25)<= 0.25);
    
    case 7  % A non-smooth function
        
    a = 2;
    b = 2.25;
        
    f = exp(-0.4*sqrt((8*x-a).^2+(9*y-b).^2)).*...
            cos(1.2*sqrt((8*x-a).^2+(9*y-b).^2));
        
    case 8  % The eye
        
    a = 0.5;
    b = 1;
        
    f = exp(-0.08*((8*x-a).^2+(12*y-b).^2)).*...
            cos(0.25*((8*x-a).^2+(12*y-b).^2));
        
    case 9  % The runge function
        
    a = 10;
    
    f = 1./(1+a*(x.^2+y.^2));
    
    case 10  % Tensor-product of two oscillating functions
        
    f = cos(17*x).*cos(13*y);
    
    otherwise
        
    error('There is no function associated to this number');

end