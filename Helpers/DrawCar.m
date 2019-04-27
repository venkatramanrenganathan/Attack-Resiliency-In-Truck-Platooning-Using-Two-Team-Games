%
% This function draws a unicycle robot.
%
% Inputs:
%          param:  1x3 array
%           - state = [a, b, delta, theta]
%           - (a,b): the center of the robot
%           - delta: the steering angle
%           - theta: the heading angle
%           - L    : wheelbase
%
% Output:
%          - hdl: plot handles
%
%
% (C) Kaveh Fathian, 2018.  Email: kaveh.fathian@gmail.com
%
function hdl = DrawCar(state, L, hdl)

if nargin == 2
    hdl = {[], [], [], [], []};
end

%%
x = state(1);
y = state(2);
delta = state(3); % Steering angle in a global coordinate frame
theta = state(4); % Heading angle in a global coordinate frame



%% Prameters
bodyLength = 1.66*L; 
bodyWidth = 0.62*bodyLength;

wheelLength = 0.5*bodyWidth;
wheelWidth = 0.15*bodyWidth;

bodyShift = 0.3*bodyLength;
xBody = x - bodyShift * cos(theta);
yBody = y - bodyShift * sin(theta);

bodyStyle.FaceColor = [255 99 0] ./ 255; % Orange
bodyStyle.EdgeColor = [244 164 96] ./ 255; % Light orange
bodyStyle.LineWidth = 10*bodyWidth;

wheelShift = (bodyWidth + wheelWidth) ./ 2; 
xWheel1 = x - wheelShift * sin(theta);
yWheel1 = y + wheelShift * cos(theta);
xWheel2 = x + wheelShift * sin(theta);
yWheel2 = y - wheelShift * cos(theta);

xWheel3 = x - L * cos(theta) - wheelShift * sin(theta);
yWheel3 = y - L * sin(theta) + wheelShift * cos(theta);
xWheel4 = x - L * cos(theta) + wheelShift * sin(theta);
yWheel4 = y - L * sin(theta) - wheelShift * cos(theta);


wheelStyle.FaceColor = [139 69 19] ./ 255; % Orange
wheelStyle.EdgeColor = [139 69 19] ./ 255; % Brown
wheelStyle.LineWidth = 3*bodyWidth;

%% Draw the robot

% Main body
param = [xBody,yBody, bodyLength, bodyWidth, theta];
hb = DrawRectangle(param, bodyStyle, hdl{1});

% Front Wheels
param = [xWheel1, yWheel1, wheelLength, wheelWidth, delta];
hw1 = DrawRectangle(param, wheelStyle, hdl{2});
param = [xWheel2, yWheel2, wheelLength, wheelWidth, delta];
hw2 = DrawRectangle(param, wheelStyle, hdl{3});


% Back Wheels
param = [xWheel3, yWheel3, wheelLength, wheelWidth, theta];
hw3 = DrawRectangle(param, wheelStyle, hdl{4});
param = [xWheel4, yWheel4, wheelLength, wheelWidth, theta];
hw4 = DrawRectangle(param, wheelStyle, hdl{5});

hdl = {hb, hw1, hw2, hw3, hw4};

end