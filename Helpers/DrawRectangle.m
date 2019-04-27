%
% This function draws a rectangle.
%
% Inputs:
%          param:  1x5 array
%           - param = [a, b, w, h, theta]
%           - (a,b): the center of the rectangle
%           - (w,h): width and height of the rectangle > 0
%           - theta: the rotation angle of the rectangle
%          style:  plot style struct
%           - style.FaceColor, style.EdgeColor, style.LineWidth
%          hdl: previous graphic handel
%
% Output:
%          - hdl: plot handle
%
% Code based on: 'DrawRectangle.m' written by Rasoul Mojtahedzadeh 
% (mojtahedzadeh@gmail.com), November, 2011.
%
% (C) Kaveh Fathian, 2018.  Email: kaveh.fathian@gmail.com
%
function hdl = DrawRectangle(param, style, hdl)

a = param(1);
b = param(2);
w = param(3);
h = param(4);
theta = param(5);

X = [-w/2 w/2 w/2 -w/2 -w/2];
Y = [h/2 h/2 -h/2 -h/2 h/2];
P = [X; Y];
ct = cos(theta);
st = sin(theta);
R = [ct, -st; st, ct];
P = R * P;

FaceColor = style.FaceColor;
EdgeColor = style.EdgeColor;
LineWidth = style.LineWidth;


%%

if isempty(hdl)
    hdl = patch(P(1,:)+a, P(2,:)+b, FaceColor);
    set(hdl, 'EdgeColor',  EdgeColor);
    set(hdl, 'LineWidth',  LineWidth);
else
    set(hdl,'XData', P(1,:)+a);
    set(hdl,'YData', P(2,:)+b);
end




end