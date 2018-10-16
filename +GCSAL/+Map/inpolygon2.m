function [in] = inpolygon2(x,y,xv,yv)
% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
%
% Function to determine if points (x,y) is inside or outside a polygon.
%
% INPUTS
%   (xv,yv): polygon is specified by (xv,yv) points
%   (x,y): point coordinates are specified as (x,y) pairs. x and y should be vectors of same size.
% OUTPUTS:
%   in : is a logical array (0 means point is outside and 1 means point is inside)
%
% Implementation is based on winding algorithm explained in http://geomalgorithms.com/a03-_inclusion.html
% Example
%       xv = rand(6,1); yv = rand(6,1);
%       xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
%       x = rand(1000,1); y = rand(1000,1);
%       in = inpolygon(x,y,xv,yv);
%       plot(xv,yv,x(in),y(in),'.r',x(~in),y(~in),'.b')


if ((xv(1) ~= xv(end)) || (yv(1) ~= yv(end)))
        xv = [xv ; xv(1)];
        yv = [yv ; yv(1)];
end

n=length(xv)-1; % number of polygon corners
np=length(x); % number of points to be evaluated
wn=zeros(np,1); % starts with all points outside

for j=1:np
    for i=1:n
        if (yv(i) <=y(j))
            if (yv(i+1) > y(j))
                if (is_point_on_left(x(j),y(j),xv(i),yv(i),xv(i+1),yv(i+1))>0)
                    wn(j)=wn(j)+1;
                end
            end
        else
            if (yv(i+1) <= y(j))
                if (is_point_on_left(x(j),y(j),xv(i),yv(i),xv(i+1),yv(i+1))<0)
                    wn(j)=wn(j)-1;
                end
            end
        end
    end
end

in=logical(wn); % convert to logical 0-1


function isleft=is_point_on_left(px,py,p1Lx,p1Ly,p2Lx,p2Ly)
% Determine if point (px,py) is on the left | right | On the line.
% points on the line is specified as (p1Lx,p1Ly) & (p2Lx,p2Ly)
% isleft=1 if p is on the left , isleft= 0 if p is along the line
% isleft=-1 if p is on the right.

p1_to_p2=[p2Lx-p1Lx;p2Ly-p1Ly];
p1_to_p=[px-p1Lx;py-p1Ly];
isleft = det([p1_to_p2 p1_to_p]);
end

end
