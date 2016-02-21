function [P U V]  = pTest(X,Y,R,L)

% y = linspace(-h / 2, h / 2, ny);
% dx = y(2) - y(1);
% nx = round((w / dx) + 1);
% %X = linspace(-width / 2, width / 2, numXCells);
% x = linspace(0, w, nx);
% %x = y + h/2;
% [X Y] = meshgrid(x,y);

dx = Y(2,1) - Y(1,1);

P = zeros(size(X));
U = zeros(size(X));
V = zeros(size(X));

%%%%%%%%%%%%%%The value of E might change%%%%%%
E = .5 * R;
%E = dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = sqrt( (X-(L+R)).^2 + Y.^2 ) - R;

id = -E <= z & z <= E;
id2 = z < -E;

P(id) = -1/(2*R)*(1 - z(id)/E - sin(pi*z(id)/E)/pi);
P(id2) = -1/(R);

%temp = zeros(ny,nx);

% for i = 1:ts
%     temp(id) = -1/(2*R)*(1 - z(id)/E - sin(pi*z(id)/E)/pi);
%     temp(id2) = -1/(R);
%     P(:,:,i) = temp;
% end

end