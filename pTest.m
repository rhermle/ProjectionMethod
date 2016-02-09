function [P U V X Y nx ny]  = pTest(h, w, R, L, ts, ny)

%%%%%%%%%%%%%%The value of E might change%%%%%%
E = .5 * R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = linspace(-h / 2, h / 2, ny);
dx = y(2) - y(1);
nx = round((w / dx) + 1);
%X = linspace(-width / 2, width / 2, numXCells);
x = linspace(0, w, nx);
[X Y] = meshgrid(x,y);

P = zeros(ny, nx, ts);
U = zeros(ny, nx, ts);
V = zeros(ny, nx, ts);

z = sqrt( (X-(L+R)).^2 + Y.^2 ) - R;

id = -E < z & z < E;
id2 = z < -E;

temp = zeros(ny,nx);

for i = 1:ts
    temp(id) = -1/(2*R)*(1 - z(id)/E - sin(pi*z(id)/E)/pi);
    temp(id2) = -1/(R);
    P(:,:,i) = temp;
end

end