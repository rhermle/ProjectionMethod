function [ f1 f2 ] = force(x, y, mu, width, p0, g, R, dx, L)

%alpha = 1;

theta = atan2(y,x - (R + L));
%alpha = 1 * (cos(2*theta));
alpha = 20 * (-cos(theta));

id2 = theta > (-pi / 2) & theta < (pi / 2);

alpha(id2) = 0;

%E = 2 * dx;
%E = 2 * .125;
E = .25 * R;

%z = x.^2 + y.^2 - R^2;
z = sqrt((x - (R + L)).^2 + y.^2) - R;

id1 = (z > -E) & (z < E);

delta = zeros(size(x));
delta (id1) = (1 + cos((pi * z(id1) / E))) / (2 * E);
%2E in denominator?

f1 = delta .* cos(theta);
f2 = delta .* sin(theta);

% f1 = alpha .* delta .* (2*(x - (R+L)));
% f2 = alpha .* delta .* (2*y);

f1 = f1';
f2 = f2';

end