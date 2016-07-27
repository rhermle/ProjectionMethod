function [ f1 f2 ] = force(x, y, mu, width, p0, g, R, dx, L)

graph = 0;

%E = dx;
E = .5 * R;
z = sqrt((x - (R + L)).^2 + y.^2) - R;

id1 = (z >= -E) & (z <= E);

delta = zeros(size(x));
delta (id1) = (1 + cos((pi * z(id1) / E))) / (2 * E);

f1 = (delta/R) .* ((x - (R + L))) ./ sqrt((x - (R + L)).^2 + y.^2 + eps);
f2 = (delta/R) .* (y ./ (sqrt((x - (R + L)).^2 + y.^2 + eps)));

if graph

    figure(12345);
    hold on;
    quiver(x,y,f1,f2);
    %title('F1 and F2');
    %quiver(x(1:quivRes:end,1:quivRes:end),Y(1:quivRes:end,1:quivRes:end),U(1:quivRes:end,1:quivRes:end, i),V(1:quivRes:end,1:quivRes:end, i));
    contour(x,y,z,[0,0]);
    hold off;
    
end


end