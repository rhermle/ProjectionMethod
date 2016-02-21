function [ P U V X Y XP YP numXCells numYCells dx] = Stokes2DPC(g, numYCells, p0, mu, height, width, R, L, numTimeSteps, debug)

numXCells = numYCells;
%numXCells = round((width / dx) + 1);
numXPCells = numXCells + 1;
numYPCells = numYCells + 1;
y = linspace(-height / 2, height / 2, numYCells);
x = linspace(0, width, numXCells);
dx = y(2) - y(1);
yp = linspace(-height / 2 - dx/2, height / 2 + dx/2, numYPCells);
xp = linspace(0 - dx/2, width + dx/2, numXPCells);
dt = .1*dx.^2;
[X Y] = meshgrid(x,y);
[XP YP] = meshgrid(xp,yp);

%ind will act like a sub2ind call, transforming 2D subscripts into a
%linear index
ind = zeros(numXCells,numYCells);
count = 1;
for i=1:numXCells
    for j=1:numYCells
        ind(i,j) = count;
        count = count + 1;
    end
end

%indp will act like a sub2ind call, transforming 2D subscripts into a
%linear index
pind = zeros(numXPCells,numYPCells);
count = 1;
for i=1:numXPCells
    for j=1:numYPCells
        pind(i,j) = count;
        count = count + 1;
    end
end

[F1 F2] = force(X, Y, mu, width, p0, g, R, dx, L);

if debug
    figure();
    surf(X,Y,F1);
    title('F1');
    figure();
    surf(X,Y,F2);
    title('F2');
end

P = zeros(numYPCells, numXPCells, numTimeSteps);
U = zeros(numYCells, numXCells, numTimeSteps);
V = zeros(numYCells, numXCells, numTimeSteps);

u_star = zeros(numXCells*numYCells, numTimeSteps);
v_star = zeros(numXCells*numYCells, numTimeSteps);
u = zeros(numXCells*numYCells, numTimeSteps);
v = zeros(numXCells*numYCells, numTimeSteps);
p = zeros(numXPCells*numYPCells, numTimeSteps);

for t = 1:numTimeSteps  %This loop will surround the entire algoritm

%Set initial condition for u, flow of 0
if(t == 1)
    uinit = zeros(numXCells*numYCells);
    for i=1:numXCells
        for j=1:numYCells
            uinit(ind(i,j), t) = 0;
            vinit(ind(i,j), t) = 0;
        end
    end
end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Step 1 of algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
    %Question:  Assume small strain tensor is symmetric, so transpose goes away?
    %Random Question:  Do they ever use FFT to solve these?
    
    %u_star = u + mu*(uxx + uyy) + F1
    %v_star = v + mu*(vxx + vyy) + F2
    
    if t == 1
    
        for i=2:numXCells-1
            for j=2:numYCells-1
                u_star(ind(i,j),t) = uinit(ind(i,j)) + dt * ((1/(dx^2)) * mu * (uinit(ind(i-1,j)) - 4 * uinit(ind(i,j)) + uinit(ind(i+1,j)) + uinit(ind(i,j+1)) + uinit(ind(i,j-1))) + F1(j,i));
                v_star(ind(i,j),t) = vinit(ind(i,j)) + dt * ((1/(dx^2)) * mu * (vinit(ind(i-1,j)) - 4 * vinit(ind(i,j)) + vinit(ind(i+1,j)) + vinit(ind(i,j+1)) + vinit(ind(i,j-1))) + F2(j,i));
            end
        end
    
    else
        
        for i=2:numXCells-1
            for j=2:numYCells-1
                u_star(ind(i,j),t) = u(ind(i,j),t-1) + dt * ((1/(dx^2)) * mu * (u(ind(i-1,j),t-1) - 4 * u(ind(i,j),t-1) + u(ind(i+1,j),t-1) + u(ind(i,j+1),t-1) + u(ind(i,j-1),t-1)) + F1(j,i));
                v_star(ind(i,j),t) = v(ind(i,j),t-1) + dt * ((1/(dx^2)) * mu * (v(ind(i-1,j),t-1) - 4 * v(ind(i,j),t-1) + v(ind(i+1,j),t-1) + v(ind(i,j+1),t-1) + v(ind(i,j-1),t-1)) + F2(j,i));
            end
        end
        
    end
    
    %Set the boundary condition for u_star to be the same as u
    
    %Left edge, excluding corners, use forward difference for uxx and vxx
    %Right edge, excluding corners, use backward difference for uxx and vxx
    for j = 1:numYCells
        i = 1;
        u_star(ind(i,j),t) = 0;
        v_star(ind(i,j),t) = 0;
        
        i = numXCells;
        u_star(ind(i,j), t) = 0;
        v_star(ind(i,j), t) = 0;
    end
    
    %Bottom edge, excluding corners, use forward difference for uyy and vyy
    %Top edge, excluding corners, use backward difference for uyy and vyy
    for i = 2:numXCells-1
        j = 1;
        %Want uy = 0
        %u_star(ind(i,j),t) = u_star(ind(i,j+1),t);
        u_star(ind(i,j),t) = 0;
        v_star(ind(i,j),t) = 0;
        
        j = numYCells;
        %u_star(ind(i,j),t) = u_star(ind(i,j-1),t);
        u_star(ind(i,j),t) = 0;
        v_star(ind(i,j),t) = 0;
    end
    
    tmp1 = u_star(:, t);
    tmp2 = v_star(:, t);

    %We must transpose here because in matlab rows will map to Y values and
    %columns will map to X values, which is the opposite of how we have
    %defined it in our code.
       
    U_Star_Temp = tmp1(ind)';
    V_Star_Temp = tmp2(ind)';
   
    if debug
        
        figure(56789)
        quiver(X,Y,U_Star_Temp,V_Star_Temp);
        title('U* Quiver');
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Step 2 of algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %p_xx + p_yy = u_star_x + v_star_y
    %p(i-1,j) - 4p_(i,j) + p(i+1,j) + p(i,j+1) + p(i,j-1) = (d/2) * (u_star(i+1,j) - u_star(i-1,j) + v_star(i,j+1) - v_star(i,j-1))  

    
%%%%%%%%%%%%%  This will calculate the analytic pressure for debugging.
%     z = sqrt( (X-(L+R)).^2 + Y.^2 ) - R;
%     E = .5 * R;
%     id = -E < z & z < E;
%     id2 = z <= -E;
% 
%     temp = zeros(numYCells,numXCells);
% 
%     temp(id) = -1/(2*R)*(1 - z(id)/E - sin(pi*z(id)/E)/pi);
%     temp(id2) = -1/(R);
%     
%     for i = 1:numXCells
%         for j = 1:numYCells
%             p(ind(i,j),t) = temp(j,i);
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%This will calculate the analytic gradient for debugging
%     z = sqrt( (X-(L+R)).^2 + Y.^2 ) - R;
%     E = .5 * R;
%     
%     delta = zeros(numYCells, numXCells);
%     id = -E <= z & z <= E;
%     delta(id) = (1 + cos(pi*(-z(id)/E)))/(2*E);
%     delta_ = zeros(numYCells, numXCells);
%     delta_(id) = -pi*sin(pi*(-z(id)/E))/(2*E^2);
%     lapP = -1/R * delta_ + (1/R * delta) ./ (z+R);
%     
%     lapPHat = zeros(numXCells * numYCells);
%     
%     for i = 1:numXCells
%         for j = 1:numYCells
%             lapPHat(ind(i,j)) = lapP(j,i);
%         end
%     end
    

%     id = -E < z & z < E;
%     id2 = z <= -E;

%     px = zeros(numYCells,numXCells);
%     py = zeros(numYCells,numXCells);
% 
%     px(id) = (1/(2*R*E)) * ((X(id) - (L+R)) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 )) + (cos(pi*z(id)/E)) .*(X(id) - (L+R)) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 + eps)));
%     px(id2) = 0;
% 
%     py(id) = (1/(2*R*E)) * (Y(id) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 )) + (cos(pi*z(id)/E)) .*Y(id) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 + eps)));
%     py(id2) = 0;
% 
%     pxHat = zeros(numXCells*numYCells,1);
%     pyHat = zeros(numXCells*numYCells,1);
% 
%     for i = 1:numXCells
%         for j = 1:numYCells
%             pxHat(ind(i,j)) = px(j,i);
%             pyHat(ind(i,j)) = py(j,i);
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = zeros(numXPCells*numYPCells,numXPCells*numYPCells);
    b = zeros(numXPCells*numYPCells,1);

%     %interior points
    for i=2:numXPCells-1
        for j=2:numYPCells-1

            A(pind(i,j), pind(i-1,j)) = 1;
            A(pind(i,j), pind(i,j)) = -4;
            A(pind(i,j), pind(i+1,j)) = 1;
            A(pind(i,j), pind(i,j+1)) = 1;
            A(pind(i,j), pind(i,j-1)) = 1;
            
            b(pind(i,j)) = (1/dt)*(dx/2) * (u_star(ind(i,j),t) - u_star(ind(i-1,j),t) + ...
                                    u_star(ind(i,j-1),t) - u_star(ind(i-1,j-1),t) + ... 
                                    v_star(ind(i,j),t ) - v_star(ind(i,j-1),t) + ...
                                    v_star(ind(i-1,j),t ) - v_star(ind(i-1,j-1),t));

%             b(ind(i,j)) = (1/2) * ((pxHat(ind(i+1,j+1)) - pxHat(ind(i,j+1)))/dx + ...
%                                             (pxHat(ind(i+1,j)) - pxHat(ind(i,j)))/dx + ... 
%                                             (pyHat(ind(i+1,j+1)) - pyHat(ind(i+1,j)))/dx + ...
%                                             (pyHat(ind(i,j+1)) - pyHat(ind(i,j)))/dx);


%         b(ind(i,j)) = lapPHat(ind(i,j));
        
        
        end
    end

    %boundary conditions
    %Left and right edges

    for j = 2:numYPCells-1
        
        A(pind(1,j), pind(1,j)) = 1;
        b(pind(1,j)) = 0;
        A(pind(numXPCells,j), pind(numXPCells,j)) = 1;
        b(pind(numXPCells,j)) = 0;
    end

    %Top and bottom edges
    for i = 1:numXPCells

        A(pind(i,1),pind(i,1)) = 1;
        b(pind(i,1)) = 0;
        A(pind(i,numYPCells),pind(i,numYPCells)) = 1;
        b(pind(i,numYPCells)) = 0;
        
    end
 
    %p(:, t) = pinv(A) * b;
    
    if debug
        
        figure();
        surf(X,Y,b(ind(:,:)));
        message = sprintf('b, t = %f',t);
        title(message);
        
    end
    
    A = sparse(A);
    
    p(:,t) = A\b;
    %rcond(A)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%Step 3 of algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
    %u_next = u_star - grad(p) 
      
        
%%%%%%%%%%%%% Analytic p gradient, for debugging
%         z = sqrt( (X-(L+R)).^2 + Y.^2 ) - R;
%         E = .5 * R;
%         id = -E < z & z < E;
%         id2 = z <= -E;
%   
%         px = zeros(numYCells,numXCells);
%         py = zeros(numYCells,numXCells);
% 
%         px(id) = (1/(2*R*E)) * ((X(id) - (L+R)) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 )) + (cos(pi*z(id)/E)) .*(X(id) - (L+R)) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 + eps)));
%         px(id2) = 0;
% 
%         py(id) = (1/(2*R*E)) * (Y(id) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 )) + (cos(pi*z(id)/E)) .*Y(id) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 + eps)));
%         py(id2) = 0;
%         
%         pxHat = zeros(numXCells*numYCells,1);
%         pyHat = zeros(numXCells*numYCells,1);
%         
%         for i = 1:numXCells
%             for j = 1:numYCells
%                 pxHat(ind(i,j)) = px(j,i);
%                 pyHat(ind(i,j)) = py(j,i);
%             end
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    for i=1:numXCells
        for j=1:numYCells 

            if j==1 || j==numYCells || i==1 || i==numXCells
                v(ind(i,j),t) = 0; 
                u(ind(i,j),t) = 0; 
            else
                 u(ind(i,j),t) = u_star(ind(i,j),t) - dt * (1/(2*dx)) * (p(pind(i+1,j+1),t) - p(pind(i,j+1),t) + p(pind(i+1,j),t) - p(pind(i,j),t));
                 v(ind(i,j),t) = v_star(ind(i,j),t) - dt * (1/(2*dx)) * (p(pind(i+1,j+1),t) - p(pind(i+1,j),t) + p(pind(i,j+1),t) - p(pind(i,j),t));

%                      u(ind(i,j),t+1) = u_star(ind(i,j),t) - dt * pxHat(ind(i,j));
%                      v(ind(i,j),t+1) = v_star(ind(i,j),t) - dt * pyHat(ind(i,j));

            end
        end    
    end        
    
    tmp1 = p(:, t);
    tmp2 = u(:, t);
    tmp3 = v(:, t);

    %We must transpose here because in matlab rows will map to Y values and
    %columns will map to X values, which is the opposite of how we have
    %defined it in our code.
    
    %P(:,:,t) = temp;
    P(:,:,t) = tmp1(pind)';
    U(:,:,t) = tmp2(ind)';
    V(:,:,t) = tmp3(ind)';
    
end
    
    
    
    