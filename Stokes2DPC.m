function [ P U V X Y numXCells numYCells dx] = Stokes2DPC(g, numYCells, p0, mu, height, width, R, L, numTimeSteps, debug)

y = linspace(-height / 2, height / 2, numYCells);
dx = y(2) - y(1);
numXCells = round((width / dx) + 1);
%X = linspace(-width / 2, width / 2, numXCells);
x = linspace(0, width, numXCells);
[X Y] = meshgrid(x,y);

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

[F1 F2] = force(X, Y, mu, width, p0, g, R, dx, L);

if debug
    figure();
    surf(X,Y,F1);
    title('F1');
    figure();
    surf(X,Y,F2);
    title('F2');
end

P = zeros(numYCells, numXCells, numTimeSteps);
U = zeros(numYCells, numXCells, numTimeSteps);
V = zeros(numYCells, numXCells, numTimeSteps);

u_star = zeros(numXCells*numYCells, numTimeSteps);
v_star = zeros(numXCells*numYCells, numTimeSteps);
u = zeros(numXCells*numYCells, numTimeSteps);
v = zeros(numXCells*numYCells, numTimeSteps);
p = zeros(numXCells*numYCells, numTimeSteps);

for t = 1:numTimeSteps  %This loop will surround the entire algoritm

    %Set initial condition for u, flow of 1 to the right everywhere
    if(t == 1)
        for i=1:numXCells
            for j=1:numYCells
                u(ind(i,j), t) = 0;
                v(ind(i,j), t) = 0;
            end
        end
    end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Step 1 of algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
    %Question:  Assume small strain tensor is symmetric, so transpose goes away?
    %Random Question:  Do they ever use FFT to solve these?
    
    %u_star = u + mu*(uxx + uyy) + F1
    %v_star = v + mu*(vxx + vyy) + F2
    
    for i=2:numXCells-1
        for j=2:numYCells-1
            u_star(ind(i,j),t) = u(ind(i,j),t) + (1/(dx^2)) * mu * (u(ind(i-1,j),t) - 4 * u(ind(i,j),t) + u(ind(i+1,j),t) + u(ind(i,j+1),t) + u(ind(i,j-1),t)) + F1(j,i);
            v_star(ind(i,j),t) = v(ind(i,j),t) + (1/(dx^2)) * mu * (v(ind(i-1,j),t) - 4 * v(ind(i,j),t) + v(ind(i+1,j),t) + v(ind(i,j+1),t) + v(ind(i,j-1),t)) + F2(j,i);
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

    A = zeros(numXCells*numYCells,numXCells*numYCells);
    b = zeros(numXCells*numYCells,1);

    %interior points
    for i=2:numXCells-1
        for j=2:numYCells-1

            A(ind(i,j), ind(i-1,j)) = 1;
            A(ind(i,j), ind(i,j)) = -4;
            A(ind(i,j), ind(i+1,j)) = 1;
            A(ind(i,j), ind(i,j+1)) = 1;
            A(ind(i,j), ind(i,j-1)) = 1;

            %b(ind(i,j)) = (dx/2) * (u_star(ind(i+1,j),t) - u_star(ind(i-1,j),t) + v_star(ind(i,j+1),t ) - v_star(ind(i,j-1),t));
            b(ind(i,j)) = (dx/2) * (u_star(ind(i+1,j+1),t) - u_star(ind(i,j+1),t) + ...
                                    u_star(ind(i+1,j),t) - u_star(ind(i,j),t) + ... 
                                    v_star(ind(i+1,j+1),t ) - v_star(ind(i+1,j),t) + ...
                                    v_star(ind(i,j+1),t ) - v_star(ind(i,j),t));
        end
    end

    %boundary conditions
    %Left and right edges

    for j = 2:numYCells-1
        
        %Neumann
        %p_x = 0
        %O(2) forward difference formula for left edge
        %px(1,j) = 0 = -3p(1,j)+4p(2,j)-1p(3,j)
        %O(2) backward difference formula for right edge
        %px(M,j) = 0 = 3p(M,j)-4p(M-1,j)+1p(M-2,j)
        
%         A(ind(1,j), ind(1,j)) = -3;
%         A(ind(1,j), ind(2,j)) = 4;
%         A(ind(1,j), ind(3,j)) = -1;
%         b(ind(1,j)) = 0;
%         
%         A(ind(numXCells,j), ind(numXCells,j)) = 3;
%         A(ind(numXCells,j), ind(numXCells-1,j)) = -4;
%         A(ind(numXCells,j), ind(numXCells-2,j)) = 1;
%         b(ind(numXCells,j)) = 0;

        A(ind(1,j), ind(1,j)) = 1;
        b(ind(1,j)) = 0;
        A(ind(numXCells,j), ind(numXCells,j)) = 1;
        b(ind(numXCells,j)) = 0;
    end

    %Top and bottom edges
    for i = 1:numXCells

%         %Neumann
%         A(ind(i,1), ind(i,1)) = -3;
%         A(ind(i,1), ind(i,2)) = 4;
%         A(ind(i,1), ind(i,3)) = -1;
%         b(ind(i,1)) = 0;
%         
%         A(ind(i,numYCells), ind(i,numYCells)) = 3;
%         A(ind(i,numYCells), ind(i,numYCells-1)) = -4;
%         A(ind(i,numYCells), ind(i,numYCells-2)) = 1;
%         b(ind(i,numYCells)) = 0;

        A(ind(i,1),ind(i,1)) = 1;
        b(ind(i,1)) = 0;
        A(ind(i,numYCells),ind(i,numYCells)) = 1;
        b(ind(i,numYCells)) = 0;
        
    end
 
    %p(:, t) = pinv(A) * b;
    
    if debug
        
        figure();
        surf(X,Y,b(ind(:,:)));
        message = sprintf('b, t = %f',t);
        title(message);
        
    end
    
    p(:,t) = A\b;
    cond(A);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%Step 3 of algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
    %u_next = u_star - grad(p) 
    
    %Question:  Why can't we just set u and v to zero at the borders?
    
    if t ~= numTimeSteps
     
        for i=1:numXCells
            for j=1:numYCells 

                if j==1 || j==numYCells || i==1 || i==numXCells
                    v(ind(i,j),t+1) = 0; 
                    u(ind(i,j),t+1) = 0; 
                else
                     u(ind(i,j),t+1) = u_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i,j),t) - p(ind(i-1,j),t) + p(ind(i,j-1),t) - p(ind(i-1,j-1),t));
                     v(ind(i,j),t+1) = v_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i,j),t) - p(ind(i,j-1),t) + p(ind(i-1,j),t) - p(ind(i-1,j-1),t));

                      %u(ind(i,j),t+1) = (1/(2*dx)) * (p(ind(i,j),t) - p(ind(i-1,j),t) + p(ind(i,j-1),t) - p(ind(i-1,j-1),t));
                      %v(ind(i,j),t+1) = (1/(2*dx)) * (p(ind(i,j),t) - p(ind(i,j-1),t) + p(ind(i-1,j),t) - p(ind(i-1,j-1),t));
                end
                
%                 if j==1
%                     v(ind(i,j),t+1) = v_star(ind(i,j),t) - (1/(dx)) * (p(ind(i,j+1),t) - p(ind(i,j),t));
%                     %v(ind(i,j),t+1) = 0;     
%                 elseif j==numYCells
%                     v(ind(i,j),t+1) = v_star(ind(i,j),t) - (1/(dx)) * (p(ind(i,j),t) - p(ind(i,j-1),t));
%                     %v(ind(i,j),t+1) = 0;
%                 elseif i==1 || i==numXCells
%                     v(ind(i,j),t+1) = v_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i,j+1),t) - p(ind(i,j-1),t));
%                 else
%                     %v(ind(i,j),t+1) = v_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i,j+1),t) - p(ind(i,j-1),t));
%                     v(ind(i,j),t+1) = v_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i,j),t) - p(ind(i,j-1),t) + p(ind(i-1,j),t) - p(ind(i-1,j-1),t));
%                 end
%                 
%                 if i==1 
%                     u(ind(i,j),t+1) = u_star(ind(i,j),t) - (1/(dx)) * (p(ind(i+1,j),t) - p(ind(i,j),t));
%                     %v(ind(i,j),t+1) = 0;
%                 elseif i==numXCells
%                     u(ind(i,j),t+1) = u_star(ind(i,j),t) - (1/(dx)) * (p(ind(i,j),t) - p(ind(i-1,j),t));
%                     %v(ind(i,j),t+1) = 0;
%                 elseif j==1 || j==numYCells
%                     u(ind(i,j),t+1) = u_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i+1,j),t) - p(ind(i-1,j),t));
%                 else
%                     %u(ind(i,j),t+1) = u_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i+1,j),t) - p(ind(i-1,j),t));
%                     u(ind(i,j),t+1) = u_star(ind(i,j),t) - (1/(2*dx)) * (p(ind(i,j),t) - p(ind(i-1,j),t) + p(ind(i,j-1),t) - p(ind(i-1,j-1),t));
%                 end
                       
                
            end    
        end        
    end
    
    tmp1 = p(:, t);
    tmp2 = u(:, t);
    tmp3 = v(:, t);

    %We must transpose here because in matlab rows will map to Y values and
    %columns will map to X values, which is the opposite of how we have
    %defined it in our code.
    
    P(:,:,t) = tmp1(ind)';
    U(:,:,t) = tmp2(ind)';
    V(:,:,t) = tmp3(ind)';
    
end
    
    
    
    