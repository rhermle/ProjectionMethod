clear all;
close all;
clc;

animate = 1;
checkL2E = 1;

%Pressure Differential (Right side)
p0 = 0;

%Viscosity
mu = 1;

%grav. constant
g = 0;

height = 20;
width = 20;
R = 5;
L = 5;
%T = .25;
timeSteps = 1;

quivRes = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Stokes2DPCwith a high discretization and use it as the "analytical" solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%num / 2 must be even
num = 32;
tic;
%[P U V X Y NUMXCELLS NUMYCELLS] = Stokes2DPC(g, num, p0, mu, height, width, R,L, timeSteps, 0);

[P U V X Y NUMXCELLS NUMYCELLS] = pTest(height, width, R, L, timeSteps, num);

toc;

% i=2;
% figure(6);
% hold on;
% quiver(X(1:quivRes:end,1:quivRes:end),Y(1:quivRes:end,1:quivRes:end),U(1:quivRes:end,1:quivRes:end, i),V(1:quivRes:end,1:quivRes:end, i));
% z = (X - (R + L)).^2 + Y.^2 - R^2;
% contour(X,Y,z,[0,0]);
% streamline(X,Y,U(:,:,i),V(:,:,i),.1,1);
% streamline(X,Y,U(:,:,i),V(:,:,i),.1,-1);
% % streamline(X,Y,U,V,-width / 2, 0.1);
% % streamline(X,Y,U,V,-width / 2, -0.1);
% drawnow;
% hold off;


%surf(X,Y,U(:,:,2));

%pcolor?
%shading interp

if animate
   
    for i = 1:timeSteps

        figure();
        clf;
        hold on;

        quiver(X(1:quivRes:end,1:quivRes:end),Y(1:quivRes:end,1:quivRes:end),U(1:quivRes:end,1:quivRes:end, i),V(1:quivRes:end,1:quivRes:end, i));
        z = (X - (R + L)).^2 + Y.^2 - R^2;
        contour(X,Y,z,[0,0]);

        streamline(X,Y,U(:,:,i),V(:,:,i),.1,1);
        streamline(X,Y,U(:,:,i),V(:,:,i),.1,-1);
        % streamline(X,Y,U,V,-width / 2, 0.1);
        % streamline(X,Y,U,V,-width / 2, -0.1);
        
        message = sprintf('UV Quiver, t = %f', i);
        title(message);
        
        drawnow;
        hold off;
        
%         figure();
%         surf(X,Y,U(:,:,i));
%         message = sprintf('U, t = %f', i);
%         title(message);
% 
%         figure();
%         surf(X,Y,V(:,:,i));
%         message = sprintf('V, t = %f', i);
%         title(message);
% 
        figure();
        surf(X,Y,P(:,:,i));
        message = sprintf('P, t = %f', i);
        title(message);
        drawnow;
        
        %pause();
    end
end

if checkL2E

    %Compute the L2E with 20 different grid spacing values
    %Outputting NAN when n is odd
    %testPoints = (num/2 - 20): 2 : (num/2);
    testPoints = (num - 20): 2 : (num);

    L2EP = zeros(10,1);
    L2EU = zeros(10,1);
    L2EV = zeros(10,1);
    dx = zeros(10,1);

    % Compute the L2 error
    %2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
    for i = 1:10
        testPoints(i)
        tic;
        [ p u v x y numXCells numYCells dx(i)] = Stokes2DPC(g, testPoints(i), p0, mu, height, width, R, L, timeSteps, 0);
        
        %d(i) = (1-0) / (i - 1);

        for j=1:numXCells
            for k = 1:numYCells
                % interp2(X,Y,F,x,y) interpolates F(X,Y) at x,y
                % Don't Want:  X is numx1 vector, Y is a numx1, P is a num^2x1 vector
                % Want: X is numxnum array, Y is numxnum array, P is a numxnum array
                L2EP(i) = L2EP(i) + (p(k,j,timeSteps) - interp2(X,Y,P(:,:,timeSteps),x(k,j),y(k,j), 'cubic'))^2;
                L2EU(i) = L2EU(i) + (u(k,j,timeSteps) - interp2(X,Y,U(:,:,timeSteps),x(k,j),y(k,j), 'cubic'))^2;
                L2EV(i) = L2EV(i) + (v(k,j,timeSteps) - interp2(X,Y,V(:,:,timeSteps),x(k,j),y(k,j), 'cubic'))^2;
            end
        end   
        L2EP(i) = sqrt(L2EP(i) / (numXCells * numYCells));
        L2EU(i) = sqrt(L2EU(i) / (numXCells * numYCells));
        L2EV(i) = sqrt(L2EV(i) / (numXCells * numYCells));
        
        toc;
    end

    figure();
    surf(x,y, p(:,:,timeSteps) - interp2(X,Y,P(:,:,timeSteps),x,y));
    title('Localized Error for P (Pressure)');
    print('LocalPError', '-djpeg');

    % dx = dx(testPoints);
    % L2EP = L2EP(testPoints);
    % L2EU = L2EU(testPoints);
    % L2EV = L2EV(testPoints);

    save('data.mat','dx','L2EP','L2EU','L2EV');
    figure()
    loglog(dx,L2EP,'-',dx,dx.^2,'--');
    title('L2 Error for P (Pressure)');
    print('_L2EP', '-djpeg');

    %loglog(d,L2EP,'-',d,L2EP(1)-d(1) + d,'r--',d,L2EP(1)-d(1)^2 + d.^2,'g--');
    %title('L2 Error for P (Pressure)');

    figure()
    loglog(dx,L2EU,'-',dx,dx.^2,'--');
    title('L2 Error for U (Horizontal Velocity)');
    print('_L2EU', '-djpeg');

    figure()
    loglog(dx,L2EV,'-',dx,dx.^2,'--');
    title('L2 Error for V (Vertical Velocity)');
    print('_L2EV', '-djpeg');
end
