clear all;
close all;
clc;

animate = 0;
checkError = 1;
Linf = 1;

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
timeSteps = 20;

quivRes = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Stokes2DPCwith a high discretization and use it as the "analytical" solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%num / 2 must be even
num = 40;
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
        
        figure();
        surf(X,Y,U(:,:,i));
        message = sprintf('U, t = %f', i);
        title(message);

        figure();
        surf(X,Y,V(:,:,i));
        message = sprintf('V, t = %f', i);
        title(message);

        figure();
        surf(X,Y,P(:,:,i));
        message = sprintf('P, t = %f', i);
        title(message);
        drawnow;
        
        %pause();
    end
end

if checkError

    %Compute the L2E with 20 different grid spacing values
    %Outputting NAN when n is odd
    %testPoints = (num/2 - 20): 2 : (num/2);
    testPoints = (num - 22): 2 : (num-2);

    L2EP = zeros(length(testPoints),1);
    L2EU = zeros(length(testPoints),1);
    L2EV = zeros(length(testPoints),1);
    dx = zeros(length(testPoints),1);

    % Compute the L2 error
    %2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
    for i = 1:length(testPoints)
        testPoints(i)
        tic;
        [ p u v x y numXCells numYCells dx(i)] = Stokes2DPC(g, testPoints(i), p0, mu, height, width, R, L, timeSteps, 0);
        
        %[p u v x y numXCells numYCells] = pTest(height, width, R, L, timeSteps, testPoints(i));
        %[p u v x y numXCells numYCells] = pTest(height, width, R, L, timeSteps, 19);
        dx(i) = y(2,1) - y(1,1);
        
        %%%%%Interp must go from coarse to fine, not fine to coarse!!!!!!%%%%%%%

        if Linf
            L2EP(i) = max(max(abs(interp2(x,y,p(:,:,timeSteps),X,Y,'cubic') - P(:,:,timeSteps))));
            L2EU(i) = max(max(abs(interp2(x,y,u(:,:,timeSteps),X,Y,'cubic') - U(:,:,timeSteps))));
            L2EV(i) = max(max(abs(interp2(x,y,v(:,:,timeSteps),X,Y,'cubic') - V(:,:,timeSteps))));
        else
            L2EP(i) = L2EP(i) + sum(sum((interp2(x,y,p(:,:,timeSteps),X,Y,'cubic') - P(:,:,timeSteps)).^2));
            L2EU(i) = L2EP(i) + sum(sum((interp2(x,y,u(:,:,timeSteps),X,Y,'cubic') - U(:,:,timeSteps)).^2));
            L2EV(i) = L2EP(i) + sum(sum((interp2(x,y,u(:,:,timeSteps),X,Y,'cubic') - V(:,:,timeSteps)).^2));

            L2EP(i) = sqrt(L2EP(i) / (numXCells * numYCells));
            L2EU(i) = sqrt(L2EU(i) / (numXCells * numYCells));
            L2EV(i) = sqrt(L2EV(i) / (numXCells * numYCells)); 
        end
        
        toc;
    end

    figure();
    surf(x,y, p(:,:,timeSteps) - interp2(X,Y,P(:,:,timeSteps),x,y,'cubic'));
    %surf(x,y, p(:,:,timeSteps) - P(1:2:end,1:2:end,timeSteps));
    title('Localized Error for P (Pressure)');
    print('LocalPError', '-djpeg');

    save('data.mat','dx','L2EP','L2EU','L2EV');
    figure()
    loglog(dx,L2EP,'-',dx,dx.^2,'--');
    title('Error for P (Pressure)');
    print('_L2EP', '-djpeg');

    %loglog(d,L2EP,'-',d,L2EP(1)-d(1) + d,'r--',d,L2EP(1)-d(1)^2 + d.^2,'g--');
    %title('L2 Error for P (Pressure)');

    figure()
    loglog(dx,L2EU,'-',dx,dx.^2,'--');
    title('Error for U (Horizontal Velocity)');
    print('_L2EU', '-djpeg');

    figure()
    loglog(dx,L2EV,'-',dx,dx.^2,'--');
    title('Error for V (Vertical Velocity)');
    print('_L2EV', '-djpeg');
end
