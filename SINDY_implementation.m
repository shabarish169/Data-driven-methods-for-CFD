% Clean workspace
clc;
clear;
close all;

addpath("C:\Users\shaba\Downloads\DATA\DATA")
load dataforsindy.mat

%% Plot chaotic attractor
patch(xsol(1,:), xsol(2,:), xsol(3,:), t, 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1); % rainbow coloured trajectory!
xlabel('$x_1$', 'Interpreter', 'Latex')
ylabel('$x_2$', 'Interpreter', 'Latex')
zlabel('$x_3$', 'Interpreter', 'Latex')
set(gca, 'Fontsize', 16)
axis tight
view(45, 45)
grid on



%% Estimate derivatives
dt = t(2)-t(1);
dxdt = (xsol(:,3:end) - xsol(:,1:end-2)) / (2*dt);
x = xsol(:,2:end-1);

%% Create Theta matrix from monomials up to degree 2
Theta = ones([1, length(x(1, :))]); % constant term
Theta(2:4, :) = x; % linear term
Theta(5, :) = x(1, :).^2; % x^2 term
Theta(6, :) = x(1, :).*x(2, :); % xy term
Theta(7, :) = x(1, :).*x(3, :); % xz term
Theta(8, :) = x(2, :).^2; % y^2 term
Theta(9, :) = x(2, :).*x(3, :); % yz term
Theta(10, :) = x(3, :).^2; % z^2 term

%% Find coefficients using pseudoinverse
Xi = dxdt * pinv(Theta);

%% Apply SINDy algorithm to obtain sparse model
lam = 0.1; % Sparsity parameter
k = 1;
Xi_new = Xi;
while sum(sum(abs(Xi - Xi_new))) > 0 || k == 1 
    Xi = Xi_new;
    
    % Find small coefficients
    smallinds = (abs(Xi) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi_new(smallinds) = 0;  
    
    % Loop over all 3 variables of the Lorenz system
    for ind = 1:3 
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(ind, :); 
        
        % Find coefficients for reduced library
        Xi_new(ind, biginds) = dxdt(ind, :) * pinv(Theta(biginds, :)); 
    end
    
    k = k + 1;
end

%% Print model

% Clear command window
clc

% Monomials up to degree 2
mons2 = ["", "x", "y", "z", "x^2", "xy", "xz", "y^2", "yz", "z^2"];
XiString = string(Xi_new);

fprintf('Discovered model using SINDy: \n\n')

% Print dx/dt model:
fprintf('dx/dt = ')
bigcoeffs = abs(Xi_new(1, :)) > 1e-5; % chosen small just to weed out zero coefficients
for jnd = 1:length(bigcoeffs)
    if bigcoeffs(jnd) == 1
        % Print the model by excluding zeroed out terms 
        term = strcat(XiString(1, jnd), mons2(jnd));
        if Xi_new(1, jnd) < 0  
            fprintf('- %s ', term);
        else
            fprintf('+ %s ', term);
        end
    end
end
fprintf('\n')

% Print dy/dt model:
fprintf('dy/dt = ')
bigcoeffs = abs(Xi_new(2, :)) > 1e-5;
for jnd = 1:length(bigcoeffs)
    if bigcoeffs(jnd) == 1
        term = strcat(XiString(2, jnd), mons2(jnd));
        if Xi_new(2, jnd) < 0  
            fprintf('- %s ', term);
        else
            fprintf('+ %s ', term);
        end
    end
end
fprintf('\n')

% Print dz/dt model:
fprintf('dz/dt = ')
bigcoeffs = abs(Xi_new(3, :)) > 1e-5;
for jnd = 1:length(bigcoeffs)
    if bigcoeffs(jnd) == 1
        term = strcat(XiString(3, jnd), mons2(jnd));
        if Xi_new(3, jnd) < 0  
            fprintf('- %s ', term);
        else
            fprintf('+ %s ', term);
        end
    end
end
fprintf('\n')
