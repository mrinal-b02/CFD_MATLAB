% In this code we will save the values at nodes at each timestep and then
% the evolution of the system in the domain is depicted

clear all
close all
clc
%% Meshing
n_points = 31;
dom_size = 1;
dt = 0.0001;
h = dom_size/(n_points - 1) ;
alpha = dt/(h*h) ;

%% Initializing
y(n_points,n_points) = 0;
y(1,:) = 1;

y_new(n_points,n_points) = 0;
y_new(1,:) = 1;

y_transient = 0; % this is the new variable that stores the unique values at mesh nodes at each timestep

error_mag = 1;
error_req = 1e-6;
iterations = 0;

error_track = 0; % tracks how error changes with time by saving values at each iteration

%% Calculations
while error_mag > error_req
    for i = 2:(n_points-1)
        for j = 2:(n_points-1)
            y_new(i,j) = y(i,j) + alpha.*(y(i+1,j) + y(i-1,j) + y(i,j+1) + y(i,j-1) - 4.*(y(i,j)));
            y_transient(iterations+1, 1:n_points, 1:n_points) = y_new;
        end
    end
    iterations = iterations + 1;


    error_mag = 0; % Reset error magnitude for the next iteration
    for i = 2:(n_points-1)
        for j = 2:(n_points-1)
            error_mag = error_mag + abs(y(i,j) - y_new(i,j));
            error_track(iterations) = error_mag;
        end    
    end

    if rem(iterations,100) == 0
        iterations
        error_mag
    end

   y = y_new;
end

%% Plotting
x_dom = ((1:n_points)-1).*h ;
y_dom = 1 - ((1:n_points) - 1).*h ;
[X,Y] = meshgrid (x_dom,y_dom);
contourf(X,Y,y,15)
colorbar

figure;
time = dt.*(1:iterations);
plot(time,error_track)

%% Subplots
figure;
subplot(2,1,1)
plot(time,error_track)
xlabel('time')
ylabel('error')
subplot(2,1,2)
contourf(X,Y,y,15)
colorbar
xlabel('x')
ylabel('y')

%% Plotting a timestep
figure;
timestep_selected = 900;
x_dom = ((1:n_points)-1).*h ;
y_dom = 1 - ((1:n_points) - 1).*h ;
[X,Y] = meshgrid (x_dom,y_dom);
y_timestep = y_transient(timestep_selected, :, :);
y_timestep = reshape(y_timestep,[n_points, n_points]);
contourf(X,Y,y_timestep,12)
colorbar
title (['Time =' num2str(timestep_selected * dt) 's']);

%% Animation at N timesteps
N = 100;
timestep_array = 1:N:iterations;
figure;
for k = 1:length(timestep_array)
    timestep_selected = timestep_array(k);
    y_timestep = y_transient(timestep_selected, :, :);
    y_timestep = reshape(y_timestep, [n_points, n_points]);
    contourf(X, Y, y_timestep, 12)
    colorbar
    title(['Time = ' num2str(timestep_selected * dt) 's'])
    pause(0.5) % Pause for a brief moment to create animation effect
end

%% PLotting the temperature variation at centere of domain
y_centre = y(:, (n_points+1)/2);
figure; 
plot(y_centre)