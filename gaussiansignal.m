function gaussiansignal

close all

x = linspace(-15,15,512);
y = x;
w = 3;
[X,Y] = meshgrid(x,y);
centergaussian = [0, 0];

x0 = centergaussian(1);
y0 = centergaussian(2);

gaussian = (1/w^2)*exp(-pi*((X-x0).^2 + (Y-y0).^2)./w^2);

imagesc(x,y,gaussian);
set(gca,'YDir','normal');
title('Gaussian')

[numericalR, numericalt] = radon(gaussian);

theta_deg = linspace(0,180,256);
theta_rad = theta_deg*(pi/180);
t = numericalt;
[THETA_RAD, T] = meshgrid(theta_rad,t);

analyticradongaussian = gaussian_projection_shifted(T,THETA_RAD,w,centergaussian);

figure

subplot(2,1,1); imagesc(T, THETA_RAD, analyticradongaussian); colorbar; title('Analytic Radon')
subplot(2,1,2); imagesc(numericalR); colorbar; title('Numerical Radon')

figure
fbp_gaussian = iradon(analyticradongaussian, theta_deg);
fbp_numerical_gaussian = iradon(numericalR,0:179);

subplot(2,1,1); imagesc(fbp_gaussian); colorbar
subplot(2,1,2); imagesc(fbp_numerical_gaussian); colorbar


function p = gaussian_projection(t, w)
%  Projection of Gaussian centered at origin
p = (1/w)*exp(-pi*(t/w).^2);

function p = gaussian_projection_shifted(t, theta, w, center_gaussian)
%  Projection of shifted disk at (x0, y0)
%  Shift property:  f(x-x0, y-y0) -> p(t - x0*cos(theta) - y0*sin(theta))

x0 = center_gaussian(1,1);
y0 = center_gaussian(1,2);

p = gaussian_projection(t - x0.*cos(theta) - y0.*sin(theta), w);



