function [flippedlumps, lumpprojections] = lumpy_phantom_validation(x,y,energy,theta_deg)

%  Generate lumps
%  Define parameters for the lumps.

numberofgaussians = 100;
insetradius = 10;
lumpcenters = -10 + (10 + 10)*rand(numberofgaussians,2); % 
r = sqrt(lumpcenters(:,1).^2 + lumpcenters(:,2).^2);

lumpcenters = lumpcenters(r<insetradius,:);
numberofgaussians = length(lumpcenters);

widths = 1 + 3*rand(numberofgaussians,1);
muvalues = 0.1*(-5 + (5 + 5))*attenuation(energy, 'Bone')*rand(numberofgaussians);%[1; 1; 1];

subplot(2,2,1)
lumps = generate_lumps(x,y, numberofgaussians, widths, lumpcenters, muvalues);

%  Numerical Radon Transform

subplot(2,2,3)
[numericalR, numericalt] = matlabtransform(x,y,theta_deg,lumps);  % Radon transform via MATLAB

lumpprojections = numericalR;
 
%  Analytic Radon Transform

subplot(2,2,2)
t  = numericalt;  %  t-values at which to evaluate analytic taken from numerical results
%t = linspace(-3,3,180);
theta_rad = theta_deg*(pi/180);
[THETA_RAD, T] = meshgrid(theta_rad,t);
%radontransform = zeros(size(theta_rad));

lumpprojections = zeros(size(T));
for i = 1:numberofgaussians
    lumpprojections = lumpprojections + gaussian_projection_shifted(T, THETA_RAD,...
        muvalues(i), widths(i), lumpcenters(i,:));
end

imagesc(theta_deg,t,lumpprojections);
xlabel('\theta');
ylabel('t');
title('Analytic');
set(gca,'YDir','normal')
axis square
colorbar

%  Matlab FBP
subplot(2,2,4)
fbp = matlab_fbp(lumpprojections,theta_deg,x,y);


%  There is a bug, and iradon returns upside-down lumps.  In this case it
%  doesn't matter so I am patching by flipping the initial lumps to match
%  the projections.

lumps = flipud(lumps);
subplot(2,2,1)
imagesc(x,y,lumps);
set(gca,'YDir','normal');
axis square;

function p = gaussian_projection(t, w)
%  Projection of Gaussian centered at origin
p = (1/w)*exp(-pi*(t/w).^2);

function p = gaussian_projection_shifted(t, theta, mu, w, center_gaussian)
%  gaussian_projection of shifted Gaussian at (x0, y0)
%  Shift property:  f(x-x0, y-y0) -> p(t - x0*cos(theta) - y0*sin(theta))

x0 = center_gaussian(1,1);
y0 = center_gaussian(1,2);

p = gaussian_projection(t - x0.*cos(theta) - y0.*sin(theta), w);

function lumpy_phantom = generate_lumps(x,y, numberofgaussians, widths, centers, muvalues)

[X,Y] = meshgrid(x,y);

lumpy_phantom = zeros(size(X));

for i = 1:numberofgaussians
    w = widths(i);
    center_gaussian = centers(i,:);
    strength_gaussian = muvalues(i);
    RminusRiSq = (X-center_gaussian(1)).^2 + (Y-center_gaussian(2)).^2;
    gaussian = (1/w^2)*exp(-pi*RminusRiSq./w^2);
    lumpy_phantom = lumpy_phantom + gaussian;
end
imagesc(x, y, lumpy_phantom); axis square; set(gca,'YDir','normal');
colorbar;

function fbp_cm = matlab_fbp(radontransform,theta_deg,x,y)
L = max(x) - min(x);
N = length(x);
cm_per_pixel = L/N;

fbp_pixels = iradon(radontransform,theta_deg);
fbp_cm = fbp_pixels / cm_per_pixel;  % Just reverse previous scaling

imagesc(x, y,  fbp_cm); axis square;...
    set(gca,'YDir','normal'), title('Filtered backprojection');
colorbar

function [R_cm,t_cm] = matlabtransform(x,y,theta_deg,phantom)

%  Matlab indexing is upside down.
%  N pixels = L cm
%  Matlab lengths (including integrals) are in pixels and have to be multiplied by L/N.

phantom = flipud(phantom);

L = max(x) - min(x);
N = length(x);
cm_per_pixel = L/N;
[R_pixels,t_pixels] = radon(phantom, theta_deg);
t_cm = t_pixels*cm_per_pixel;
R_cm = R_pixels*cm_per_pixel;
imagesc(theta_deg,t_cm,R_cm);
title('Numerical')
xlabel('\theta (degrees)');
ylabel('t');
colormap(hot);
set(gca,'YDir','normal')
colorbar
axis square;
