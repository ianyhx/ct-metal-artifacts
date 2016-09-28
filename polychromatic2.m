function [polyenergetic_phantom, fbp] = polychromatic2

 
%  Bugs and to dos
%{
  
From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3426508/#R36

Note also that the values of mu are to be interpreted as
attenuation factors relative to that of water, so that 1000 (mu - 1)
is in Hounsfield units.

Why does monoenergetic code flip the phantom upside down?  (Temp fixed with extra flipud statement before iradon)
Polyenergetic code does not do this.

Remove all the teeth and see if the intensities of phantom and
reconstruction match...if not, why not?

Find a better model for teeth.

Why do numerical and analytic projections yield different results?

Is there a material less attenuating than iron where we get effects but we
do not break the numerics?  Other cross-sections from NIST.


% Radon transform of Gaussian:
% http://www.aprendtech.com/blog/P15Ctprojsim/P15Ctprojsim.html

%}

load spectrum2.txt

i0 = 19362005.47;  %  Mean number of photons/(mm^2) (integrated over energy)
scatter = 1000;

energies_all = spectrum2(:,1);
weights_all = spectrum2(:,2);
theta_deg = linspace(0,180,300);
theta_deg = theta_deg(1:end-1);  % Make sure 0 = 180 is not repeated.

%selected_spectrum_indices = [20:50];

selected_spectrum_indices = [10, 30, 80, 100];

%selected_spectrum_indices = [10, 15, 20, 25, 30, 52, 60, 80, 100];
plot(energies_all,weights_all); hold on;

energies = energies_all(selected_spectrum_indices);
weights = weights_all(selected_spectrum_indices);

plot(energies, weights, 'ro');

weights = weights/sum(weights);

for energyindex = 1:length(energies)
    
    energyindex
    
    energy = energies(energyindex);
    weight = weights(energyindex);
    
    %    weight = weights(energyindex)./ sum(weights(1:30:length(energies)));
    
    %  Set parameters defining the phantom
    phantom_parameters = set_jaw_phantom_parameters(energy);
    colormap bone;
    
    %  Generate a numerical phantom
    subplot(2,2,1)
    phantom = generate_phantom(phantom_parameters);
    
    
    %  Generate numerical projections via Matlab's radon function
    [numericalRadon, numericalt] = matlabtransform(phantom_parameters.x, phantom_parameters.y, theta_deg, phantom);  % Radon transform via MATLAB
    %radontransform = numericalRadon;
    radontransform = generate_analytic_projections(theta_deg, numericalt, phantom_parameters);
    %   radontransform = generate_analytic_projections(theta_deg, t, phantom_parameters);
    
    
    %  Add lumpy background
    
    [lumps, lump_projections] = addlumps(phantom_parameters.x, phantom_parameters.y,energy,numericalt, theta_deg);
    
    %lumps = 0.1*lumps;
    %lump_projections = 0.1*lump_projections;
    
    phantom = phantom + lumps;
    radontransform = radontransform + lump_projections;
    
    
    if ~exist('polyenergetic_phantom')
        polyenergetic_phantom = zeros(size(phantom));
    end
    
    polyenergetic_phantom = polyenergetic_phantom + weight*phantom;
    
    if ~exist('ioveri0')
        ioveri0 = zeros(size(radontransform));
    end
    
    ioveri0 = ioveri0 + poissrnd(i0*weight*exp(-radontransform)+scatter)/i0;
    
    %    No noise or scatter:
    %    ioveri0 = ioveri0 + weight*exp(-radontransform);
    
end

i = ioveri0*i0 ;

noisy_ioveri0 = poissrnd(i+scatter)/i0;

ioveri0 = noisy_ioveri0;

measured_projections = -log(ioveri0);

idx = isinf(measured_projections);
measured_projections(idx) = -log(eps/i0);

%  Matlab FBP

fbp = matlab_fbp(measured_projections,theta_deg,phantom_parameters.x,phantom_parameters.y);

close all;

figure; colormap bone;
subplot(1,2,1)
imagesc(phantom_parameters.x, phantom_parameters.y, polyenergetic_phantom)
set(gca,'YDir','normal');
axis square;
ax = get(gca,'CLim');
colorbar;
title('Phantom')
colorscale = get(gca, 'Clim');
set(gca,'CLim',[0,.3])


subplot(1,2,2)
imagesc(phantom_parameters.x, phantom_parameters.y, fbp)
set(gca,'CLim',ax);
axis square;
colorbar;
set(gca, 'Clim', colorscale);
title('FBP with polychromatic')
set(gca,'CLim',[0,.3])


figure
plot(energies, weights, '.')
xlabel('keV');
ylabel('Photons/area/time (normalized)');

function radontransform = monoenergetic(energy)

%  Set parameters defining the phantom
phantom_parameters = set_phantom_parameters(energy);
colormap bone;

%  Generate a numerical phantom
subplot(2,2,1)
phantom = generate_phantom(phantom_parameters);
colorscale = get(gca, 'Clim');

%  Generate numerical projections via Matlab's radon function
subplot(2,2,2)
[numericalR, numericalt] = matlabtransform(phantom_parameters.x, phantom_parameters.y,phantom);  % Radon transform via MATLAB

%  Generate analytic projections for comparison
subplot(2,2,3)
theta_deg = 0:179;
%t = linspace(-22,22,1024);
radontransform = generate_analytic_projections(theta_deg, numericalt, phantom_parameters);

% Use Matlab iradon to perform FBP
subplot(2,2,4)
fbp = matlab_fbp(radontransform,theta_deg,phantom_parameters.x,phantom_parameters.y);
set(gca, 'Clim', colorscale);

function phantom = generate_phantom(phantom_parameters)

x = phantom_parameters.x;
y = phantom_parameters.y;
radii = phantom_parameters.radii;
centers = phantom_parameters.centers;
numberofdisks = length(centers);
muvalues = phantom_parameters.mu_additive_values;

[X,Y] = meshgrid(x,y);

phantom = zeros(size(X));
phantomdisk = zeros(size(X));

for i = 1:numberofdisks
    radius_disk = radii(i);
    center_disk = centers(i,:);
    mu_disk = muvalues(i);
    RminusRi = sqrt((X-center_disk(1)).^2 + (Y-center_disk(2)).^2);
    phantomdisk(RminusRi <= radius_disk) =  mu_disk;
    %phantom(phantomdisk ~= 0) = phantomdisk(phantomdisk ~= 0);
    phantom = phantom + phantomdisk;
    phantomdisk = zeros(size(X));
end

function radontransform = generate_analytic_projections(theta_deg, t, phantom_parameters)

% For comparison with numerical results, use the t-values generated by MATLAB's Radon
% t = linspace(-3,3,180);

radii = phantom_parameters.radii;
centers = phantom_parameters.centers;
numberofdisks = length(centers);
muvalues = phantom_parameters.mu_additive_values;

theta_rad = theta_deg*(pi/180);
[THETA_RAD, T] = meshgrid(theta_rad,t);
radontransform = zeros(size(T));

for i = 1:numberofdisks
    %radontransform = projection(T, muvalues(1), radii(1));  %  Disk At origin
    radontransform = radontransform +...
        projection_shifted(T, THETA_RAD, muvalues(i), radii(i), centers(i,:));
end

imagesc(theta_deg,t,radontransform);
xlabel('\theta');
ylabel('t');
title('Analytic');
set(gca,'YDir','normal')
axis square
colorbar

function p = projection(t, mu, radius_disk)
%  Projection of disk centered at origin
p = radius_disk^2 - t.^2;
p(p < 0) = 0;
p = sqrt(p);
p = 2*mu*p;

function p = projection_shifted(t, theta, mu, radius_disk, center_disk)
%  Projection of shifted disk at (x0, y0)
%  Shift property:  f(x-x0, y-y0) -> p(t - x0*cos(theta) - y0*sin(theta))

x0 = center_disk(1,1);
y0 = center_disk(1,2);

p = projection(t - x0.*cos(theta) - y0.*sin(theta), mu, radius_disk);

function [R_cm,t_cm] = matlabtransform(x,y,theta_deg,phantom)

%  Matlab indexing is upside down.
%  N pixels = L cm
%  Matlab lengths (including integrals) are in pixels and have to be multiplied by L/N.

phantom = flipud(phantom);

L = max(x) - min(x);
N = length(x);
cm_per_pixel = L/N;
[R_pixels,t_pixels] = radon(phantom, theta_deg);
%[R_pixels,t_pixels] = radon(phantom,theta_deg);
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

function phantom_parameters = set_jaw_phantom_parameters(energy)

mu_water = attenuation(energy,'Water');
mu_bone = attenuation(energy,'Bone');
mu_aluminum = attenuation(energy,'Aluminum');
mu_iron = attenuation(energy,'Iron');

x = linspace(-15,15,512);
%x = ((0:399)-199.5)*0.075;  %  Dimensions taken from Wunderlich, Noo, PMB paper
y = x;

outerdiskradius = 12;
outerdiskcenter = [0,0];
innerdiskradius = 10;
innerdiskcenter = [0,0];
insetradius = 4.5;
insetcenter = [0,0];

xteeth = -4.5:4.5;
yteeth = 7-(.5*xteeth).^2;
teethcenters = [xteeth.', yteeth.'];
teethradii = 0.5*ones(size(xteeth));
numberofteeth = length(xteeth);

fillingcenters = [teethcenters(2,:); teethcenters(10,:)];
numberoffillings = length(fillingcenters);
fillingradii = 0.2*ones(1,numberoffillings);

%
% %{  Uncomment to remove teeth
% teethcenters = [];
% teethradii = [];
% numberofteeth = 0;
% %}
%
%
% %{  Uncomment to remove fillings
% fillingcenters = [];
% fillingradii = [];
% numberoffillings = 0;
% %}

centers = [outerdiskcenter; innerdiskcenter; insetcenter; teethcenters; fillingcenters];
radii = [outerdiskradius, innerdiskradius, insetradius, teethradii, fillingradii];

mu_outer_disk = mu_bone;
mu_inner_disk = mu_water;
mu_inset      = 1.05*mu_water;
mu_teeth      = mu_bone;
mu_fillings   = mu_iron;

mu_absolute_values =...
    [mu_outer_disk,...
    mu_inner_disk,...
    mu_inset,...
    mu_teeth*ones(1,numberofteeth),...
    mu_fillings*ones(1,numberoffillings)];


mu_additive_disks =...
    [mu_outer_disk,...
    mu_inner_disk - mu_outer_disk,...
    mu_inset - (mu_inner_disk - mu_outer_disk) - mu_outer_disk];

mu_additive_teeth = ones(1,numberofteeth)*(mu_teeth - sum(mu_additive_disks));
mu_additive_fillings = ones(1,numberoffillings)*(mu_fillings - mu_teeth - sum(mu_additive_disks));

mu_additive_values = [mu_additive_disks, mu_additive_teeth, mu_additive_fillings];

phantom_parameters =...
    struct(...
    'x', x,...
    'y', y,...
    'energy', energy,...
    'centers', centers,...
    'radii', radii,...
    'mu_absolute_values', mu_absolute_values,...
    'mu_additive_values', mu_additive_values...
    );


function mu = attenuation(energy,material)

%  This function is from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3426508/#R36
%  "Note that ? and Ep are expressed in cm?1 and keV, respectively."

switch material
    case 'Bone'
        p=[-.179564 2.851439 -16.055087 35.24159 -23.704935];
    case 'Iron'
        p=[0 .238719 -2.479780 5.725218 3.971565];
    case 'Water'
        p=[-.014027 -.0459590 2.366105 -13.683202 21.867818];
    case 'Aluminum'
        p=[-.226414 3.57626 -20.22871 46.52298 -33.57025];
        
    otherwise
        p=[];
end
lep=log(energy);
mu=exp(p(1)*lep.^4+p(2)*lep.^3+p(3)*lep.^2+p(4)*lep+p(5));

function fbp_cm = matlab_fbp(radontransform,theta_deg,x,y)
L = max(x) - min(x);
N = length(x);
cm_per_pixel = L/N;

fbp_pixels = iradon(radontransform,theta_deg);
fbp_cm = fbp_pixels / cm_per_pixel;  % Just reverse previous scaling

imagesc(x, y,  fbp_cm); axis square;...
    set(gca,'YDir','normal'), title('Filtered backprojection');
colorbar

function [lumps, lumpprojections] = addlumps(x,y,energy,t,theta_deg)

numberofgaussians = 15;
insetradius = 4;  
lumpcenters = -insetradius + (insetradius+insetradius)*rand(numberofgaussians,2); %

r = sqrt(lumpcenters(:,1).^2 + lumpcenters(:,2).^2);
lumpcenters = lumpcenters(r<insetradius,:);

numberofgaussians = length(lumpcenters);

widths = 1 + 2*rand(numberofgaussians,1);
muvalues = 1.5*(-1 + (1 + 1))*attenuation(energy, 'Water')*rand(numberofgaussians, 1);

%  Add Gaussian signal

lumpcenters = [lumpcenters; 0, 0];
numberofgaussians = numberofgaussians + 1;
widths = [widths; 1];  % Signal
muvalues = [muvalues; attenuation(energy, 'Bone')];


subplot(2,2,1)
lumps = generate_lumps(x,y, numberofgaussians, widths, lumpcenters, muvalues);

theta_rad = theta_deg*(pi/180);
[THETA_RAD, T] = meshgrid(theta_rad,t);

lumpprojections = zeros(size(T));
for i = 1:numberofgaussians
    lumpprojections = lumpprojections + gaussian_projection_shifted(T, THETA_RAD,...
        muvalues(i), widths(i), lumpcenters(i,:));
end

%  There is a bug, and iradon returns upside-down lumps.  In this case it
%  doesn't matter so I am patching by flipping the initial lumps to match
%  the projections.

lumps = flipud(lumps);

function p = gaussian_projection(t, w)
%  Projection of Gaussian centered at origin
p = (1/w)*exp(-pi*(t/w).^2);

function p = gaussian_projection_shifted(t, theta, mu, w, center_gaussian)
%  gaussian_projection of shifted Gaussian at (x0, y0) with amplitude mu
%  Shift property:  f(x-x0, y-y0) -> p(t - x0*cos(theta) - y0*sin(theta))

x0 = center_gaussian(1,1);
y0 = center_gaussian(1,2);

p = mu*gaussian_projection(t - x0.*cos(theta) - y0.*sin(theta), w);

function lumpy_phantom = generate_lumps(x,y, numberofgaussians, widths, centers, muvalues)

[X,Y] = meshgrid(x,y);

lumpy_phantom = zeros(size(X));

for i = 1:numberofgaussians
    w = widths(i);
    center_gaussian = centers(i,:);
    strength_gaussian = muvalues(i);
    RminusRiSq = (X-center_gaussian(1)).^2 + (Y-center_gaussian(2)).^2;
    gaussian = strength_gaussian*(1/w^2)*exp(-pi*RminusRiSq./w^2);
    lumpy_phantom = lumpy_phantom + gaussian;
end
%imagesc(x, y, lumpy_phantom); axis square; set(gca,'YDir','normal');
colorbar;


