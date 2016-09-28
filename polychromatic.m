function [phantom, fbp] = polychromatic

load spectrum.txt
energies = spectrum(:,1); weights = spectrum(:,2);
phantom_parameters = set_phantom_parameters(50);
phantom = generate_phantom(phantom_parameters);

phantom = phantom + lumpy


x = ((0:399)-199.5)*0.075;  %  Dimensions taken from Wunderlich, Noo, PMB paper
y = x;

for energyindex = 20:20
    %1:length(energies)
    
    energyindex
    
    %for energyindex = 1:length(energies)
    
    energy = energies(energyindex);
    weight = weights(energyindex);
    
    radontransform = monoenergetic(energy);
     
    
    if ~exist('ioveri0')
        ioveri0 = zeros(size(radontransform));
    end
    
    ioveri0 = ioveri0 + weights(energyindex)*exp(-radontransform);
    
end

measured_projections = -log(ioveri0);

%measured_projections = addnoise(measured_projections);



%  Matlab FBP
subplot(2,2,4)
theta_deg = 0:179;
fbp = matlab_fbp(measured_projections,theta_deg,x,y);

close all;

figure; colormap bone;
subplot(1,2,1)
imagesc(x, y, phantom)
set(gca,'YDir','normal');
axis square;
ax = get(gca,'CLim');
colorbar;
title('Phantom')

subplot(1,2,2)
imagesc(x, y, fbp)
set(gca,'CLim',ax);
axis square;
colorbar;
title('FBP with polychromatic')


figure
plot(energies, weights, '.')
xlabel('keV');
ylabel('Photons/area/time (normalized)');



function phantom_parameters = set_phantom_parameters(energy)

x = ((0:399)-199.5)*0.075;  %  Dimensions taken from Wunderlich, Noo, PMB paper
y = x;

ADD_NOISE = 0;
TRUNCATE_PROJECTIONS = 0;

%energy = 110;  %  In keV?

outerdiskradius = 12;
outerdiskcenter = [0,0];
innerdiskradius = 0;
innerdiskcenter = [0,0];

disk_angles = linspace(0, 2*pi,5).';
disk_angles = disk_angles(1:length(disk_angles)-1);
numberofsmalldisks = length(disk_angles);

smalldiskcenters = mean([outerdiskradius, innerdiskradius]).*[cos(disk_angles), sin(disk_angles)];
centers = [outerdiskcenter; innerdiskcenter; smalldiskcenters];
radii = [outerdiskradius,innerdiskradius, .4*ones(1,numberofsmalldisks)];

numberofdisks = length(centers);

mu_water = attenuation(energy,'Water');
mu_bone = attenuation(energy,'Bone');
mu_aluminum = attenuation(energy,'Aluminum');

%  Outer ring, inner ring, small disks

% The phantom is set up to sum disks, because the analytic radon transform
% is set up that way.  Be careful with disk overlap, etc.

muvalues = [mu_bone - mu_water, 0, mu_aluminum*ones(1,numberofsmalldisks) - mu_water];

phantom_parameters =...
    struct(...
    'x', x,...
    'y', y,...
    'energy', energy,...
    'centers', centers,...
    'radii', radii,...
    'muvalues', muvalues...
    );

function phantom = generate_phantom(phantom_parameters)

x = phantom_parameters.x;
y = phantom_parameters.y;
radii = phantom_parameters.radii;
centers = phantom_parameters.centers;
numberofdisks = length(centers);
muvalues = phantom_parameters.muvalues;

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

imagesc(x, y, phantom); axis square; set(gca,'YDir','normal');
colorbar;



function noisy_projections = addnoise(radontransform)
noisy_projections = radontransform + 0.05*rand(size(radontransform)).*sqrt(radontransform);

function [R_cm,t_cm] = matlabtransform(x,y,phantom)

%  Matlab indexing is upside down.
%  N pixels = L cm
%  Matlab lengths (including integrals) are in pixels and have to be multiplied by L/N.

phantom = flipud(phantom);

L = max(x) - min(x);
N = length(x);
cm_per_pixel = L/N;
theta_deg = 0:179;
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

function phantom = generate_phantom(x,y, numberofdisks, radii, centers, muvalues)

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
imagesc(x, y, phantom); axis square; set(gca,'YDir','normal');
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

function mu = attenuation(energy,material)

%  This function is from Wunderlich, Noo, PMB paper

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