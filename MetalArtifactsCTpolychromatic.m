close all
clear all
 
clc
load spectrum.txt
energies = spectrum(:,1); weights = spectrum(:,2);

radonpoly = [];

for energyindex = 1:length(energies)
    
    energyindex
    
    energy = energies(energyindex);
    weight = weights(energyindex);
    
    mu_water = attenuation(energy,'Water');
    mu_bone = attenuation(energy,'Bone');
    mu_iron = attenuation(energy,'Iron');
    mu_aluminum = attenuation(energy,'Aluminum');
    
    e= [  1   .9    .9   	    0      	0         90
        -.8  .8    .8   	    0      	0         90
        
        % Small metal circles
        
        mu_aluminum  	.03 	.03   	0     	.35         0 % metal inside write region
        mu_aluminum      .03 	.03    -.22    0          90 % metal inside back region
        mu_aluminum      .03 	.03     .22    0          90 % metal inside back region
        mu_aluminum  	.03  	.03     0     -.35         0]; % metal on bottom
    
    
    theta1 =0:179;
    
    phantom_original = double(phantom(e,512));
    
    %figure,imshow(phantom_original),title('Phantom with metal');%('original with metal non recon')
    
    [Radonphantom,t]=radon(phantom_original,theta1);
    
    if ~exist('ioveri0')
        ioveri0 = zeros(size(Radonphantom));
    end
    
    ioveri0 = ioveri0 + weights(energyindex)*exp(-Radonphantom);
      
end

measured_projections = -log(ioveri0);


figure,imagesc(measured_projections)
set(gca,'YDir','normal')
xlabel('\theta (degrees)')
ylabel('t''')
colormap(gray),title('Measured Projections')


%https://ruiminpan.wordpress.com/2016/03/10/the-curious-case-of-poisson-noise-and-matlab-imnoise-command/
%noisy_measured_projections = 1e11*imnoise(1e-11*measured_projections,'poisson');

noisy_measured_projections = imnoise(measured_projections,'gaussian');


polyfbp = iradon(measured_projections, theta1,'linear','shepp-logan',0.9);
polyfbp2 = iradon(noisy_measured_projections, theta1,'linear','shepp-logan',0.9);

% doublephantom=double(phantom_original);
% MetalPixels=doublephantom >1.3;
% NonMetalPixels=doublephantom <1.3;
% 
% MetalPart=MetalPixels.*doublephantom;  % metal parts of phantom
% NonMetalPart=NonMetalPixels.*doublephantom;  % not metal
% 
% [RadonMetal,t2]=radon(MetalPart,theta1);
% %[RadonNonmetal,t3]=radon(NonMetalPart,theta1);
% 
% figure,imshow(MetalPart),title('Metal Part of Phantom');
% 
% maxRadonPhantom=max(max(Radonphantom));
% 
% corruptedSinogram = Radonphantom;
% replaceMetalWithValue = 0.5*maxRadonPhantom;
% 
% corruptedSinogram(RadonMetal ~= 0) =  replaceMetalWithValue;
% 
% 
% %Image generate artifacts
% 
% artifacts = iradon(corruptedSinogram, theta1,'linear','shepp-logan',0.9);

figure,imshow(polyfbp),title(' metal artifact reconstruction');
set(gca,'CLim',[0,0.5])

figure,imshow(polyfbp2),title(' metal artifact reconstruction');
set(gca,'CLim',[0,0.5])

