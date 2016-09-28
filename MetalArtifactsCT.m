close all
clear all
clc

e= [ 1  .9    .9   	    0      	0         90
    -.8  .8    .8   	    0      	0         90
    %1    .920    .690   	    0      	0         90
    %-.8 	.874 	.6624   	0      -.0184     90
    %     -.2 	.310 	.1100     .22    	  0       72
    %     -.2 	.410 	.1600    -.22    	  0      108
    %     .1  	.250 	.2100   	0     	.35       90
    %     .1  	.0460 	.046   	    0     	.1         0
    %     .1  	.0460 	.046  	    0      -.1         0
    %     .1  	.0460 	.023     -.08      -.605       0
    %     .1  	.0230 	.023    	0      -.605       0
    %     .1  	.046  	.023   	  .06      -.605      90
    
      
    % Small metal circles
    
         10  	.018 	.018   	0     	.35         0 % metal inside write region
         10      .018 	.018    -.22    0         90 % metal inside back region
         10      .018 	.018     .22    0         90 % metal inside back region
         10  	.018  	.018     0     -.35          0]; % metal on bottom
    
    
    %{
% for big circle
    10  	.05 	.05   	0     	.35         0 % metal inside write region
    10      .05 	.05    -.22    0         90 % metal inside back region
    10      .05 	.05     .22    0         90 % metal inside back region
    10  	.05  	.05     0     -.35          0]; % metal on bottom
    %}
    %{
 % for small oval
    10  	.056 	.023   	0     	.35         180 % metal inside write region
    10      .056 	.023    -.22    0         90 % metal inside back region
    10      .056 	.023     .22    0         90 % metal inside back region
    10  	.056  	.023     0     -.45          180]; % metal on bottom
    %}
    %{
 % for big oval
    10  	.076 	.023   	0     	.35         180 % metal inside write region
    10      .076 	.023    -.22    0         90 % metal inside back region
    10      .076 	.023     .22    0         90 % metal inside back region
    10  	.076  	.023     0     -.45          180]; % metal on bottom
    %}   
    %{
 % for big oval
 
    10      .076 	.026    -.22    0         90 % metal inside back region

    10  	.076  	.026     0     -.75          180]; % metal on bottom
    
    %}
    
    
    
    theta1 =0:1:179;
    
    phantom_original = phantom(e,512);
    
    figure,imshow(phantom_original),title('Phantom with metal');%('original with metal non recon')
    
    %theta1 =0:1:359;%(180 projection) if used different in 2 degree ; 0:2:178(90 projection)
    
    [Radonphantom,t]=radon(phantom_original,theta1);
    
    figure,imshow(Radonphantom,[],'Xdata',theta1,'Ydata',t,...
    'InitialMagnification','fit')
    xlabel('\theta (degrees)')
    ylabel('t''')
    colormap(gray),title('Phantom')
    
   
    doublephantom=double(phantom_original);
    MetalPixels=doublephantom >1.3;
    NonMetalPixels=doublephantom <1.3;
   
   
    MetalPart=MetalPixels.*doublephantom;  % metal parts of phantom
    NonMetalPart=NonMetalPixels.*doublephantom;  % not metal
    
    [RadonMetal,t2]=radon(MetalPart,theta1);
    %[RadonNonmetal,t3]=radon(NonMetalPart,theta1);
    
    figure,imshow(MetalPart),title('Metal Part of Phantom');
    
    
    maxRadonPhantom=max(max(Radonphantom));
           
    corruptedSinogram = Radonphantom; 
    replaceMetalWithValue = 0.5*maxRadonPhantom;
    
    corruptedSinogram(RadonMetal ~= 0) =  replaceMetalWithValue;
    
   
    %Image generate artifacts
    
    artifacts = iradon(corruptedSinogram, theta1,'linear','shepp-logan',0.9);
    
    figure,imshow(artifacts),title(' metal artifact reconstruction');
    
    
