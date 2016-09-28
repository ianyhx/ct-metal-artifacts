%close all
%clear all
%clc

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
    
    D = floor(sqrt(size(phantom_original,1)^2 + size(phantom_original,2)^2))+10;
        
    [fanbeamphantom, Fpos, Fangles] = fanbeam(phantom_original,D,'FanSensorSpacing',0.2);
    
    imshow(fanbeamphantom,[],'XData',Fangles,'YData',Fpos,...
        'InitialMagnification','fit')
    axis normal
    xlabel('Rotation Angles (degrees)')
    ylabel('Sensor Positions (degrees)')

    
    I = 512;
    D = floor(sqrt(size(phantom_original,1)^2 + size(phantom_original,2)^2)) + 10;
    
    doublephantom=double(phantom_original);
    MetalPixels=doublephantom >1.3;
    NonMetalPixels=doublephantom <1.3;
    
    
    MetalPart=MetalPixels.*doublephantom;  % metal parts of phantom
    NonMetalPart=NonMetalPixels.*doublephantom;  % not metal
    
    [RadonMetal,t2]=radon(MetalPart,theta1);

    [fanbeammetal, F2pos, F2angles] = fanbeam(MetalPart,D,'FanSensorSpacing',0.2);
            
    figure,imshow(MetalPart),title('Metal Part of Phantom');
    
    
    maxfanbeamphantom=max(max(fanbeamphantom));
    
    corruptedSinogram = fanbeamphantom;
    replaceMetalWithValue = 0.5*maxfanbeamphantom;
    
    corruptedSinogram(fanbeammetal ~= 0) =  replaceMetalWithValue;
    
    
    %Image generate artifacts
    
    artifacts = ifanbeam(corruptedSinogram, D, 'FanSensorSpacing',0.2);
    
    figure,imshow(artifacts),title(' metal artifact reconstruction');
    
    
