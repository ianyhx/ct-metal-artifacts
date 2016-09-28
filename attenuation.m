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