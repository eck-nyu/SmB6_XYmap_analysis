

function kax = theta2kx( thax, azi0, tilt0, pol0, hv, WF, EB, innerPot )
    % Enter all angles in degrees, all energies in eV
    
    ec = 1.6e-19; me=9.109e-31; hbar=1.0546e-34; 

    eul = [azi0, tilt0, pol0]*pi/180; 
    rotm = eul2rotm(eul,'ZYX');

    kMap = zeros(length(thax),3); % Starting vecs for kx,ky,kz before adding rotations. 
    kMap(:,1) = sin(thax*pi/180); % Init kx is sin of thetax
    kMap(:,2) = zeros(length(thax),1); % Init ky is 0
    kMap(:,3) = cos(thax*pi/180); % Init kz is cos of thetax
    kMap = (rotm*kMap'); % Apply the total rotation matrix to kMap. End up with new [kx,ky,kz] values, only need kx&kz. 

    kMap(1:2,:)=kMap(1:2,:)*(sqrt(((hv-WF-EB)*ec)*2*me)/hbar)*1e-10;
    kMap(3,:)=(sqrt((((hv-WF-EB)*kMap(3,:).*kMap(3,:)+innerPot)*ec)*2*me)/hbar)*1e-10;
    
    kax = kMap(1,:);
end
