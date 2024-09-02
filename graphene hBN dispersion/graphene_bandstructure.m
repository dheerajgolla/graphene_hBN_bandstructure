t = 2.7; %eV
tdash = -1*0.2*2.8; %eV
a = 1.42; %angstorm
kx = -2:0.05:2;ky = -2:0.05:2;
[Kx,Ky] = meshgrid(kx,ky);
f =  t*sqrt(3+(2*cos(sqrt(3)*Ky*a)+4*cos((sqrt(3)/2)*Ky*a).*cos((3/2)*Kx*a)))- tdash*(2*cos(sqrt(3)*Ky*a)+4*cos((sqrt(3)/2)*Ky*a).*cos((3/2)*Kx*a));
fm = -1*f - tdash*(2*cos(sqrt(3)*Ky*a)+4*cos((sqrt(3)/2)*Ky*a).*cos((3/2)*Kx*a));
figure;surf(Kx,Ky,(f));hold on;surf(Kx,Ky,(fm))


