function GBN_Dispersion_absorption
clear
close all

%c = 3*10^(10); %speed of light gaussian
c=1;
m = 1;
n = 1;
points = 100; %k space mesh
bound = 2.0 ; %bound in k-space in nm^(-1) measure from original Dirac point

Ex = 10^(-8); %electric field of light (Gaussian units)
Ey = Ex;
Wi = (c/4*pi)*(Ex^2 + Ey^2); %incident energy 
%Wi = Wi*6.24*10^(11);%incident energy in eV


kx = linspace(-bound,bound,points);
ky = linspace(-bound,bound,points);
d = zeros(14,points,points);

for i = 1:points
    for j = 1:points
        H = mymatrix(kx(i),ky(j),0,0);
        d(:,m,n) = eig(H);        
        m = m + 1;
    end
    m = 1;
    n = n + 1;
end
energy1(:,:) = d(1,:,:);
energy2(:,:) = d(2,:,:);
energy3(:,:) = d(3,:,:);
energy4(:,:) = d(4,:,:);
energy5(:,:) = d(5,:,:);
energy6(:,:) = d(6,:,:);
energy7(:,:) = d(7,:,:);
energy8(:,:) = d(8,:,:);
energy9(:,:) = d(9,:,:);
energy10(:,:) = d(10,:,:);
energy11(:,:) = d(11,:,:);
energy12(:,:) = d(12,:,:);
energy13(:,:) = d(13,:,:);
energy14(:,:) = d(14,:,:);
figure;
 hold on
%surf(energy1)
%surf(energy2)
%surf(energy3)
%surf(energy4)
%surf(energy5)
surf(energy6)
surf(energy7)
surf(energy8)
surf(energy9)
%surf(energy10)
%surf(energy11)
%surf(energy12)
%surf(energy13)
%surf(energy14)
xlabel('k_x')
ylabel('k_y')
zlabel('Energy (eV)')
shading interp
 hold off
E_min = 0.6; %min and max of photon energy
E_max = 1.5;
delta_E = 0.01;
bins = ceil((E_max-E_min)/delta_E);
energy = linspace(E_min,E_max,bins);
%energy = energy*1000;
% numstates1 = numstates(energy1,E_min,E_max,delta_E,points);
% numstates2 = numstates(energy2,E_min,E_max,delta_E,points);
% numstates3 = numstates(energy3,E_min,E_max,delta_E,points);
% numstates4 = numstates(energy4,E_min,E_max,delta_E,points);
% numstates5 = numstates(energy5,E_min,E_max,delta_E,points);
abs78 = absorption(energy7,energy8,E_min,E_max,delta_E,points,7,8,kx,kx,Ex,Ey,Wi);
abs79 = absorption(energy7,energy9,E_min,E_max,delta_E,points,7,9,kx,ky,Ex,Ey,Wi);
abs68 = absorption(energy6,energy8,E_min,E_max,delta_E,points,6,8,kx,ky,Ex,Ey,Wi);
abs69 = absorption(energy6,energy9,E_min,E_max,delta_E,points,6,9,kx,ky,Ex,Ey,Wi);
% numstates10 = numstates(energy10,E_min,E_max,delta_E,points);
% numstates12 = numstates(energy12,E_min,E_max,delta_E,points);
% numstates13 = numstates(energy13,E_min,E_max,delta_E,points);
% numstates14 = numstates(energy14,E_min,E_max,delta_E,points);
% statesdist = numstates1+numstates2+numstates3+numstates4+numstates5+numstates6+numstates7+numstates8+numstates9+numstates10+numstates11+numstates12+numstates13+numstates14;
abso = abs78+abs79+abs68+abs69;
figure; 
hold on
plot(energy,abso,'r')
xlabel('Energy (eV)')
ylabel('absorption (a.u.)')
Band_1 = zeros(1,points);Band_2 = zeros(1,points);Band_3 = zeros(1,points);
Band_4 = zeros(1,points);Band_5 = zeros(1,points);Band_6 = zeros(1,points);
Band_7 = zeros(1,points);Band_8 = zeros(1,points);Band_9 = zeros(1,points);
Band_10 = zeros(1,points);Band_11 = zeros(1,points);Band_12 = zeros(1,points);
Band_13 = zeros(1,points);Band_14 = zeros(1,points);
for i = 1:points
    Band_1(1,i) = d(1,points/2,i);Band_2(1,i) = d(2,points/2,i);Band_3(1,i) = d(3,points/2,i);
    Band_4(1,i) = d(4,points/2,i);Band_5(1,i) = d(5,points/2,i);Band_6(1,i) = d(6,points/2,i);
    Band_7(1,i) = d(7,points/2,i);Band_8(1,i) = d(8,points/2,i);Band_9(1,i) = d(9,points/2,i);
    Band_10(1,i) = d(10,points/2,i);Band_11(1,i) = d(11,points/2,i);Band_12(1,i) = d(12,points/2,i);
    Band_13(1,i) = d(13,points/2,i);Band_14(1,i) = d(14,points/2,i);
end
% figure; 
% hold on
% plot(kx,1000*Band_1)
% plot(kx,1000*Band_2)
% plot(kx,1000*Band_3)
% plot(kx,1000*Band_4)
% plot(kx,1000*Band_5)
% plot(kx,1000*Band_6)
% plot(kx,1000*Band_7)
% plot(kx,1000*Band_8)
% plot(kx,1000*Band_9)
% plot(kx,1000*Band_10)
% plot(kx,1000*Band_11)
% plot(kx,1000*Band_12)
% plot(kx,1000*Band_13)
% plot(kx,1000*Band_14)
% axis([-bound bound 1000*E_min 1000*E_max])
% pbaspect([0.7 1 1])
% xlabel('k')
% ylabel('Energy (meV)')
% hold off
% dlmwrite('conductionband_lowres.txt', energy8);
% type conductionband_lowres.txt;
% dlmwrite('valenceband_lowres.txt', energy7);
% type valenceband_lowres.txt;
% dlmwrite('newvalenceband_lowres.txt', energy6);
% type newvalenceband_lowres.txt;
end

function abso = absorption(energy_valv, energy_valc,E_min,E_max,delta_E,points,r,t,kx,ky,Ex,Ey,Wi)
difference_energy = energy_valc - energy_valv;
%min(min(difference_energy));
%size(difference_energy);
r;
t;
energy = E_min;
c=1;
%c = 3*10^(10);
energy_index = 1;
bins = ceil((E_max-E_min)/delta_E);
abso = zeros(bins,1);
%hbar = 1;
hbar = 6.582*10^(-16); %hbar in units of eV.s
factor = 2*pi/(hbar); %2*pi/hbar (eV.s)^(-1)
while energy < E_max - delta_E
    omega = energy/hbar; %units of s^(-1) same in CGS and SI units
    Ax = 1i*c*Ex/omega; %Gaussian units
    Ay = 1i*c*Ey/omega;
    for i = 1:points
        for j = 1:points
            if difference_energy(i,j) > energy && difference_energy(i,j) <= energy + delta_E
                Hlight = mymatrix(kx(i),ky(j),Ax,Ay);  % Hamiltonian under light influence (note the Ax/energy and Ay/energy)
                H = mymatrix(kx(i),ky(j),0,0); % Hamiltonian unperturbed by light
                [V,D] = eig(H); %eigenvals and eigen states of unperturbed hamiltonian (already calculated in main function - so can be shortened later)
                Hpert = Hlight-H; %Perturbation 
                Msq = (abs((V(:,1))' * Hpert * V(:,2))).^2; %square of perturbation matrix element between initial and final states
                Wa = factor*(energy*Msq)/delta_E;
                abso(energy_index) = (Wa/Wi) + abso(energy_index) ; %adds to absorption - note the energy.*(Msq) term
            end
        end
    end
    energy = energy + delta_E;
    energy_index = energy_index + 1;
end
end

function H = mymatrix(kx,ky,Ax,Ay) %calculates Hamiltonian (centered around gamma point)
% e = 4.8*10^(-10); %charge of an electron in gaussian units
% c = 3*10^(10);
%vf = 1*10^8;
%hbar = 1.05*10^(-27); %hbar in CGS units

c=1;
e=1;
ap= 0.246; %lattice constant in nm
delta = 0.018; %lattice mismatch
phip = (pi/180)*3.5; %increase to decrease moire wavelength and push the new Dirac points to higher energy
theta = atan((sin(phip))/(1+ delta - cos(phip)));
lambda = ((1+delta)*ap)/(sqrt((2*(1+delta)*(1-cos(phip))) + delta^2)); %moire wavelength
Gp = (4*pi/(sqrt(3)*lambda));
vf = 1.1*10^15; %nm/s
hbar = 6.58*10^(-16); %eV.s
Vg = hbar*Gp*vf/2;%eV
%Vg = 0.206202; 
%theta = 0.766;
%Vg = 0.06;
%G = 0.711;  %change to move dip in DOS
G=Gp;
%G=0;
 % 0.678 is hcross*v0 (the coeff in the normal dispersion relation)
 % The lower right 2*2 matric is the Hamiltonian of graphene linearized about Ddirac point
 H = [0 0.687*(kx - (e*Ax/(hbar*c)) +G*cos(theta)-1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta))) Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0;
     0.687*(kx - (e*Ax/(hbar*c)) + G*cos(theta)+1i*(-(ky - (e*Ay/(hbar*c)))- G*sin(theta))) 0 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 ;
     Vg/2 0 0 0.687*(kx - (e*Ax/(hbar*c)) + G*cos(theta+pi/3)-1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+pi/3))) Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0;
     0 Vg/2 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+pi/3)+1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+pi/3))) 0 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 0 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+2*pi/3)-1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+2*pi/3))) Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0;
     0 Vg/2 0 Vg/2 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+2*pi/3)+1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+2*pi/3))) 0 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+pi)-1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+pi))) Vg/2 0 Vg/2 0 Vg/2 0;
     0 Vg/2 0 Vg/2 0 Vg/2 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+pi)+1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+pi))) 0 0 Vg/2 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+4*pi/3)-1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+4*pi/3))) Vg/2 0 Vg/2 0;
     0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+4*pi/3)+1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+4*pi/3))) 0 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+5*pi/3)-1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+5*pi/3))) Vg/2 0;
     0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0.687*(kx - (e*Ax/(hbar*c))+G*cos(theta+5*pi/3)+1i*(-(ky - (e*Ay/(hbar*c)))-G*sin(theta+5*pi/3))) 0 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx - (e*Ax/(hbar*c)) + 1i*(ky - (e*Ay/(hbar*c))));
     0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0.687*( kx - (e*Ax/(hbar*c)) - 1i*(ky - (e*Ay/(hbar*c)))) 0];
 
end
