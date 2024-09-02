function GBN_Dispersion_Bandstructure_one
clear
close all
m = 1;
n = 1;
points = 100;
bound = 0.5;
kx = linspace(-bound,bound,points);
ky = linspace(-bound,bound,points);
d = zeros(14,points,points);
for i = 1:points
    for j = 1:points
        H = mymatrix(kx(i),ky(j));
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
hold off
E_min = -0.4;
E_max = 0.4;
delta_E = 0.01;
bins = (E_max-E_min)/delta_E;
energy = linspace(E_min,E_max,bins);
energy = energy*1000;
numstates1 = numstates(energy1,E_min,E_max,delta_E,points);
numstates2 = numstates(energy2,E_min,E_max,delta_E,points);
numstates3 = numstates(energy3,E_min,E_max,delta_E,points);
numstates4 = numstates(energy4,E_min,E_max,delta_E,points);
numstates5 = numstates(energy5,E_min,E_max,delta_E,points);
numstates6 = numstates(energy6,E_min,E_max,delta_E,points);
numstates7 = numstates(energy7,E_min,E_max,delta_E,points);
numstates8 = numstates(energy8,E_min,E_max,delta_E,points);
numstates9 = numstates(energy9,E_min,E_max,delta_E,points);
numstates10 = numstates(energy10,E_min,E_max,delta_E,points);
numstates11 = numstates(energy11,E_min,E_max,delta_E,points);
numstates12 = numstates(energy12,E_min,E_max,delta_E,points);
numstates13 = numstates(energy13,E_min,E_max,delta_E,points);
numstates14 = numstates(energy14,E_min,E_max,delta_E,points);
%statesdist = numstates1+numstates2+numstates3+numstates4+numstates5+numstates6+numstates7+numstates8+numstates9+numstates10+numstates11+numstates12+numstates13+numstates14;
statesdist = numstates6+numstates7+numstates8+numstates9;
figure; plot(energy,statesdist,'black')
xlabel('Energy (meV)')
ylabel('DOS (a.u.)')
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
figure; 
hold on
%plot(kx,1000*Band_1)
%plot(kx,1000*Band_2)
%plot(kx,1000*Band_3)
%plot(kx,1000*Band_4)
%plot(kx,1000*Band_5)
plot(kx,1000*Band_6)
plot(kx,1000*Band_7)
plot(kx,1000*Band_8)
plot(kx,1000*Band_9)
%plot(kx,1000*Band_10)
%plot(kx,1000*Band_11)
%plot(kx,1000*Band_12)
%plot(kx,1000*Band_13)
%plot(kx,1000*Band_14)
axis([-bound bound 1000*E_min 1000*E_max])
pbaspect([0.7 1 1])
xlabel('k')
ylabel('Energy (meV)')
hold off
% dlmwrite('conductionband_lowres.txt', energy8);
% type conductionband_lowres.txt;
% dlmwrite('valenceband_lowres.txt', energy7);
% type valenceband_lowres.txt;
% dlmwrite('newvalenceband_lowres.txt', energy6);
% type newvalenceband_lowres.txt;
end

function statesdist = numstates(energy_vals,E_min,E_max,delta_E,points)
energy = E_min;
energy_index = 1;
bins = (E_max-E_min)/delta_E;
statesdist = zeros(bins,1);
while energy < E_max
    for i = 1:points
        for j = 1:points
            if energy_vals(i,j) > energy && energy_vals(i,j) <= energy + delta_E
                statesdist(energy_index) = statesdist(energy_index) + 1;
            end
        end
    end
    energy = energy + delta_E;
    energy_index = energy_index + 1;
end
end

function H = mymatrix(kx,ky)
Vg = 0.206202;
G = 0.711;
theta = 0.766; %angle between planes
H = [0 0.687*(kx+G*cos(theta)-1i*(-ky-G*sin(theta))) Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0;
     0.687*(kx+G*cos(theta)+1i*(-ky-G*sin(theta))) 0 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 ;
     Vg/2 0 0 0.687*(kx+G*cos(theta+pi/3)-1i*(-ky-G*sin(theta+pi/3))) Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0;
     0 Vg/2 0.687*(kx+G*cos(theta+pi/3)+1i*(-ky-G*sin(theta+pi/3))) 0 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 0 0.687*(kx+G*cos(theta+2*pi/3)-1i*(-ky-G*sin(theta+2*pi/3))) Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0;
     0 Vg/2 0 Vg/2 0.687*(kx+G*cos(theta+2*pi/3)+1i*(-ky-G*sin(theta+2*pi/3))) 0 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx+G*cos(theta+pi)-1i*(-ky-G*sin(theta+pi))) Vg/2 0 Vg/2 0 Vg/2 0;
     0 Vg/2 0 Vg/2 0 Vg/2 0.687*(kx+G*cos(theta+pi)+1i*(-ky-G*sin(theta+pi))) 0 0 Vg/2 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx+G*cos(theta+4*pi/3)-1i*(-ky-G*sin(theta+4*pi/3))) Vg/2 0 Vg/2 0;
     0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0.687*(kx+G*cos(theta+4*pi/3)+1i*(-ky-G*sin(theta+4*pi/3))) 0 0 Vg/2 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx+G*cos(theta+5*pi/3)-1i*(-ky-G*sin(theta+5*pi/3))) Vg/2 0;
     0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0.687*(kx+G*cos(theta+5*pi/3)+1i*(-ky-G*sin(theta+5*pi/3))) 0 0 Vg/2;
     Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 0 0.687*(kx + 1i*ky);
     0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0 Vg/2 0.687*(kx - 1i*ky) 0];
end
