function Bilayer_Dispersion_Bandstructure_Matt
clear
close
m = 1;
n = 1;
points = 100;
bound = 1.5e9;
kx = linspace(-bound,bound,points);
ky = linspace(-bound,bound,points);
d = zeros(12,points,points);
for i = 1:points
    for j = 1:points
        H = bilayermatrix(kx(i),ky(j));
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
hold on
%surf(energy1)
%surf(energy2)
%surf(energy3)
%surf(energy4)
surf(energy5)
surf(energy6)
surf(energy7)
surf(energy8)
%surf(energy9)
%surf(energy10)
%surf(energy11)
%surf(energy12)
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
statesdist = numstates1+numstates2+numstates3+numstates4+numstates5+numstates6+numstates7+numstates8+numstates9+numstates10+numstates11+numstates12;
figure; plot(energy,statesdist,'black')
xlabel('Energy (meV)')
ylabel('DOS (a.u.)')
Band_1 = zeros(1,points);Band_2 = zeros(1,points);Band_3 = zeros(1,points);
Band_4 = zeros(1,points);Band_5 = zeros(1,points);Band_6 = zeros(1,points);
Band_7 = zeros(1,points);Band_8 = zeros(1,points);Band_9 = zeros(1,points);
Band_10 = zeros(1,points);Band_11 = zeros(1,points);Band_12 = zeros(1,points);
for i = 1:points
    Band_1(1,i) = d(1,points/2,i);Band_2(1,i) = d(2,points/2,i);Band_3(1,i) = d(3,points/2,i);
    Band_4(1,i) = d(4,points/2,i);Band_5(1,i) = d(5,points/2,i);Band_6(1,i) = d(6,points/2,i);
    Band_7(1,i) = d(7,points/2,i);Band_8(1,i) = d(8,points/2,i);Band_9(1,i) = d(9,points/2,i);
    Band_10(1,i) = d(10,points/2,i);Band_11(1,i) = d(11,points/2,i);Band_12(1,i) = d(12,points/2,i);
end
figure; 
hold on
%plot(kx,1000*Band_1)
%plot(kx,1000*Band_2)
%plot(kx,1000*Band_3)
%plot(kx,1000*Band_4)
plot(kx,1000*Band_5)
plot(kx,1000*Band_6)
plot(kx,1000*Band_7)
plot(kx,1000*Band_8)
%plot(kx,1000*Band_9)
%plot(kx,1000*Band_10)
%plot(kx,1000*Band_11)
%plot(kx,1000*Band_12)
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


function H = bilayermatrix(kx,ky)
phi = 2*pi/180;
a = 2.46*1E-10;
vf = (1.054*1E-34*1.1*1E6)/(1.6*1E-19);
tbar = 0.11;
lambda = a/(2*sin(phi/2));
theta = atan(sin(phi)/(1-cos(phi)));
G1x = (4*pi/(sqrt(3)*lambda))*cos(theta);
G1y = (4*pi/(sqrt(3)*lambda))*sin(theta);
G2x = cos(120*pi/180)*G1x - sin(120*pi/180)*G1y;
G2y = sin(120*pi/180)*G1x + cos(120*pi/180)*G1y;
phase = exp(2*pi*1i/3);
phaseC = conj(phase);
DKx = 4*pi*(cos(phi)-1)/(3*a);
DKy = 4*pi*sin(phi)/(3*a);
k1 = (kx+DKx/2)+1i*(ky+DKy/2);
k2 = (kx-DKx/2)+1i*(ky-DKy/2);
k1_G1 = kx-G1x+DKx/2+1i*(ky-G1y+DKy/2);
k1_G1G2 = kx-G1x-G2x+DKx/2-1i*(ky-G1y-G2y+DKy/2);
k2_G1 = kx+G1x-DKx/2-1i*(ky+G1y-DKy/2);
k2_G1G2 = kx+G1x+G2x-DKx/2-1i*(ky+G1y+G2y-DKy/2);
H = [
    0             vf*conj(k1)   0           0              0           0                tbar         tbar         tbar*phase tbar*phaseC   tbar*phaseC tbar*phase;
    vf*k1         0             0           0              0           0                tbar         tbar         tbar       tbar*phase    tbar        tbar*phaseC;
      
    0             0             0           vf*conj(k1_G1) 0           0                tbar*phase  tbar*phaseC   0          0              0          0;
    0             0             vf*k1_G1    0              0           0                tbar        tbar*phase    0          0              0          0;
 
    0             0             0           0              0           vf*conj(k1_G1G2) tbar*phaseC  tbar*phase   0          0              0          0;                       
    0             0             0           0              vf*k1_G1G2  0                tbar         tbar*phaseC  0          0              0          0;                       
    
    tbar          tbar          tbar*phaseC tbar           tbar*phase  tbar             0             vf*conj(k2) 0          0              0          0;
    tbar          tbar          tbar*phase  tbar*phaseC    tbar*phaseC tbar*phase       vf*k2         0           0          0              0          0;
    
    tbar*phaseC   tbar          0           0              0           0                0             0           0          vf*conj(k2_G1) 0          0;
    tbar*phase    tbar*phaseC   0           0              0           0                0             0           vf*k2_G1   0              0          0;
    
    tbar*phase    tbar          0           0              0           0                0             0           0          0              0          vf*conj(k2_G1G2) ;
    tbar*phaseC   tbar*phase    0           0              0           0                0             0           0          0              vf*k2_G1G2 0];

%H = [
%    0             vf*conj(k1)   0           0              0           0                tbar*phase   tbar*phaseC  tbar        tbar           tbar*phaseC tbar*phase;
%    vf*k1         0             0           0              0           0                tbar         tbar*phase   tbar        tbar           tbar        tbar*phaseC;
%      
%    0             0             0           vf*conj(k1_G1) 0           0                tbar         tbar         tbar*phase  tbar*phaseC    tbar*phaseC tbar*phase;
%    0             0             vf*k1_G1    0              0           0                tbar         tbar         tbar        tbar*phase     tbar        tbar*phaseC;
% 
%    0             0             0           0              0           vf*conj(k1_G1G2) tbar*phaseC  tbar*phase   tbar*phaseC tbar*phase     0           0;                       
%    0             0             0           0              vf*k1_G1G2  0                tbar         tbar*phaseC  tbar        tbar*phaseC    0           0;                       
%    
%    tbar*phaseC   tbar          tbar        tbar           tbar*phase  tbar             0             vf*conj(k2) 0           0              0           0;
%    tbar*phase    tbar*phaseC   tbar        tbar           tbar*phaseC tbar*phase       vf*k2         0           0           0              0           0;
%    
%    tbar          tbar          tbar*phaseC tbar           tbar*phase  tbar             0             0           0           vf*conj(k2_G1) 0           0;
%    tbar          tbar          tbar*phase  tbar*phaseC    tbar*phaseC tbar*phase       0             0           vf*k2_G1    0              0           0;
%    
%    tbar*phase    tbar          tbar*phase  tbar           0           0                0             0           0           0              0           vf*conj(k2_G1G2) ;
%    tbar*phaseC   tbar*phase    tbar*phaseC tbar*phase     0           0                0             0           0           0              vf*k2_G1G2  0];

end