a = 0.246; %lattice vector for graphene in nm
del = 0.017; %lattice mismatch
phi = 0:0.05:30; %angle between lattices
phir = (pi/180) .* (phi); %radians
lambda = ((1+del)*a)./sqrt((2*(1+del).*(1-cosd(phi)))+(del^2));
plot(phi,lambda)
%%%%%%%%%%%%%%%%%%%%
G = 4*pi./(sqrt(3)*lambda);
vf = 1.1*10^15; %nm/s
hbar = 6.58*10^(-16); %eV.s
Edip = hbar*G*vf/2;
plot(phi,Edip)
axis tight
xlabel('angle \phi')
ylabel('Energy of SDP (eV)')
%print -depsc2 anglevsphi.jpg