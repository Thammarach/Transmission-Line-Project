% Magnitude of the reflection coefficient v.s. frequency (Figure 5.3 in
% the textbook), page 233
% EIE/ENE 450 Applied Communications and Transmission Lines
% Instructor: Watcharapan Suwansantisuk

%-----------------------%
% adjustable parameters %
%-----------------------%
clear all;
f0 = 2 * 10^(9);

B1 = 0.01309307341;
X1 = 42.91287847;

B2 = -0.01309307341;
X2 = -2.912878468;

% a series of RC load
R_load = 35; % (Ohm) resistor at the load
C_load = 1./(2*pi*f0*20); % (F) capacitor in the load

% circuit elements in the matching network
L1 = X1./(2*pi*f0);  % (H) inductor at the matching network (solution 1)
C1 = B1./(2*pi*f0); % (F) capactor at the matching network (solution 1)

C2 = -1./(2*pi*f0*X2); % (F) capactor at the matching network (solution 2)
L2 = -1./(2*pi*f0*B2);  % (H) inductor at the matching network (solution 2)

Z0 = 50; % (Ohm) characteristic impedance
f = linspace(0, 6 * 10^(9));  % (Hz) range of frequencies to plot

%-----------------------%
% program starts here   %
%-----------------------%

ZL = R_load + 1./(j*2*pi*f*C_load);  % load impedance 

% Solution 1
ZC1 = 1./(j*2*pi*f*C1);  % impedance of C in the matching network 
ZL1 = j*2*pi*f*L1;
Zin1 = ( ( ZL + ZL1 ) .* ( ZC1 ) ) ./ ( ZL + ZL1 + ZC1 );  % input impedance
Gamma1 = ( Zin1 - Z0 ) ./ ( Zin1 + Z0 ); % reflection coefficient at the matching network

% Solution 2
ZL2 = j*2*pi*f*L2;
ZC2 = 1./(j*2*pi*f*C2);
Zin2 = ( ( ZL + ZC2 ) .* ( ZL2 ) ) ./ ( ZL + ZL2 + ZC2 );  % impedance of L in the matching network
Gamma2 = ( Zin2 - Z0 ) ./ ( Zin2 + Z0 ); % reflection coefficient at the matching network

plot( f/10^9, abs(Gamma1), ...
      f/10^9, abs(Gamma2), 'Linewidth', 2 );


xlabel('frequency (GHz)');
ylabel('|\Gamma|');
legend('Solution 1', 'Solution 2' );
grid on