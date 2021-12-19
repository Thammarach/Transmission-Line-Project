% EIE/ENE 450 Applied Communications and Transmission Lines
% Instructor: Watcharapan Suwansantisuk

clear all;
RS = 45;        % source impedance (ohm)
f0 = 5e9;    % center frequency (Hz)
w0 = 2*pi*f0;   % center frequency (rad/sec)
fbw = 0.16;     % fractional bandwidth = (w2- w1)/w0
g0 = 1;
g = [0.6180 1.6180 2.0000 1.6180 0.6180 1.0000];  % table of prototype elements
f = linspace(0,10*10^9,3000);  % frequency to plot (Hz)

% transformation (the symbols here are the the symbols in the solution
% with prime signs)
diff = w0*fbw;  % w2-w1 

R0 = g0 / RS

L1 = g(1) * RS * diff / ( w0 ^ (2) )
C1 = 1 / ( g(1) * RS * diff )

L2 = RS / ( g(2) * diff )
C2 = g(2) * diff / ( RS * w0 ^ (2) )

L3 = g(3) * RS * diff / ( w0 ^ (2) )
C3 = 1 / ( g(3) * RS * diff )

L4 = RS / ( g(4) * diff )
C4 = g(4) * diff / ( RS * w0 ^ (2) )

L5 = g(5) * RS * diff / ( w0 ^ (2) )
C5 = 1 / ( g(5) * RS * diff )


GL = g(6)/ RS
RL = 1 / GL

% find input impedance
w = 2*pi*f;  % the angular frequency (rad/sec) to plot

% input impedance looking into C5
Zin5 = parallel( j * w * L5 , 1 ./ ( j * w * C5 ) ) + RL ;

% input impedance looking into C4
Zin4 = parallel( ( j * w * L4 ) + ( 1 ./ ( j * w * C4 ) ) , Zin5 );

% input impedance looking into C2
Zin3 = parallel( j * w * L3 , 1 ./ ( j * w * C3 ) ) + Zin4 ;

% input impedance looking into C3
Zin2 = parallel( ( j * w * L2 ) + ( 1 ./ ( j * w * C2 ) ) , Zin3 );

% input impedance 
Zin = parallel( j * w * L1 , 1 ./ ( j * w * C1 ) ) + Zin2 ;

% alternative code for P_LR:
%   Gam = (zin-RS)./(zin+RS); % reflection coefficient
%   P_LR = 1 ./ ( 1 - abs(Gam).^2 ); % power loss ratio
P_LR = abs( Zin + RS ).^2 ./ ( 2 * RS * ( conj( Zin ) + Zin) ); % power loss ratio
IL = 10 * log10( P_LR ); % insertion loss

subplot (2,1,2);
semilogx(f/10^9, IL, 'Linewidth', 2);
xlabel('frequency (dB)');
ylabel('Attenuation, i.e., IL (dB)');
title('Frequency response (dB) of a maximally flat Bandstop filter (N=5)');
grid on;

subplot (2,1,1);
plot(f/10^9, IL, 'Linewidth', 2);
xlabel('frequency (Hz)');
ylabel('Attenuation, i.e., IL (dB)');
title('Frequency response (Hz) of a maximally flat Bandstop filter (N=5)');
grid on;

% return the equivalent impedance of two impedances in parallel
% Input:
%  Z1, Z2 - two vectors of the same size
% Output:
%  Z      - a vector the same size as Z1, where 
%           Z(i,j) = the parellel impedanec of Z1(i,j) and Z2(i,j)
function Z = parallel(Z1, Z2)
    Z = Z1 .* Z2 ./ (Z1+Z2);
end