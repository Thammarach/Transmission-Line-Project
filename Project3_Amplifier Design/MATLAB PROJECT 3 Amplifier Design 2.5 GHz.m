% Single-stage amplifier design for maximum gain, at 2.5 GHz
% EIE/ENE 450 Applied Communications and Transmission Lines
function amplifier_design_for_class
%-----------------------%
% adjustable parameters %
%-----------------------%
clear all;
Z0 = 50; % line impedance (Ohms)
% scattering parameter at 2.5GHz
f = 2.5e9; % frequency (Hz)
S_abs = [0.314 0.482; 1.109 0.423]; % magnitudes of the scattering parameters
S_deg = [57.3 28; 11.7 -135.9]; % angles of the scattering parameters (degrees)
fprintf('S_abs = %d', S_abs(2,1));
% scattering parameters at other frequencies, to check the stability
% and gain
check{1}.f = 1e9; % frequency to check (Hz), i.e., 1 GHz
check{1}.S_abs = [0.352 0.239; 2.254 0.210]; % magnitudes of the scattering parameters at 1 GHz
check{1}.S_deg = [132.4 59.2; 54.4 -124.2]; % angles of the scattering parameters (degrees) at 1 GHz
%---------------------------%
% amplifier design at 2.5 GHz %
%---------------------------%
S = S_abs .* exp( j*pi/180*S_deg );

%---- Part 1: check stability ------
Delta = det(S);
abs_Delta = abs( Delta ); % magnitude of Delta
K = (1 - ( S_abs(1,1) )^2 - ( S_abs(2,2) )^2 + ...
 abs_Delta^2 ) / ( 2*abs( S(1,2)*S(2,1) ) ); % paramter K
% print out
fprintf('\n--Step 1: K-Delta stability test\n');
fprintf('Delta: |Delta| = %f, deg of Delta = %f\n', abs(Delta), angle(Delta)*180/pi );
fprintf('K = %g\n', K );
% print out the result
if ( K > 1 && abs_Delta < 1 )
  fprintf('Unconditionally stable\n');
else
  fprintf('Conditionally stable\n');
end

%---- Part 2: GammaS ------
B1 = 1 + (S_abs(1,1))^2 - (S_abs(2,2))^2 - abs_Delta^2;
C1 = S(1,1) - Delta*conj(S(2,2));
GammaS = ( B1 + [1 -1]*sqrt( B1^2 - 4*(abs(C1))^2) ) / (2*C1);
% print out answer
fprintf('\n--Step 2: GammaS\n');
for g = GammaS
 fprintf('GammaS: |GammaS| = %f, deg of GammaS = %f\n', abs(g), angle(g)*180/pi );
end
fprintf('choose the 2nd solution because it is inside the input stability region\n');
gammaS = GammaS(2);
fprintf('gammaS = %f', gammaS);

%---- Part 3: Zs ------
Zs = Z0*( 1 + gammaS ) / ( 1 - gammaS );
fprintf('\n--Step 3: Zs\n');
fprintf('Zs = %f + j(%f)\n', real(Zs), imag(Zs) );
fprintf('|Zs| = %f, deg of Zs = %f\n', abs(Zs), angle(Zs)*180/pi );

%---- Part 4: input matching circuit ------
tuner_in = single_shunt_stub( conj(Zs), Z0 );
% printout only open-circuited shunt stub
nsol = length(tuner_in);
fprintf('\n--Step 4: input matching circuit\n');
fprintf(1, '[Single-stub shunt tuner] %d solution(s):', nsol );
for k=1:nsol
 fprintf(1, '\nSolution #%d\n', k );
 fprintf(1, ' Distance of the stub: d/lambda = %g\n', tuner_in{k}.d );
 fprintf(1, ' Short circuit: ls/lambda = %g\n', tuner_in{k}.ls );
 fprintf(1, ' Open circuit: lo/lambda = %g\n', tuner_in{k}.lo );
end

%---- Part 5: GammaL ------
B2 = 1 + (S_abs(2,2))^2 - (S_abs(1,1))^2 - abs_Delta^2;
C2 = S(2,2) - Delta*conj(S(1,1));
GammaL = ( B2 + [1 -1]*sqrt( B2^2 - 4*(abs(C2))^2) ) / (2*C2);
% print out answer
fprintf('\n--Step 5: GammaL\n');
for g=GammaL
 fprintf('GammaL: |GammaL| = %f, deg of GammaL = %f\n', abs(g), angle(g)*180/pi );
end
fprintf('choose the 2nd solution because it is inside the input stability region\n');
gammaL = GammaL(2);

%---- Part 6: Zl ------
Zl = Z0*( 1 + gammaL ) / ( 1 - gammaL );
fprintf('\n--Step 6: Zl\n');
fprintf('Zl = %f + j(%f)\n', real(Zl), imag(Zl) );
fprintf('|Zl| = %f, deg of Zl = %f\n', abs(Zl), angle(Zl)*180/pi );

%---- Part 7: output matching circuit ------
tuner_out = single_shunt_stub1( conj(Zl), Z0 );
% printout only open-circuited shunt stub
nsol = length(tuner_out);
fprintf('\n--Step 7: output matching circuit\n');
fprintf(1, '[Single-stub shunt tuner] %d solution(s):', nsol );
for k=1:nsol
 fprintf(1, '\nSolution #%d\n', k );
 fprintf(1, ' Distance of the stub: d/lambda = %g\n', tuner_out{k}.d );
 fprintf(1, ' Short circuit: ls/lambda = %g\n', tuner_out{k}.ls );
 fprintf(1, ' Open circuit: lo/lambda = %g\n', tuner_out{k}.lo );
end

%---- Part 8: Check stability
gamma_out = S(2,2) + S(1,2)*S(2,1)*gammaS/( 1-S(1,1)*gammaS );
gamma_in = S(1,1) + S(1,2)*S(2,1)*gammaL/( 1-S(2,2)*gammaL );
fprintf('\n--Step 8: Careful check of stability\n');
fprintf('The input port of transistor is...');
if ( abs(gammaS) < 1 && abs(gamma_out) < 1 )
 fprintf('stable\n');
else
 fprintf('unstable\n');
end
fprintf('The output port of transistor is...');
if ( abs(gammaL) < 1 && abs(gamma_in) < 1 )
 fprintf('stable\n');
else
 fprintf('unstable\n');
end

%---- Part 9: Compute the tranducer gain
Gs = 1/(1-(abs(gammaS))^2);
G0 = (S_abs(2,1))*(S_abs(2,1));
GL = (1-(abs(gammaL))^2) / ( abs(1-S(2,2)*gammaL) )^2;
GTmax = Gs * G0 * GL;
fprintf('\n--Step 9: Gains\n');
fprintf('Gs = %f = %f dB\n', Gs, 10*log10(Gs) );
fprintf('G0 = %f = %f dB\n', G0, 10*log10(G0) );
fprintf('GL = %f = %f dB\n', GL, 10*log10(GL) );
fprintf('GTmax = %f = %f dB\n', GTmax, 10*log10(GTmax) );
%-----------------------------------------%
% gain and stability at other frequencies %
%-----------------------------------------%
ncheck = length(check); % number of frequencies to check
for k=1:ncheck

 fprintf('\n\nFrequency = %g Hz\n', check{k}.f );

 % scattering parameter at frequency check{k}.f
 SS = check{k}.S_abs .* exp( j*pi/180*check{k}.S_deg );

 % stability at the input port of transistor
 % input impedance seen looking into the input matching network
 % {2} for the second solution of the single shunt stub
 Z = Zin_shunt_stub( f, tuner_in{2}.d, tuner_in{2}.lo, Z0, Z0, 'Short', check{k}.f );
 gamS = (Z - Z0)/(Z + Z0); % reflection coeffient at the source side
 gam_out = SS(2,2) + SS(1,2)*SS(2,1)*gamS/( 1-SS(1,1)*gamS );
 % stability at the ouput port of transistor
 % {2} for the second solution of the single shunt stub
 Z = Zin_shunt_stub( f, tuner_out{2}.d, tuner_out{2}.lo, Z0, Z0, 'Open', check{k}.f );
 gamL = (Z - Z0)/(Z + Z0); % reflection coeffient at the load side
 gam_in = SS(1,1) + SS(1,2)*SS(2,1)*gamL/( 1-SS(2,2)*gamL );

 fprintf(' gamS: |gamS| = %f, deg of gamS = %f\n', abs(gamS), angle(gamS)*180/pi );
 fprintf(' gam_out: |gam_out| = %f, deg of gam_out = %f\n', abs(gam_out), angle(gam_out)*180/pi );
 fprintf(' gamL: |gamL| = %f, deg of gamL = %f\n', abs(gamL), angle(gamL)*180/pi );
 fprintf(' gam_in: |gam_in| = %f, deg of gam_in = %f\n', abs(gam_in), angle(gam_in)*180/pi );
 fprintf('\n The input port of transistor is...');
 if ( abs(gamS) < 1 && abs(gam_out) < 1 )
 fprintf('stable\n');
 else
 fprintf('unstable\n');
 end
 fprintf(' The output port of transistor is...');
 if ( abs(gamL) < 1 && abs(gam_in) < 1 )
 fprintf('stable\n');
 else
 fprintf('unstable\n');
 end

 %---- Compute the tranducer gain
 % Note: This calculation comes from equation (12.13) of the textbook
 % Equation (12.37) will not work here because (12.37) is only for
 % the maximum gain (conjugate matching), at 2.5 GHz
 Gs = (1-(abs(gamS))^2) / ( abs(1-gamS*gam_in) )^2;
 G0 = (check{k}.S_abs(2,1))^2;
 GL = (1-(abs(gamL))^2) / ( abs(1-SS(2,2)*gamL) )^2;
 GTmax = Gs * G0 * GL;
 fprintf('\n Gains:\n');
 fprintf(' Gs = %f = %f dB\n', Gs, 10*log10(Gs) );
 fprintf(' G0 = %f = %f dB\n', G0, 10*log10(G0) );
 fprintf(' GL = %f = %f dB\n', GL, 10*log10(GL) );
 fprintf(' GT = %f = %f dB\n', GTmax, 10*log10(GTmax) );
end
end
% Matches a load Zload to the line with impedance Z0 using a
% shunt-stub tuner
% Input:
% ZL - a complex number for the load impedance (Ohm)
% Z0 - a real number for the characteristic impednace of the
% transmission line (Ohm)
% Output:
% tuner - a cell, where
% tuner{k}.d is (the length from the line to the load)/lambda
% tuner{k}.ls is (the length of the short-circuited stub)/lambda
% tuner{k}.lo is (the length of the open-circuited stub)/lambda
% Here, length(tuner) is the number of solutions
function tuner = single_shunt_stub( ZL, Z0 )
 RL = real(ZL); % load resistance
 XL = imag(ZL); % load reactance
 Y0 = 1/Z0; % characteristic admittance of the line
 % obtain t
 if ( RL ~= Z0 )
 % two solution. Put them in a vector 't'
 t = ( XL + (-1).^[0 1] * sqrt( RL*( (Z0-RL)^2 + XL^2 )/Z0 ) ) ...
 / ( RL - Z0 );
 else
 % one solution
 t = -XL/(2*Z0);
 end
 % obtain B
 B = ( RL^2*t - (Z0-XL*t).*(XL + Z0*t) ) ...
 ./ ( Z0*(RL^2 + (XL + Z0*t).^2 ) );
 % obtain normalized length, norm_ls = ls/lambda, of the short-circuit stub
 norm_ls = atan( Y0./B ) / (2*pi);
 norm_ls( norm_ls < 0 ) = norm_ls( norm_ls < 0 ) + 1/2;
 % obtain normalized length, norm_lo = lo/lambda, of the open-circuit stub
 norm_lo = -atan( B/Y0 ) / (2*pi);
 norm_lo( norm_lo < 0 ) = norm_lo( norm_lo < 0 ) + 1/2;
 % obtain the normalized distance, norm_d = d/lambda, of the stub
 % Note that t can be a vector, so we can have multiple solutions of norm_d
 norm_d = atan( t ) / (2*pi);
 norm_d( t<0 ) = norm_d( t< 0 ) + 1/2;
 % prepare the output
 nsol = length( norm_d ); % number of solution
 tuner = cell( nsol );
 for k=1:nsol
 tuner{k}.d = norm_d(k);
 tuner{k}.ls = norm_ls(k);
 tuner{k}.lo = norm_lo(k);
 end
 %---------------------%
% Gamma vs frequency %
%---------------------%
f0 = 2.5*10^9; % (Hz) frequency for which the load is ZL as given
f = linspace( 1.5*10^9, 4*10^9, 500 ); % (Hz) range of frequencies to plot
Z0 = 50; % (Ohm) characteristic impedance
% Suppose ZL is a series of RC or of RL
if ( XL < 0 ) % series of R and C
 C = 1/(-XL*2*pi*f0); % capacitance (F) in the series
 ZL_f = RL + 1./(j*2*pi*f*C); % load impedance, as a function of frequency
 % show the components in the load
else
 L = XL / (2*pi*f0); % indunctance (H) in the series
 ZL_f = RL + j*2*pi*f*L; % load impedance, as a function of frequency
 % show the components in the load
end
% compute |Gamma| vs frequency
figure(1);
clf;
hold all;
leg = cell( 1, 2*nsol ); % legend
for k=1:nsol
 % impedance down a length d from the load, equation (5.7) in the
 % textbook
 tan_term = tan( 2*pi*norm_d(k)*f/f0 ); % tan (beta*d)
 Z = Z0*( ZL_f + j*Z0*tan_term ) ./ (Z0 + j*ZL_f.*tan_term);

 % impedances of the stub transmission line
 % circuit
 Z_stub_o = -j*Z0*cot( 2*pi*norm_lo(k)*f/f0 ); % open stub
 Z_stub_s = j*Z0*tan( 2*pi*norm_ls(k)*f/f0 ); % short-circuited stub

 % the over-all input impedance looking into the matching network:
 % A parellel of Z and Z_sub
 Zin_o = Z .* Z_stub_o ./ (Z + Z_stub_o); % open stub
 Zin_s = Z .* Z_stub_s ./ (Z + Z_stub_s); % short-circuited stub

 % Gamma
 Gamma_o = (Zin_o - Z0) ./ (Zin_o + Z0); % open stub
 Gamma_s = (Zin_s - Z0) ./ (Zin_s + Z0); % short-circuited stub

 % plot and set legend
 figure(1);
 plot( f/10^9, [ abs(Gamma_o) ; abs(Gamma_s) ], 'Linewidth', 2 );
 leg{2*k-1} = sprintf('Sol #%d: open', k );
 leg{2*k} = sprintf('Sol #%d: shorted', k );
 grid on;
end
legend( leg );
title('Magnitude of the reflection coefficient at various frequencies');
xlabel('f (GHz)');
ylabel('|\Gamma|');
end
function tuner = single_shunt_stub1( ZL, Z0 )
 RL = real(ZL); % load resistance
 XL = imag(ZL); % load reactance
 Y0 = 1/Z0; % characteristic admittance of the line
 % obtain t
 if ( RL ~= Z0 )
 % two solution. Put them in a vector 't'
 t = ( XL + (-1).^[0 1] * sqrt( RL*( (Z0-RL)^2 + XL^2 )/Z0 ) ) ...
 / ( RL - Z0 );
 else
 % one solution
 t = -XL/(2*Z0);
 end
 % obtain B
 B = ( RL^2*t - (Z0-XL*t).*(XL + Z0*t) ) ...
 ./ ( Z0*(RL^2 + (XL + Z0*t).^2 ) );
 % obtain normalized length, norm_ls = ls/lambda, of the short-circuit stub
 norm_ls = atan( Y0./B ) / (2*pi);
 norm_ls( norm_ls < 0 ) = norm_ls( norm_ls < 0 ) + 1/2;
 % obtain normalized length, norm_lo = lo/lambda, of the open-circuit stub
 norm_lo = -atan( B/Y0 ) / (2*pi);
 norm_lo( norm_lo < 0 ) = norm_lo( norm_lo < 0 ) + 1/2;
 % obtain the normalized distance, norm_d = d/lambda, of the stub
 % Note that t can be a vector, so we can have multiple solutions of norm_d
 norm_d = atan( t ) / (2*pi);
 norm_d( t<0 ) = norm_d( t< 0 ) + 1/2;
 % prepare the output
 nsol = length( norm_d ); % number of solution
 tuner = cell( nsol );
 for k=1:nsol
 tuner{k}.d = norm_d(k);
 tuner{k}.ls = norm_ls(k);
 tuner{k}.lo = norm_lo(k);
 end
 %---------------------%
% Gamma vs frequency %
%---------------------%
f0 = 2.5*10^9; % (Hz) frequency for which the load is ZL as given
f = linspace( 1.5*10^9, 4*10^9, 500 ); % (Hz) range of frequencies to plot
Z0 = 50; % (Ohm) characteristic impedance
% Suppose ZL is a series of RC or of RL
if ( XL < 0 ) % series of R and C
 C = 1/(-XL*2*pi*f0); % capacitance (F) in the series
 ZL_f = RL + 1./(j*2*pi*f*C); % load impedance, as a function of frequency
 % show the components in the load
else
 L = XL / (2*pi*f0); % indunctance (H) in the series
 ZL_f = RL + j*2*pi*f*L; % load impedance, as a function of frequency
 % show the components in the load
end
% compute |Gamma| vs frequency
figure(2);
clf;
hold all;
leg = cell( 1, 2*nsol ); % legend
for k=1:nsol
 % impedance down a length d from the load, equation (5.7) in the
 % textbook
 tan_term = tan( 2*pi*norm_d(k)*f/f0 ); % tan (beta*d)
 Z = Z0*( ZL_f + j*Z0*tan_term ) ./ (Z0 + j*ZL_f.*tan_term);

 % impedances of the stub transmission line
 % circuit
 Z_stub_o = -j*Z0*cot( 2*pi*norm_lo(k)*f/f0 ); % open stub
 Z_stub_s = j*Z0*tan( 2*pi*norm_ls(k)*f/f0 ); % short-circuited stub

 % the over-all input impedance looking into the matching network:
 % A parellel of Z and Z_sub
 Zin_o = Z .* Z_stub_o ./ (Z + Z_stub_o); % open stub
 Zin_s = Z .* Z_stub_s ./ (Z + Z_stub_s); % short-circuited stub

 % Gamma
 Gamma_o = (Zin_o - Z0) ./ (Zin_o + Z0); % open stub
 Gamma_s = (Zin_s - Z0) ./ (Zin_s + Z0); % short-circuited stub

 % plot and set legend
 figure(2);
 plot( f/10^9, [ abs(Gamma_o) ; abs(Gamma_s) ], 'Linewidth', 2 );
 leg{2*k-1} = sprintf('Sol #%d: open', k );
 leg{2*k} = sprintf('Sol #%d: shorted', k );
 grid on;
end
legend( leg );
title('Magnitude of the reflection coefficient at various frequencies');
xlabel('f (GHz)');
ylabel('|\Gamma|');
end
% Obtain the input impedance seen looking into a single-stub shunt tuner at
% a given frequency
% Input:
% f - frequency where the shunt tuner is designed for (Hz)
% d - (length of the line)/lambda, where lambda = wavelength on the
% transmission line for frequency f
% l - (length of the stub)/lambda
% Zload - the impedance of the load, a complex number (Ohm)
% Z0 - the impedance of the transmission line,
% a complex number although usually be real (Ohm)
% type - type of the stub, either 'Open' (for open-circuited) or 'Close'
% for (for closed-circuit)
% f_given - frequency to find the input impedance (Hz)
% Output:
% Z - the input impedance seen
function Z = Zin_shunt_stub( f, d, l, Zload, Z0, type, f_given )
 % input impedance of the stub at the given frequency
 beta_l = 2*pi*f_given/f * l; % beta (at the given frequency) * (stub length)
 if ( type(1) == 'O' )
 Zstub = -j*Z0*cot( beta_l );
 else
 Zstub = j*Z0*tan( beta_l );
 end

 % Zload || Zstub
 Zeq = Zload * Zstub / (Zload + Zstub);

 % input impedance of a line of length d terminated by a load Zeq
 beta_d = 2*pi*f_given/f * d; % beta (at the given frequency) * (line length)
 t = tan( beta_d );
 Z = Z0*( Zeq + j*Z0*t)/( Z0 + j*Zeq*t );

 fprintf(' >>> debugging: Zstub = %f + j(%f) Ohm\n', real(Zstub), imag(Zstub) );
 fprintf(' >>> debugging: Zeq = %f + j(%f) Ohm\n', real(Zeq), imag(Zeq) );
 fprintf(' >>> debugging: Z = %f + j(%f) Ohm\n', real(Z), imag(Z) );
 grid on
end
