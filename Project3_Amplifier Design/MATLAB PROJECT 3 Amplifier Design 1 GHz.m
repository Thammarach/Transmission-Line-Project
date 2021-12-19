clear all;
% S parameter at 1GHz
S_abs = [0.352 0.239; 2.254 0.210]; % magnitudes of the scattering parameters
S_deg = [132.4 59.2; 54.4 -124.2]; % angles of the scattering parameters (degrees)
S = S_abs .* exp( j*pi/180*S_deg );

%input matching network
%find ZA
ls_in = 0.0855975;
d_in = 0.333558;
f_old = 2.5e9;
f_new = 1.0e9;
Zo=50;
Za =j*Zo*tan((2*pi*ls_in*f_new)/f_old); %short ckt
%find Zeq_in
Zeq_in =(Zo*Za)/(Zo+Za);
%find Zs at 2.5GHz
Zs =Zo*((Zeq_in+j*Zo*tan((2*pi*d_in*f_new)/f_old))...
 /((Zo+j*Zeq_in*tan((2*pi*d_in*f_new)/f_old))));
%find gamma_s at 2.5GHz
gamma_s =(Zs-Zo)/(Zs+Zo);
%find gamma_out at 2.5GHz
gamma_out = S(2,2)+((S(2,1)*S(1,2)*gamma_s)/(1-S(1,1)*gamma_s));
fprintf('Za is %f%+fj ohm\n', real(Za), imag(Za));
fprintf('Zeq_in is %f%+fj ohm\n', real(Zeq_in), imag(Zeq_in));
fprintf('Zs is %f%+fj ohm\n', real(Zs), imag(Zs));
fprintf('Gamma_s is %fexp(j%f)\n',abs( gamma_s),(180/pi)*angle( gamma_s));
fprintf('Gamma_out is %fexp(j%f)\n',abs( gamma_out),(180/pi)*angle( gamma_out));

% output matching network
%find ZB
lo_out = 0.177888;
d_out = 0.158936;
Zb =-j*Zo*cot((2*pi*lo_out*f_new)/f_old); %open ckt
%find Zeq_out
Zeq_out =(Zo*Zb)/(Zo+Zb);
%find ZL at 2.5GHz
ZL =Zo*((Zeq_out + j*Zo*tan((2*pi*d_out*f_new)/f_old))...
 /((Zo + j*Zeq_out*tan((2*pi*d_out*f_new)/f_old))));
%find gamma_L at 2.5GHz
gamma_L =(ZL-Zo)/(ZL+Zo);
%find gamma_in at 2.5GHz
gamma_in = S(1,1)+((S(2,1)*S(1,2)*gamma_L)/(1-S(2,2)*gamma_L));
fprintf('Zb is %f%+fj ohm\n', real(Zb), imag(Zb));
fprintf('Zeq_out is %f%+fj ohm\n', real(Zeq_out), imag(Zeq_out));
fprintf('ZL is %f%+fj ohm\n', real(ZL), imag(ZL));
fprintf('Gamma_L is %fexp(j%f)\n',abs( gamma_L),(180/pi)*angle( gamma_L));
fprintf('Gamma_in is %fexp(j%f)\n',abs( gamma_in),(180/pi)*angle( gamma_in));

%find transducer power gain by gamma_L, gamma_in, gamma_s
GT = (1-abs(gamma_s)^2)*abs(S(2,1))*abs(S(2,1))*(1-abs(gamma_L)^2);
GT = GT/(abs(1 - gamma_in*gamma_s)*abs(1 - gamma_in*gamma_s));
GT = GT/(abs(1 - S(2,2)*gamma_L)*abs(1 - S(2,2)*gamma_L));
fprintf('Transducer power gain %f\n', GT);
GT = 10*log10(GT);
fprintf('Transducer power gain in dB is %f dB\n', GT);
