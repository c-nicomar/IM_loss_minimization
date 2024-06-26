%% IM model parameters

IM_model_param.Rs=0.0074; %Stator resistance ohms
IM_model_param.Lds=3.8503e-4; %Stator leakage inductance in H
IM_model_param.Rr=0.0084; %Rotor resistance ohms
IM_model_param.Ldr=3.8503e-4; %Rotor leakage inductance in H
IM_model_param.Lm1=0.0128; %Magnetizing inductance for paralell Rfe model

IM_model_param.p=2; %Pole pairs
IM_model_param.Vdc=400*sqrt(2); %DC voltage at inverter

%IM_model_param.coeff_low_Rfe=7.1e-5; % Coefficients for Rfe=fn(f)in ohms for f<50 Hz: coeff_low_Rfe(1)+coeff_low_Rfe(2)*f+coeff_low_Rfe(3)*f^2 (paralel model)
%IM_model_param.coeff_high_Rfe=2.2e-7; % Coefficients for Rfe=fn(f) in ohms for f>50 Hz: coeff_high_Rfe(1)+coeff_high_Rfe/f (paralel model)

IM_model_param.Rfe_coeff=[1.79665e-3, 17.94871e-6, 2.88913e-8]; %Rfe=(1.79665e-3 + 17.94871e-6.*abs(ws)+2.88913e-8.*ws.^2);

IM_model_param.maxRotorSpeed=pi/30*4000; % Maximum speed that the IM is expected to function with (rad/s)
IM_model_param.currentLimit=600; %Maximum current that the IM can withstand (A)
IM_model_param.maxTorque=900;