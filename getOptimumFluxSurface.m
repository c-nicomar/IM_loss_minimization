function  [surface_coefficients,fitresult_poly, st, v, r, gof]=getOptimumFluxSurface(IM_model_param)
tic
%This function provides the parameters p00, p10... for a poly 3-3 surface
%function providing the optimum stator flux reference (in terms of loss
%minimization) in steady state for a squirrel cage induction motor,
%provided reference values for rotor speed (wr in rad/s) and
%electromagnetic torque (Te in Nm).

%Inputs: Inputs for this model are expected to be stored inside a Matlab
%structure and contain the specific values for the IM model:


%Outputs: The output is a vector containing the 10 coefficients for the
%poly 3-3 function. The function will be built as:

%optFluxValue_manual=surface_coefficients(1)+surface_coefficients(2)*x+surface_coefficients(3)*y+surface_coefficients(4)*x^2+surface_coefficients(5)*x*y+surface_coefficients(6)*y^2 ...
%+surface_coefficients(7)*x^3+surface_coefficients(8)*x^2*y+surface_coefficients(9)*x*y^2+surface_coefficients(10)*y^3;

%Where x is the rotor speed and y the electromagnetic torque


%% PART 1: Solving steady state model problem for all possible working points

%Parameter assignment

Rs=IM_model_param.Rs; %Stator resistance ohms
Lds=IM_model_param.Lds; %Stator leakage inductance in H
Rr=IM_model_param.Rr; %Rotor resistance ohms
Ldr=IM_model_param.Ldr; %Rotor leakage inductance in H
Lm1=IM_model_param.Lm1; %Magnetizing inductance for paralell Rfe model

p=IM_model_param.p; %Pole pairs
Vdc=IM_model_param.Vdc; %DC voltage at inverter

% coeff_low_Rfe=IM_model_param.coeff_low_Rfe; % Coefficients for Rfe=fn(f)in ohms for f<50 Hz: coeff_low_Rfe(1)+coeff_low_Rfe(2)*f+coeff_low_Rfe(3)*f^2 (paralel model)
% coeff_high_Rfe=IM_model_param.coeff_high_Rfe; % Coefficients for Rfe=fn(f) in ohms for f>50 Hz: coeff_high_Rfe(1)+coeff_high_Rfe/f (paralel model)

Rfe_coeff=IM_model_param.Rfe_coeff;

maxRotorSpeed=IM_model_param.maxRotorSpeed; % Maximum speed that the IM is expected to function with (rad/s)
currentLimit=IM_model_param.currentLimit; %Maximum current that the IM can withstand (A)
maxTorque=IM_model_param.maxTorque; %Maxium torque considered

% Sweep of input values


%Voltage Vrms in stator
maxVrms=Vdc/sqrt(2);
stepVrms=round(maxVrms/400);
Vrms_in=0:stepVrms:maxVrms;  

%Slip (fixed range)
s_in=-0.05:0.002:0.05;

%Rotor speed (rpm)
max_n=30/pi*maxRotorSpeed;
step_n=round(max_n/3000,1);
n_in=0:step_n:max_n;

%Determine sizes
size_V=length(Vrms_in);
size_s=length(s_in);
size_n=length(n_in);

% Build input vectors for all combinations of 3 inputs
n=repmat(n_in,1,size_V*size_s)';
s=repelem(s_in,size_n*size_V)';
Vrms=repmat(repelem(Vrms_in,size_n),1,size_s)';


%Compute other variables for all cases
wr=n*2*pi/60;    % Rotor speed in rad/s       
we=wr*p;            %Electric speed in rad/s

ws=we./(1-s);       %Syncronous speed in rad/s
f=ws./(2*pi);       %Syncronous frequency in Hz

%Rfe
Rfe=Rfe_coeff(1)+ Rfe_coeff(2).*abs(ws)+Rfe_coeff(3).*ws.^2;
Rfe(f>50)=25.5425e-3-4.82058./ws(f>50);

%Rfe=coeff_low_Rfe.*ws+coeff_high_Rfe.*ws.^2;

Lm=Lm1;

Ls=Lds+Lm;
Lr=Ldr+Lm;


% Prepare input vector for the IM model (voltage vector)

Vls=Vrms*sqrt(2);                       
PhaseV=0;                               
PhaseVr=PhaseV*pi/180;                  
Vs=Vls/sqrt(3)*exp(j*PhaseVr);          
                          

v=nan(length(s),2);
st=nan(length(s),2);
r=nan(length(s),2);

% Compute values for each vector of the model, st containing stator flux
% and current and r containing rotor flux and current

for i=1:length(s)
   
    A=[j*ws(i) Rs+Rfe(i); 0 Rfe(i)];
    B=[0 Rfe(i); j*(ws(i)-we(i)) Rr+Rfe(i)];
%     C=[0 0; 0 0];
%     D=[0 0; 0 0];
    E=[Lr/Lm (Lm^2-Lr*Ls)/Lm; 1/Lm -Ls/Lm];
    
    aux=A+B*E;
    
    v(i,:)=[Vs(i); 0];
    
    st(i,:)=inv(aux)*v(i,:).';
    
    r(i,:)=E*st(i,:).';
    
   if i==round(length(s)/4)
      disp('Cuarto')
      toc
   end
   
      if i==round(length(s)/2)
      disp('Mitad')
      toc
      end
      
      if i==round(3*length(s)/4)
      disp('3/4')
      toc
      end
    
   
    
end


Is=st(:,2);
Flux_s=st(:,1);
Flux_r=r(:,1);
Ir=r(:,2);

Is_modulo=abs(Is);                     
Ir_modulo=abs(Ir);

Flux_s_mod=abs(Flux_s); % Nominal stator flux

Flux_r_mod=abs(Flux_r);

% Electromagnetic torque calculation 

Te=3/2*p*imag(conj(Flux_s).*Is);


%Losses calculation
% 
% Ploss_cu=Is_modulo.^2*Rs+Ir_modulo.^2*Rr;
% I_m=abs(Is+Ir);
% Ploss_iron=Rfe.*I_m.^2;
% Ploss_total=Ploss_cu+Ploss_iron;

Ploss_cu=3/2*(Is_modulo.^2*Rs+Ir_modulo.^2*Rr);
I_m=abs(Is+Ir);
Ploss_iron=3/2*(Rfe.*I_m.^2);
Ploss_total=Ploss_cu+Ploss_iron;

%% PART 2: Store minimum values for flux reference for specific values of wr and Te

discrete_Te_in=50:50:maxTorque;
discrete_wr_in=1:5:250;


discrete_Te=repmat(discrete_Te_in,size(discrete_wr_in))';
discrete_wr=repelem(discrete_wr_in,length(discrete_Te_in) )';

condition=Is_modulo<currentLimit & Te>=0;

optFlux=nan(size(discrete_wr));

for i=1:length(discrete_wr)
   
    %Set condition for speed
    
    condition_speed=wr>(discrete_wr(i)-0.2)&wr<(discrete_wr(i)+0.2);
    condition_Te=Te>(discrete_Te(i)-3)&Te<(discrete_Te(i)+3);
    
    if sum(condition&condition_speed&condition_Te)<3
        optFlux(i)=nan;
    else
    
    %Get spline curve for set of points flux vs Ploss
    
    spline_fit=fit(Flux_s_mod(condition&condition_speed&condition_Te),Ploss_total(condition&condition_speed&condition_Te),'smoothingspline','SmoothingParam',0.99999);
    
    %Find optimum flux (to minimize loss)
    
    datax=min(Flux_s_mod(condition&condition_speed&condition_Te)):0.01:max(Flux_s_mod(condition&condition_speed&condition_Te));
    datay=spline_fit(datax);
   
    [c I]=min(datay);
    min_flux=datax(I);
    
    end
    
    optFlux(i)=min_flux;
    
      if i==round(length(discrete_wr)/4)
      disp('Cuarto')
      toc
      end
   
      if i==round(length(discrete_wr)/2)
      disp('Mitad')
      toc
      end
      
      if i==round(3*length(discrete_wr)/4)
      disp('3/4')
      toc
      end
    
end

%% PART 3: Obtain best fitting surface to the optimum flux points

[xData, yData, zData] = prepareSurfaceData( discrete_wr, discrete_Te, optFlux );

% Set up fittype and options.
ft = fittype( 'poly33' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';

% Fit model to data.
[fitresult_poly, gof] = fit( [xData, yData], zData, ft, opts );


figure( 'Name', 'Optimum flux reference surface' );
% h = plot( fitresult, [xData, yData], zData);
h = plot(fitresult_poly);
hold on
scatter3(discrete_wr, discrete_Te, optFlux,'.');
% Label axes
xlabel( 'Rotor speed (rad/s)','FontSize',20)
ylabel( 'Torque (Nm)','FontSize',20)
zlabel( 'Optimum stator flux (Wb)','FontSize',20)
grid on
grid minor

h.FaceColor='none';
h.EdgeColor='interp';
colormap jet
% colorbar

set(gca,'FontSize',20)

surface_coefficients=coeffvalues(fitresult_poly);


end