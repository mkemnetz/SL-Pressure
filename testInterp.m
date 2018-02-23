%%
% startPoint = 595258/238;
% overflow   = mod(length(forcingSignal), (startPoint));
% 
% x = 1:length(forcingSignalInterp); 
% xq = startPoint:startPoint/1000:length(forcingSignal)-overflow;
% 
% 
% 
% 
% v = forcingSignalInterp.*44.5300;
% forcingSignalInterp_b = interp1(x,v,xq, 'spline');

startPoint = 49405/1898;
wavelength = 107393414/906295;
overflow   = mod(length(forcingSignal), wavelength);

x = 1:length(forcingSignal); 
xq = startPoint:wavelength:length(forcingSignal)-overflow;




v = forcingSignal;
forcingSignalInterp_b = interp1(x,v,xq, 'spline');
forcingSignalInterp_b_scaled = forcingSignalInterp_b.*44.5300;

figure();
plot(forcingSignalInterp_b_scaled);