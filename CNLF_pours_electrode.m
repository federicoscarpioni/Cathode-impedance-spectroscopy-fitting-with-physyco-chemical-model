%% NONLINEAR FITTING OF POURUS ELECTRODES IMPEDANCE
% Ipedance of porus electrodes is fitted using a Complex Non-Linear
% Regression.


clear all
close all
options = optimset('Display','iter','PlotFcns',@optimplotfval);

%% Import data from text file.
% Notes: - Write down your own data file path and name
%        - Decimals must be point separated
%        - Data file must not have headers
%
% Structure of the data file:
% First coloumn must be frequency in Hz.
% Second coloumn must be real impedence in Ohm.
% Third coloumn must be minus imaginary impedance in Ohm.

path='C:Data\'; % Path of your data file
fileName='fileName'; % Data filename without extention

Z_d = load([path,fileName,'.txt']);
% Save frequency in a saperate vecotor for later use
omega=Z_d(:,1);  


%% Initial estimated parameters
%  Manually written by the user, for the first minimization step

R_HFR =  0.2 ;        
R_cont = 4;     
Q_cont = 1e-3;   
a_cont = 0.85;    
R_pore =  6 ;   
R_el = 1e-8;         
R_CT = 4 ;       
Q_CT = 1e-3;      
a_CT = 0.85;        
W = 1.24e-02 ;       


%% Vector containing estimated parameters
% Use previous defined parameters:
par0 = [R_HFR;...
        R_cont;...
        Q_cont;...
        a_cont;...
        R_pore;...
        R_el;...
        R_CT;...
        Q_CT;...
        a_CT;...
        W]; 
    
% Import parameters from a file (uncomment after the first minimization step):
% par0 =  load(['par_min\fit_par_',fileName,'.txt']);

% Upper and lower bounds of parameters
UB   = [ -inf,...
         -inf,...
         -inf,...
         1,... %a_cont
         -inf,...
         -inf,...
         -inf,...
         -inf, ...
         1,... %a_CT
         -inf];
LB   = [ 0,...
         0,...
         0,...
         0.5,... %a_cont
         0, ...
         0,...
         0,...
         0,...
         0.5,... %a_CT
         0];


%% Minimization of sum of the sqaures
% Sum of sqares with a modulus weighting
sum_of_sqr= @(omega,par) ...
            sum( ((real(fun(omega,par))-Z_d(:,2)).^2 ... % Real
                 +(-imag(fun(omega,par))-Z_d(:,3)).^2) ... % Imaginary
                 ./(Z_d(:,2).^2+Z_d(:,3).^2) ... % Weighting factor
                );

% Minimization function
[par_min, chisq] = fminsearchbnd(@(par) ...
    sum_of_sqr(omega,par),par0,LB,UB,options);

% Normalize chi^2 on number of data point minus the number of free parameters
chisq=chisq/(length(omega)-10); 

% Save minimum parameters in a txt file
dlmwrite([path,fileName,'.txt'],par_min,'delimiter','\t')

% Calculate teoretical function with minimizzed parameters
Z_t(:,1)=omega;
F=fun(omega,par_min);
Z_t(:,2)=real(F);
Z_t(:,3)=-imag(F);


%% Compute and plot fitting error residue 
% Residues of real part
res_re=100.*( Z_d(:,2)-Z_t(:,2))./Z_d(:,2); 
% Residues of imaginary part
res_im=100.*( Z_d(:,3)-Z_t(:,3))./Z_d(:,3); 
figure('Position',[0 0 1500 400])
hold on
plot((omega),res_re,'-','Color',[0 0.4470 0.7410],'LineWidth',1.5) %'Color',[0 0.4470 0.7410]
plot((omega),res_im,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5) %'Color',[0.4660 0.6740 0.1880]
hold off
set(gca,'xscale','log','FontSize',18)
xlabel('\omega (Hz)','Interpreter','tex','FontSize',24)
ylabel('Residuals %','FontSize',24)
xlim([omega(length(omega)) omega(1)])
ylim([-1 1])
legend('Real','Imaginary')
box on


%% Plotting fitting
figure('Position',[0 0 800 800])
hold on
plot(Z_t(:,2),Z_t(:,3),'k','LineWidth',1.5) % Impedance from model
scatter(Z_d(:,2),Z_d(:,3),'o') % Impedance from data
hold off
set(gca,'FontSize',24)
m=Z_d(length(omega),2);
axis equal
xlim([0 ceil(m)])
ylim([0 ceil(m)])
xlabel('Z_{Re} (\Omega)','Interpreter','tex','FontSize',24)
ylabel('-Z_{Im} (\Omega)','Interpreter','tex','FontSize',24)
legend('Fit','Data','FontSize',24)
box on

    
%% Print in command line minimum parameters
fprintf([ 'Minimized parameters:\n'...
    'R_HFR  = %.2e Ohm \n' ...
    'R_cont = %.2e Ohm\n' ...
    'Q_cont = %.2e F*s\n' ...
    'a_cont = %.2e\n' ...
    'R_pore = %.2e Ohm\n' ...
    'R_el   = %.2e Ohm\n' ...
    'R_CT   = %.2e Ohm\n' ...
    'Q_CT   = %.2e F*s\n' ...
    'a_CT   = %.2e\n' ...
    'W      = %.2e Ohm*s^(-1/2)\n' ...
    'Chi^2  = %f \n'],...
    par_min(1), par_min(2), par_min(3),par_min(4),par_min(5),par_min(6),...
    par_min(7),par_min(8),par_min(9),par_min(10),chisq)


%% Function for the regression
function y = fun(omega,par)
    RHFR  = par(1); % High frequency intercept
    Rcont = par(2);
    Qcont = par(3);
    acont = par(4);
    Rpore = par(5);
    Rel   = par(6); % Electric resistance of metal collector
    Rct   = par(7);
    Qct   = par(8);
    act   = par(9);
    W      = par(10);
    
    ni = sqrt((Rct*Qct*(1i*omega).^act +1)*(Rpore+Rel)/Rct);
    p = Rpore/(Rpore+Rel);
    s = 1-p;
    Z_str = sqrt(Rct*(Rpore+Rel)./(Rct*Qct*(1i*omega).^act +1)); % Z*
    Zpore = s*p*(Rpore+Rel) + Z_str.*(1+2*p*s.*(sqrt(1-(tanh(ni)).^2)-1))./(tanh(ni));
    Zcon = Rcont*(Rcont*Qcont*(1i*omega).^acont + 1).^-1;
    Zw = W./sqrt(omega)-1i.*(W./sqrt(omega));
    y = RHFR + Zpore + Zcon + Zw ; % Z_tot

end