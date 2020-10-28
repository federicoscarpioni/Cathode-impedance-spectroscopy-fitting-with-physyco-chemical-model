%% IMPEDANCE DATA VALLIDATION

%  ATTENTION: before lounching this script, make sure to have random
%  parametrs generated with the funstion rnd_par_gen.m

clear all
close all

options = optimset('Display','iter','PlotFcns',@optimplotfval);


%% Import data from text file.
% Notes: - Write down your own data file path and name
%        - Decimals must be point separated
%        - Data file must not have headers
%        - A file containing starting paramaters must be present
%          (use rnd_par_gen.m)
%
% Coposition of the data file:
% First coloumn must be frequency in Hz.
% Second coloumn must be real impedence in Ohm.
% Third coloumn must be minus imaginary impedance in Ohm.

path='C:Data\'; % Path of the data file
fileName='fileName'; % Data filename without extention

Z_d = load([path,fileName,'.txt']);
% Save frequency in a saperate vecotor for later use
omega=Z_d(:,1);  


%% Vector containing estimated parameters

% Import parameters value from txt file
par0 =  load([path,fileName,'.txt']);

% Upper and lower bounds of parameters
UB   = [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,...
    -inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,...
    -inf,-inf,-inf,-inf,inf,-inf,-inf,-inf,-inf];
LB   = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];


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
chisq=chisq/(length(omega)-17);

% Save minimum parameters in a txt file
dlmwrite(['par_min\valid_par_',fileName,'.txt'],par_min,'delimiter','\t')

% Calculate teoretical function with minimizzed parameters
Z_v(:,1)=omega;
F=fun(omega,par_min);
Z_v(:,2)=real(F);
Z_v(:,3)=-imag(F);
 
fprintf('Chi^2 %e \n',chisq)


%% Compute and plot fitting error residue 
% Residues of real part
res_re=100.*( Z_d(:,2)-Z_v(:,2))./Z_d(:,2); 
% Residues of imaginary part
res_im=100.*( Z_d(:,3)-Z_v(:,3))./Z_d(:,3); 
figure('Position',[0 0 1500 400])
hold on
plot((omega),res_re,'-k','LineWidth',1.5) %'Color',[0 0.4470 0.7410]
plot((omega),res_im,'--k','LineWidth',1.5) %'Color',[0.4660 0.6740 0.1880]
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
plot(Z_v(:,2),Z_v(:,3),'k','LineWidth',1.5) % Impedance from model
scatter(Z_d(:,2),Z_d(:,3),'ok') % Impedance from data
hold off
set(gca,'FontSize',24)
m=Z_d(length(omega),2);
axis equal
xlim([0 ceil(m)])
ylim([0 ceil(m)])
xlabel('Z_{Re} (\Omega)','Interpreter','tex','FontSize',24)
ylabel('-Z_{Im} (\Omega)','Interpreter','tex','FontSize',24)
legend('Voigt fit','Data','FontSize',24)
box on

      
%% Define the parametric function for general circuit with n Voigt elements
function y=fun(omega,par)
    v1=par(2)./(1+1i*omega*par(2)*par(3));
    v2=par(4)./(1+1i*omega*par(4)*par(5));
    v3=par(6)./(1+1i*omega*par(6)*par(7));
    v4=par(8)./(1+1i*omega*par(8)*par(9));
    v5=par(10)./(1+1i*omega*par(10)*par(11));
    v6=par(12)./(1+1i*omega*par(12)*par(13));
    v7=par(14)./(1+1i*omega*par(14)*par(15));
    v8=par(16)./(1+1i*omega*par(16)*par(17));
    v9=par(18)./(1+1i*omega*par(18)*par(19));
    v10=par(20)./(1+1i*omega*par(20)*par(21));
    v11=par(22)./(1+1i*omega*par(22)*par(23));
    v12=par(24)./(1+1i*omega*par(24)*par(25));
    v13=par(26)./(1+1i*omega*par(26)*par(27));
    v14=par(28)./(1+1i*omega*par(28)*par(29));
    v15=par(30)./(1+1i*omega*par(30)*par(31));
    v16=par(32)./(1+1i*omega*par(32)*par(33));
    v17=par(34)./(1+1i*omega*par(34)*par(35));
    
    y=par(1)+v1+v2+v3+v4+v5+v6+v7+v8+v9+v10+v11+v12+v13+v14+v15+v16+v17;
end