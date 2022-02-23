% =========================================================================
%   This program is to examine A-stability for transport ODE
%   Numerical Scheme: trapezoidal time difference
% =========================================================================
set(0,'defaulttextfontsize',14);
set(0,'defaultaxesfontsize',14);
set(0,'DefaultAxesTickDir', 'out')
set(0,'DefaultFigureColormap',feval('jet'));

%%% define constants %%%
tau = 1;                                                    % wave period
omega = 2*pi/(tau);                                         % wave frequency
%%% define true solution parameters %%%
phi0 = 1;
deltt = 0.011; 
t0  = 0; t1 = 4*tau;
tt_true = t0:deltt:t1;
%%% select different delt %%%
col = {'b','g','r'};

lambda0 = [0.5,1];
figure(1);clf;
for ilambda = 1:numel(lambda0)
    lambda = lambda0(ilambda);
    %%% construct true solution %%%
    phitrue = phi0*exp((lambda+1i*omega)*tt_true);                  % true solution
    gamma = (lambda+1i*omega);
    delt_shred = abs(-2*lambda/(lambda^2+omega^2));
    delt0 = [delt_shred*1.1,delt_shred*0.11,delt_shred*0.011];  % choose different delt

    %%% plot true solution %%%
    figure(1);
    subplot(1,2,ilambda)
    plot(tt_true,real(phitrue),'k','linewidth',1);hold on;    
    %%% solve ODE %%%
    for idtt = 1:numel(delt0)
        delt = delt0(idtt);
        tt = t0:delt:t1;
        
        phi = nan(size(tt));
        phi(1) = phitrue(1);
        for itt = 2:numel(tt)
            phi(itt) = (phi(itt-1)*(1+0.5*gamma*delt))/(1-0.5*gamma*delt);
        end
        plot(tt,real(phi),col{idtt},'linewidth',1);
    end
    xlabel('Time');ylabel('\phi');title(['\lambda = ',num2str(lambda)]);grid on;
    lg = legend('True','\Deltat = 1.1\Deltat(thres)','\Deltat = 0.11\Deltat(thres)',...
        '\Deltat = 0.011\Deltat(thres)','location','northwest');
    set(lg,'box','off','fontsize',12);
end