%%% Mini Matlab Model for Acoustic-Gravity wave Interactions and Coupling
% (MiniMatMAGIC) Based on EP711 Example Code by J.B. Snively
% Tsunami-Generated Gravity Wave Example: Laughman et al. [JGR, 2017]

clc;
clear all;
close all;


fsave = 'C:\Class\PHYS8715\final project\duct1.mat\';

if(~exist(fsave,'file'))
    %%% Set Domain
    
    % Domain size (recommend MKS units)
    xmin=0;
    xmax=1200000;
    ymin=0;
    ymax=220000;
    xdomain=xmax-xmin;
    ydomain=ymax-ymin;
    
    % Cell/grid scale
    dx=1500;
    dy=1000;
    
    % Time Parameters...
    tmin=0;     % Initial Time
    tmax=8000;  % Final Time in seconds
    t=tmin;     % First time t=Tmin
    n=0;        % First Step n=0
    nframe=1;   % First Output
    skipT=90;   % Number of seconds to skip
    dCFL=0.8;   % Desired max CFL number (in either direction!)
    
    % Define positions at cell centers
    x_c=[xmin-3*dx/2:dx:xmax+3*dx/2];
    y_c=[ymin:dy:ymax+3*dy/2];
    % Construct Reference Mesh
    [X,Y]=meshgrid(x_c,y_c);
    [J,I]=size(X);
    
    % Initialize Variables, just because...
    F=cat(3,0.*X,0.*X,0.*X,0.*X);   % Fluxes for x-split
    G=F;                            % Fluxes for y-split
    Q=F;                            % Solution Variables
    S=F;                            % Source Terms
    
    %%% Set Problem
    
    % Initial Condition Parameters...
    
    % Physical Parameters for Euler or Acoustics Equations
    global g R P0 rho0 gamma C;
    g=9.55;
    R=287;
    Cp=(7/2)*R;
    P0=1E5;
    rho0=1.8;
    gamma=1.4;
    T0=P0./(rho0.*R);
    C=sqrt(gamma*P0/rho0);
    
    % Temperature and Bouyancy Profile
    [zz0,TT0] = read_msis;
    TT = interp1(zz0*1000,TT0,y_c);
    
    [~,Tp]=meshgrid(x_c,TT);
    
    
    NB2=(g^2/C^2)*(gamma-1);
    N2n=(g./(0.5.*(Tp(1:end-1,:)+Tp(2:end,:)))).*(diff(Tp)./dy+g./Cp);
    
    figure(111);
    subplot(1,2,1);
    plot(sqrt(N2n(:,1)),0.5.*(y_c(1:end-1)+y_c(2:end))./1000);
    xlabel('N (rad/s)'); ylabel('Altitude (km)');
    title('Buoyancy Frequency');
    subplot(1,2,2);
    plot(Tp(:,1),y_c./1000);
    xlabel('T (K)'); ylabel('Altitude (km)');
    title('Temperature Profile');
    
    % Pressure and Density Profiles
    H=R.*Tp./g;
    %rho0=rho0.*exp(-Y./H);
    P0=P0.*exp(-Y./H);
    for j=2:J
        for i=1:I
            P0(j,i)=P0(j-1,i).*exp(trapz(Y(j-1:j,i),-1./H(j-1:j,i),1));
        end
    end
    rho0=P0./(R.*Tp);
    C=max(max(sqrt(gamma.*P0./rho0)));
    
    % Plot Pressure, Density
    figure(222);
    semilogx(P0(:,1),y_c,rho0(:,1),y_c);
    legend('Pressure','Density');
    xlabel('Pressure (Pa) or Density (kgm^-3)');
    ylabel('Altitude (km)');% Viscosity and Thermal Diffusivity Profiles
    dynvisc=1.3e-5;
    Pr=0.7;
    kinvisc=dynvisc./rho0;
    tdiffus=gamma.*kinvisc./Pr;
    
    %Force-balance for gravity (i.e., hide errors in g!)
    g=(P0(2:end,1)-P0(1:end-1,1))./(-0.5*dy*(rho0(2:end,1)+rho0(1:end-1,1)));
    g=repmat(g,1,size(X,2)-1);
    figure(333);
    plot(g(:,1),y_c(1:end-1)./1000);
    xlabel('g (m/s^2)'); ylabel('Altitude (km)');
    
    % Source and other state parameters
    global wind bsource;
    
    % Wind Profile
    wind=0.*(Y.^0);
    
    % Boundary Source
    % Boundary Source
    bsource.yes=true;
    bsource.amp=1e-7;
    bsource.x=x_c;
    bsource.y=y_c;
    bsource.sigmat=2000;
    bsource.sigmax=30000;
    bsource.sigmay=3000;
    bsource.t0=600;
    bsource.x0=600000;
    bsource.y0=12000;
    bsource.omega  = 2*pi/(6.05*60);           % period 6.05 minutes
    bsource.kxx    = x_c.*2*pi/(31.4*1000); % lamda x, 31.4 km
    
    
    
    % Set Initial Condition
    Q(:,:,1)=rho0.*X.^0;
    Q(:,:,2)=rho0.*wind;
    Q(:,:,3)=0;
    Q(:,:,4)=(P0.*X.^0)/(gamma-1)+0.5*rho0.*wind.^2;
    
    
    %%% Evolve Over Time
    
    nframe=1;                   % First Frame
    T_arr(nframe)=tmin;         % Record Simulation Start Time
    dt=dCFL.*min(dx,dy)./(C);   % Set initial dt based on CFL Constraints
    Q=bc(Q,0);                  % Set BCs for first time
    Q_save(:,:,:,nframe)=Q;     % Store initial state
    
    % Set domain ranges for evaluations (extending into first boundary cells)
    iD=2:I-1;
    jD=2:J-1;
    
    tic;                        % Record CPU Time
    
    while t<tmax
        
        %
        % Euler Equations
        %
        
        % Euler x-split
        F=Fflux(Q);
        Qs(jD,iD,:)=0.5*(Q(jD,iD,:)+Q(jD,iD+1,:))-(dt/(2*dx))*(F(jD,iD+1,:)-F(jD,iD,:));
        F=Fflux(Qs);
        Q(jD,iD,:)=Q(jD,iD,:)-(dt/dx)*(F(jD,iD,:)-F(jD,iD-1,:));
        % apply BCs
        Q=bc(Q,t);
        % Euler y-split + Gravitational Source Terms
        G=Gflux(Q);
        S=Source(0.5*(Q(jD,iD,:)+Q(jD+1,iD,:)),g(jD,iD));
        Qs(jD,iD,:)=0.5*(Q(jD,iD,:)+Q(jD+1,iD,:))-(dt/(2*dy))*(G(jD+1,iD,:)-G(jD,iD,:))+(dt/2)*S(:,:,:);
        G=Gflux(Qs);
        S=Source(Qs,g);
        Q(jD,iD,:)=Q(jD,iD,:)-(dt/dy)*(G(jD,iD,:)-G(jD-1,iD,:))+dt*0.5*(S(jD,iD,:)+S(jD-1,iD,:));
        % apply BCs
        Q=bc(Q,t);
        
        %
        % Simplified Viscosity & Conduction -- Laplacian only.
        %
        
        % Substep
        dxymin=min(dx,dy);
        difmax=max(max(kinvisc));
        difCFL=0.15;    % Might need smaller for some problems
        % Substep Timestep
        Ndtr=max(int32((dt*difmax/(dxymin^2))/difCFL),1);
        dtr=dt/double(Ndtr);
        % Main Substepping Loop
        for ndtr=1:Ndtr
            % Viscosity and Conduction
            
            Q(:,:,2:3)=Q(:,:,2:3)./Q(:,:,1);
            Q(jD,iD,2:3)=Q(jD,iD,2:3)+kinvisc(jD,iD).*(dtr/dx^2).*(Q(jD,iD+1,2:3)-2.*Q(jD,iD,2:3)+Q(jD,iD-1,2:3));
            Q(:,:,2:3)=Q(:,:,2:3).*Q(:,:,1);
            % apply BCs
            Q=bc(Q,t+double(ndtr)*dtr);
            
            Q(:,:,2:3)=Q(:,:,2:3)./Q(:,:,1);
            Q(jD,iD,2:3)=Q(jD,iD,2:3)+kinvisc(jD,iD).*(dtr/dy^2).*(Q(jD+1,iD,2:3)-2.*Q(jD,iD,2:3)+Q(jD-1,iD,2:3));
            Q(:,:,2:3)=Q(:,:,2:3).*Q(:,:,1);
            % apply BCs
            Q=bc(Q,t+double(ndtr)*dtr);
            
            Q(:,:,4)=Q(:,:,4)-0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1);
            Q(:,:,4)=(Q(:,:,4).*(gamma-1))./(R.*Q(:,:,1))-Tp;
            Q(jD,iD,4)=Q(jD,iD,4)+tdiffus(jD,iD).*(dtr/dx^2).*(Q(jD,iD+1,4)-2.*Q(jD,iD,4)+Q(jD,iD-1,4));
            Q(jD,iD,4)=Q(jD,iD,4)+tdiffus(jD,iD).*(dtr/dy^2).*(Q(jD+1,iD,4)-2.*Q(jD,iD,4)+Q(jD-1,iD,4));
            Q(:,:,4)=((Q(:,:,4)+Tp).*R.*Q(:,:,1))./(gamma-1)+0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1);
            % apply BCs
            Q=bc(Q,t+double(ndtr)*dtr);
        end
        
        % Update time & timestep
        t=t+dt;
        n=n+1;
        dt=dCFL.*min(dx,dy)./(C+max(max(max(abs(Q(:,:,2:3)./Q(:,:,1)))))); % Not Fool-Proof (!)
        if ((isnan(dt)) || (dt==0))
            break;
        end
        
        % Store results
        if (abs(mod(int32(t),int32(skipT)))<dt) && ((int32(T_arr(nframe))<int32(t-2.*dt)))
            nframe=nframe+1;
            Q_save(:,:,:,nframe)=Q;
            T_arr(nframe)=t;
            disp(['dt=',num2str(dt),'(s); Time Step n=',num2str(n),'; Time t=',num2str(t),'(s)']);
        end
        
    end
    
    toc;                % Report CPU Time
    save(fsave);
else
    load(fsave);
end



%%% Plot Outputs

field_output=figure('rend','painters','pos',[10 10 1280 724]);
Nmax=size(T_arr,2); % Save Final Number of Timesteps

% Set ranges for plotting (including only physical domain cells)
iD=3:I-2;
jD=3:J-2;

mov1= VideoWriter('C:\Class\PHYS8715\final project\Duct1','MPEG-4');
mov1.Quality = 95;
mov1.FrameRate = 2;
open(mov1)
is_savemovie = 1;

for n=1:Nmax
    scale=sqrt(rho0./rho0(3,3));
    KE=0.5*(Q_save(jD,iD,2,n).^2+Q_save(jD,iD,3,n).^2)./Q_save(jD,iD,1,n);
    PR=(Q_save(jD,iD,4,n)-KE).*(gamma-1)-P0(jD,iD);
    TE=PR./(R.*Q_save(jD,iD,1,n));
    PR=PR./P0(jD,3:end-2);  % Relative Pressure Pert.
    rhoP=squeeze(100.*(Q_save(jD,iD,1,n)-rho0(jD,iD))./rho0(jD,iD));
    % plot results: Pressure
    subplot(4,1,1);
    imagesc(x_c(iD)/1000,y_c(jD)/1000,rhoP);
    axis xy; axis equal; axis([xmin/1000 xmax/1000 ymin/1000 ymax/1000]);
    caxis([-max(max(abs(rhoP))) max(max(abs(rhoP)))]);
    colorbar;
    xlabel('x (km)'); ylabel('Altitude (km)'); title(['Density Perturbation (%)']);
    % plot results: Temperature
    subplot(4,1,2);
    imagesc(x_c(iD)/1000,y_c(jD)/1000,TE);
    axis xy; axis equal; axis([xmin/1000 xmax/1000 ymin/1000 ymax/1000]);
    caxis([-max(max(abs(TE))) max(max(abs(TE)))]);
    colorbar;
    xlabel('x (km)'); ylabel('Altitude (km)'); title(['Temperature Perturbation (K)']);
    % plot results: Velocity
    subplot(4,1,3);
    us=scale(jD,iD).*(Q_save(jD,iD,2,n)./Q_save(jD,iD,1,n));
    imagesc(x_c(iD)/1000,y_c(jD)/1000,us);
    axis xy; axis equal; axis([xmin/1000 xmax/1000 ymin/1000 ymax/1000]);
    caxis([-max(max(abs(us))) max(max(abs(us)))]);
    colorbar;
    xlabel('x (km)'); ylabel('Altitude (km)'); title(['Scaled Horizontal Velocity (m/s)']);
    % plot results: Velocity
    subplot(4,1,4);
    ws=scale(jD,iD).*(Q_save(jD,iD,3,n)./Q_save(jD,iD,1,n));
    imagesc(x_c(iD)/1000,y_c(jD)/1000,ws);
    axis xy; axis equal; axis([xmin/1000 xmax/1000 ymin/1000 ymax/1000]);
    caxis([-max(max(abs(ws))) max(max(abs(ws)))]);
    colorbar;
    xlabel('x (km)'); ylabel('Altitude (km)'); title(['ScaledVertical Velocity (m/s)']);
    %suptitle(['Model Output Fields at Time t=',num2str(T_arr(n)),'s']);
    drawnow;
    if is_savemovie == 1
        F = getframe(field_output);
        writeVideo(mov1,F);
    end
end;

if is_savemovie == 1
    close(mov1);
end

%%% Functions

% Flux and Source Term Definitions

function F = Fflux(Q)
global g R P0 rho0 gamma C;
P = (gamma-1)*(Q(:,:,4)-(0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1)));
F(:,:,1) = Q(:,:,2);
F(:,:,2) = Q(:,:,2).*Q(:,:,2)./Q(:,:,1)+P;
F(:,:,3) = Q(:,:,2).*Q(:,:,3)./Q(:,:,1);
F(:,:,4) = (Q(:,:,4)+P).*Q(:,:,2)./Q(:,:,1);
end

function G = Gflux(Q)
global g R P0 rho0 gamma C;
P = (gamma-1)*(Q(:,:,4)-(0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1)));
G(:,:,1) = Q(:,:,3);
G(:,:,2) = Q(:,:,2).*Q(:,:,3)./Q(:,:,1);
G(:,:,3) = Q(:,:,3).*Q(:,:,3)./Q(:,:,1)+P;
G(:,:,4) = (Q(:,:,4)+P).*Q(:,:,3)./Q(:,:,1);
end

function S = Source(Q,g)
global bsource;
S(:,:,1) = 0.*Q(:,:,1);
S(:,:,2) = 0.*Q(:,:,1);

if bsource.yes
    phase  = bsource.omega.*t-bsource.kxx;
    velz   = bsource.amp.*cos(phase).*...
        exp(-(t-bsource.t0)^2./(2*bsource.sigmat^2)...
        -(bsource.x-bsource.x0).^2./(2*bsource.sigmax^2)...
        -(bsource.y-bsource.y0).^2./(2*bsource.sigmay^2));
    Q(1:2,:,3) = velz.*rho0(1:2,:);
else
    Q(1:2,:,3) = -Q(3,:,3).*(rho0(1:2,:)./rho0(3,:)).^(0.5);
end


S(:,:,3) = -Q(:,:,1).*g;
S(:,:,4) = -Q(:,:,3).*g;
end

% Boundary Conditions

function Q = bc(Q,t)
global g R P0 rho0 gamma C;
global wind bsource;
% Closed Bottom
Q(1:2,:,1) = rho0(1:2,:)+(Q(3,:,1)-rho0(3,:)).*(rho0(1:2,:)./rho0(3,:)).^(0.5);
Q(1:2,:,2) = rho0(1:2,:).*wind(1:2,:)+(Q(3,:,2)-rho0(3,:).*wind(3,:)).*(rho0(1:2,:)./rho0(3,:)).^(0.5);
if bsource.yes
    phase  = bsource.omega.*t-bsource.kxx;
    velz   = bsource.amp.*cos(phase).*...
        exp(-(t-bsource.t0)^2./(2*bsource.sigmat^2)...
        -(bsource.x-bsource.x0).^2./(2*bsource.sigmax^2));
    Q(1:2,:,3) = velz.*rho0(1:2,:);
else
    Q(1:2,:,3) = -Q(3,:,3).*(rho0(1:2,:)./rho0(3,:)).^(0.5);
end
Q(1:2,:,4) = P0(1:2,:)./(gamma-1)+(Q(3,:,4)-P0(3,:)./(gamma-1)-0.5*rho0(3,:).*wind(3,:).^2).*(rho0(1:2,:)./rho0(3,:)).^(0.5)+0.5*rho0(1:2,:).*wind(1:2,:).^2;
% Open Top
Q(end-1:end,:,1) = rho0(end-1:end,:)+(Q(end-2,:,1)-rho0(end-2,:)).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5);
Q(end-1:end,:,2) = rho0(end-1:end,:).*wind(end-1:end,:)+(Q(end-2,:,2)-rho0(end-2,:).*wind(end-2,:)).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5);
Q(end-1:end,:,3) = Q(end-2,:,3).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5);
Q(end-1:end,:,4) = P0(end-1:end,:)./(gamma-1)+(Q(end-2,:,4)-P0(end-2,:)./(gamma-1)-0.5*rho0(end-2,:).*wind(end-2,:).^2).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5)+0.5*rho0(end-1:end,:).*wind(end-1:end,:).^2;
% Periodic Sides
Q(:,2,:) = Q(:,end-2,:);
Q(:,1,:) = Q(:,end-3,:);
Q(:,end,:)  = Q(:,4,:);
Q(:,end-1,:)= Q(:,3,:);
end

function [z,T] = read_msis
fname = 'C:\Class\PHYS8715\final project\2001061512.msis';

fid = fopen(fname,'r');

for i=1:19
    tline = fgetl(fid);  % skip 4 lines
end

nmax = 300;
T = nan(1,nmax);
z = nan(1,nmax);
ii=1;
while(~feof(fid))
    tline = fgetl(fid);
    a = sscanf(tline,'%f %f %f %f %f %f %f');
    z(ii) = a(1);
    T(ii) = a(6);
    ii = ii+1;
end
fclose(fid);
z = z(1:ii-1);
T = T(1:ii-1);

end