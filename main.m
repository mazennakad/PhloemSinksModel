clear all
close all
addpath('./fcn/');
addpath('./solver/');
addpath('./initial/');
path = "./";
data = readtable(path + "BCI_hydro_Cohort_hourly.csv");
duration = 24; % duration of the simulation, 24 hours

% Get a c-nu relation at ambient temperature
c = linspace(0,2000,1000)';
[~,bb,rrs] = MLR_Viscosity(c,293.15,1000);

%% Grid and time step
n   = 100;   % nb of cells      
dz  = 1/n;   % size of the cell
dt  = 1e0;   % timestep, nondimensional
nts = 14400; % nb of time steps within an hour (this makes it 1 hour)

%% Constants
k     = 5e-14;                       % membrane permeability m/(Pa s)
a     = 1e-5;                        % phloem thickness m
T     = 293.15;                      % temperature K
Rg    = 8.3145;                      % Ideal gas constant J/(mol K)
Mw    = 342.2965;                    % Molecular mass of sucrose g/mol
rho   = 1e3;                         % Density of water kg/m3
g     = 9.81;                        % Gravitational constant m/s2
D     = 4e-10;                       % Sucrose diffusivity in water m2/s

%%% Declare variables
c    = zeros(n,nts);   % Sucrose concentration at the center of the cell
u    = zeros(n-1,nts); % Axial velocity at the face of the cell
p    = zeros(n,nts);   % Dynamic fluid pressure at the center of the cell
v    = zeros(n,nts);   % Radial velocity at the center of the cell
nu   = zeros(n,nts);   % viscosity at the center of the cell

%%% For saving
cs         = zeros(n,duration);   % Sucrose mass at the center of the cell
ps         = zeros(n,duration);   % Dynamic fluid pressure at the center of the cell
vs         = zeros(n,duration);   % Radial velocity at the center of the cell
us         = zeros(n-1,duration); % Axial velocity at the face of the cell
nus        = zeros(n,duration);   % Dynamic viscosity
cc         = zeros(n,duration);   % Sucrose concentration
pt         = zeros(n,duration);   % Total pressure (pd + ph)
cm         = zeros(1,duration);   % Average concentration in the phloem 
dp         = zeros(1,duration);   % Dynamic pressure gradient along the phloem
psi_xylem  = zeros(1,duration);   % Xylem water potential
psi_phloem = zeros(1,duration);   % Phloem water potential
psi_ratio  = zeros(1,duration);   % Ratio
psi_diff   = zeros(1,duration);   % Difference
maxU       = zeros(1,duration);  % Maximum longitudinal velocity
cfl        = zeros(1,duration);  % Courant Number within each hour

% Save biomeE input
sinks = zeros(1,duration);  % total stem sink
sinksg = zeros(1,duration); % Growth stem sink
sinksr = zeros(1,duration); % Respiration stem sink
sinkr = zeros(1,duration);  % total root sink
source = zeros(1,duration); % total source at the leaf
leafw = zeros(1,duration);  % leaf water potential
rootw = zeros(1,duration);  % root water potential



fname = "./../results/paper_new1/c0ss.xlsx"; % Filename to save

m     = 0;        % k = 0 constant viscosity, 
                  % k = 1 variable visvosity
k1    = 0;        % k1 = 0 constant sink, 'CS'
                  % k1 = 1 linear sink (highest at bottom), 'LSb'
                  % k1 = 2 linear sink (highest at top), 'LSt'
                  % k1 = 3 cts respiration and linear growth (highest at
                  % bottom), 'LGSb' 
                  % k1 = 4 cts respiration and linear growth (highest at
                  % top), 'LGSt'

rootonly = 'False'; % for the 'RS' case
reduced  = 'False'; % for the reduced sink case

for j = 1:duration

    if j == 1
        dx    = data.dbh(j);                 % diameter at breast height 
        L     = data.height(j);              % tree length 
        psi0  = - data.Psi_L(j);             % Leaf water potential Pa
    end

    if strcmpi(rootonly, 'true')
        Sstemr = 0; Sstemg = 0; Sstem = 0;
        Sroot = ( data.Resp_r(j) ...
        + data.Grow_Resp_r(j) ...
        + data.Growth_r(j) ...
        + data.Resp_s(j) ...
        + data.Grow_Resp_s(j) ...
        + data.Growth_s(j) )/(3600*a*pi*dx); % Root sucrose sink g/(m s)
        k1 = 0;
    else
        Sstem = ( data.Resp_s(j) ...
            + data.Grow_Resp_s(j) ...
            + data.Growth_s(j) )/(3600*L*pi*dx); % Stem sucrose sink g/(m s)
        Sstemr = ( data.Resp_s(j) ...
            + data.Grow_Resp_s(j) )/(3600*L*pi*dx); % Stem sucrose sink g/(m s)
        Sstemg = data.Growth_s(j)/(3600*L*pi*dx); % Stem sucrose sink g/(m s)
        Sroot = ( data.Resp_r(j) ...
            + data.Grow_Resp_r(j) ...
            + data.Growth_r(j) )/(3600*a*pi*dx); % Root sucrose sink g/(m s)
    end

    if strcmpi(reduced, 'true')
        Sstem = Sstem/10;
        Sstemr = Sstemr/10;
        Sstemg = Sstemg/10;
        Sroot = Sroot/10;
    end
   
    Sleaf = ( Sstem*L + Sroot*a )/a;     % Sucrose source g/(m s)
    

    %%% Non-dimensional and scaling quantities
    c0 = psi0*Mw*pi*dx/(Rg*T);           % Sucrose scaling g/m2 
                                         % (assumed by balancing osmotic potential to leaf water potential)
                                         % it wont matter how you balance it
    cw = c0/(Mw*pi*dx);                  % mol/m3
    nu0 = exp(bb(1) + bb(2)*cw ...
        + bb(3)*(cw^2) + bb(4)*(cw^3) ...
        + bb(5)*(cw^4) );                % Viscosity of sucrose Pa*s
    os = Rg*T*c0/(Mw*pi*dx);             % Osmosis Pa
    es = a/L;                            % aspect ratio
    v0 = k*os;                           % Radial velocity scaling m/s
    u0 = v0/es;                          % Axial velocity scaling m/s
    p0 = L*nu0*u0/(a^2);                 % Pressure scaling Pa
    t0 = (a^2)/D;                        % diffusive timescale s
    G  = rho*g*L/os;                     % Gravity / osmosis
    X0 = psi0/os;                        % Xylem water potential / osmosis
    Mu = k*nu0*(L^2)/(a^3);              % Munch nb, axial / radial resistance
    Pe = a*v0/D;                         % Peclet number, advection / diffusion (radial)
    Ss = Sstem/(u0*c0);                  % Stem sink / advection
    Sl = Sleaf/(u0*c0);                  % Leaf source / advection
    Sr = Sroot/(u0*c0);                  % Root sink / advection

    Ss1 = Sstemr/(u0*c0);
    Ss2 = Sstemg/(u0*c0);

    % Save the biomeE input
    sinks(j) = Sstem;
    sinksg(j) = Sstemg;
    sinksr(j) = Sstemr;
    sinkr(j) = Sroot;
    source(j) = Sleaf;
    leafw(j) = data.Psi_L(j);
    rootw(j) = data.Psi_W(j);


    %% Initial condition
    Psi = linspace(data.Psi_L(j),data.Psi_W(j),n)'./psi0; % Xylem water potential at the center of the cell

    if j==1
        [c(:,1),u(:,1),nu(:,1),v(:,1),p(:,1)] ...
            = initial(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,Psi, ...
            cw,bb,m,k1,Ss1,Ss2,dt);
        % checkmass = dz*( sum(c(2:n-1,1)) + (1/2)*( c(1,1) + c(n,1) ) )

    else
        c(:,1)  = co;
        u(:,1)  = uo;
        v(:,1)  = vo;
        p(:,1)  = po;
        nu(:,1) = nuo;

    end

    %% Solve in time
    for i = 2:nts
        [c(:,i),u(:,i),v(:,i),p(:,i),nu(:,i)] = ...
            iterate(n,dz,Mu,X0,Pe,Sl,Ss,Sr,es,...
            Psi,c(:,i-1),u(:,i-1),v(:,i-1),p(:,i-1),nu(:,i-1),...
            cw,bb,m,k1,1,dt,c(:,i-1),Ss1,Ss2);
        % checkmass = dz*( sum(c(2:n-1,i)) + (1/2)*( c(1,i) + c(n,i) ) )
    end


    % Save the data at previous time
    co              = c(:,end);
    uo              = u(:,end);
    vo              = v(:,end);
    po              = p(:,end);
    nuo             = nu(:,end);

    % make the variables dimensional
    cs(:,j)         = (1e-3).*co.*c0;   % Kg/m2       
    ps(:,j)         = (1e-6).*po.*p0;   % MPa
    vs(:,j)         = vo.*v0;           % m/s
    us(:,j)         = uo.*u0;           % m/s
    nus(:,j)        = (1./nuo).*nu0;    % Pa s


    % variables for saving
    cc(:,j)         = (1e3).*cs(:,j)./(Mw*pi*dx);    % mol/m3
    cm(j)           = dz*( sum(cc(2:n-1,j)) + (1/2)*( cc(1,j) + cc(n,j) ) ) ; % Average mol/m3
    pt(:,j)         = (1e-6).*(p0.*po(:) + dz*L*rho*g.*((1:n)' - 1/2)); % total presusre in MPa (dynamic + hydrostatic)
    pp              = (1e-6).*(p0.*po(:) - os.*co(:));  % phloem water potential MPa
    px              = (1e-6)*Psi.*psi0; % xylem water potential in MPa
    dp(j)           = (1e-6)*p0*(po(1) - po(end)); % pressure gradient in MPa
    psi_phloem(1,j) = dz*( sum(pp(2:n-1)) + (1/2)*( pp(1) + pp(n) ) ); % average phloem water potential MPa
    psi_xylem(1,j)  = dz*( sum(px(2:n-1)) + (1/2)*( px(1) + px(n) ) ); % average xylem water potential MPa
    psi_diff(1,j)   = psi_xylem(j) - psi_phloem(j); % average water potential disequilibrium Pa
    psi_ratio(1,j)  = psi_xylem(j)/psi_phloem(j);
    maxU(j)         = dz.*(sum(uo(2:n-2)) + (1/2).*( uo(1) + uo(n-1) )).*u0;
    cfl(j)          = max(max(u))*dt/dz;

end

% saveData(cs,us,vs,ps,nus,dp,psi_phloem,psi_xylem,pt,maxU,cc,cm,psi_diff,cfl,fname)

figure
plot(1:24,sinkr.*3600,'-','Color',[0 0 0],'Linewidth',3)
ax = findobj(gcf,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$s_r (g m^{-1} hr^{-1})$}','Interpreter','latex')


figure
plot(1:24,sinks.*3600,'-','Color',[0 0 0],'Linewidth',3)
ax = findobj(gcf,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$s_s (g m^{-1} hr^{-1})$}','Interpreter','latex')

figure
plot(1:24,sinksr.*3600,'-','Color',[0 0 0],'Linewidth',3)
ax = findobj(gcf,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$f_r (g m^{-1} hr^{-1})$}','Interpreter','latex')

figure
plot(1:24,leafw,'-','Color',[0 0 0],'Linewidth',3)
ax = findobj(gcf,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$\psi_l (MPa)$}','Interpreter','latex')





