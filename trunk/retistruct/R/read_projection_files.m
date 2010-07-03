%% READ_PROJECTION_FILES - read in files produced from the folding procedure
%% DCS - 26/3/2010

%% Change this to the directory in which
%% SCGRIDCOO.csv, SCRED.csv, and SCGREEN.csv reside
datadir = 'data/Anatomy/ALU/M643-4/CONTRA';

%% Load in the grid coordinates
dat = csvread([datadir '/SCGRIDCOO.csv'], 1, 1);

PHIGRID    = dat(:,1);
LAMBDAGRID = dat(:,2);
XGRIDCOO = dat(:,3);
YGRIDCOO = dat(:,4);
COMPLETE = dat(:,5);
TOTALCEL = dat(:,6);
TOTALRED = dat(:,7);
TOTALGRE = dat(:,8);
TOTALDOU = dat(:,9);

%% Draw circular gridlines
phi0 = 50;
dl = 2*pi/90;                           % Resolution of circle
lambdas = 0:dl:2*pi;                    % Locations of longidtude lambda
phis = (-80:10:phi0) * pi/180;                     % Locations of lattitude phi
rs = sqrt(2 * (1 + sin(phis)));  % Phi transformed to r in
                                       % area-preserving polar coordinates
plot(cos(lambdas)' * rs, sin(lambdas)' * rs, 'b-')
hold on

%% Apply the area preserving transformation from lattitude to radius
%% in a polar plot
RGRID = sqrt(2 * (1 + sin(PHIGRID)));
THETAGRID = LAMBDAGRID;

% Plot the gridlines
plot(RGRID.*cos(THETAGRID), RGRID.*sin(THETAGRID), '.')

%% Load and plot red cell bodies
dat = csvread([datadir '/SCRED.csv'], 1, 1);
PHIRED = dat(:,1)
LAMBDARED = dat(:,2)

RRED = sqrt(2 * (1 + sin(PHIRED)));
THETARED = LAMBDARED;

plot(RRED.*cos(THETARED), RRED.*sin(THETARED), 'r.')

%% Load and plot green cell bodies
dat = csvread([datadir '/SCGREEN.csv'], 1, 1);
PHIGREEN = dat(:,1)
LAMBDAGREEN = dat(:,2)

RGREEN = sqrt(2 * (1 + sin(PHIGREEN)));
THETAGREEN = LAMBDAGREEN;

plot(RGREEN.*cos(THETAGREEN), RGREEN.*sin(THETAGREEN), 'g.')

