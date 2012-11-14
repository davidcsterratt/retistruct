function plot_datapoints_polararea(Dss, phi0)
save_hold = ishold;

%% Draw circular gridlines
dl = 2*pi/90;                           % Resolution of circle
lambdas = 0:dl:2*pi;                    % Locations of longidtude lambda
phis = (-80:10:phi0) * pi/180;                     % Locations of lattitude phi
rs = sqrt(2 * (1 + sin(phis)));  % Phi transformed to r in
                                       % area-preserving polar coordinates
plot(cos(lambdas)' * rs, sin(lambdas)' * rs, 'b-')
hold on

%% Plot data points
%% First column is latitude (phi); second column is longitude
%% (lambda)

polararea(Dss.green, '.g')
polararea(Dss.red, '.r')

if ~save_hold
    hold off
end

