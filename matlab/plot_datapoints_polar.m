function plot_datapoints_polar(Dss)
save_hold = ishold;
polar(Dss.green(:,2), Dss.green(:,1)*180/pi+90, '.g')
hold on
polar(Dss.red(:,2),   Dss.red(:,1)*180/pi  +90, '.r')
if ~save_hold
    hold off
end
