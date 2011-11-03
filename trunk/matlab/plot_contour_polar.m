function plot_contour_polar(cont)
save_hold = ishold;
polar(cont(:,2), cont(:,1)*180/pi+90, '.g')
hold on
if ~save_hold
    hold off
end
