function plot_landmarks_polar(Sss) 
save_hold = ishold;
hold on
f = fieldnames(Sss);

for n = 1:length(fieldnames(Sss))
    phi = Sss.(f{n})(:,1);
    lambda = Sss.(f{n})(:,2);
    polar(lambda, (phi .* 180/pi) + 90)
end
if ~save_hold
    hold off
end


