function plot_outline_polar(Tss)
save_hold = ishold;
hold on
f = fieldnames(Tss)

for n = 1:length(fieldnames(Tss))
    disp(n)
    phi = Tss.(f{n})(:,1);
    lambda = Tss.(f{n})(:,2);
    polar(lambda, (phi .* 180/pi) + 90)
end
if ~save_hold
    hold off
end
