function polararea(dat, lty)
R = sqrt(2 * (1 + sin(dat(:,1))));
theta = dat(:,2);
plot(R.*cos(theta), R.*sin(theta), lty)

