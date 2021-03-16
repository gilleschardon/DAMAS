function plothist(estp, p, bounds, bins, color, names, namem)

err = 10*log10(estp) - 10*log10(p);

pc16 = quantile(err, 0.16);
m = median(err);
pc84 = quantile(err, 0.84);

hold on

histogram( err, bins, 'FaceColor', color, 'LineStyle', 'none')
grid on
xlim(bounds)

Y = ylim();
plot([pc16 pc16], Y, 'k--', 'linewidth', 2)
plot([m m], Y, 'k', 'linewidth', 2)
plot([pc84 pc84], Y, 'k--', 'linewidth', 2)



xlabel(names)
ylabel(namem)

end
