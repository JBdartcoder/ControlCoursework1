x = [0.3 : 0.00001: 0.7];
V = (774896.55*x.^3 - 1333971.97*x.^2 + 751982.01*x - 137991.22).^0.5;
plot(x, V), xlabel('x_1^{eq} (m)'), ylabel('Voltage (V)'), title('Problem B1 Graph'),
grid on

