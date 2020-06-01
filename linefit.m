function [p,R] = linefit(x_data, y_data)
[p,S] = polyfit(x_data, y_data,1);
R = 1 - (S.normr/norm(y_data - mean(y_data)))^2;
f = polyval(p,x_data);
plot(x_data, y_data, 'o', x_data, f, '-')
legend('data','linear fit','FontSize',20)
end