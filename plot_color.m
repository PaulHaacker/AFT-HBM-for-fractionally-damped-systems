function plot_color(x,y,flag)
plot(x,y)
hold on
for i = 1:length(x)
 if flag(i)
     plot(x(i),y(i),'or')
 else
     plot(x(i),y(i),'+b')
 end
end
hold off