function [I] = defIntegral()

f = @(x) 4 - x.^2;
a = 0; b = 2;
N = 100;

c = -inf;
for i = a:b
    if (f(i) > c)
        c = f(i);
    end
end

X = a:0.1:b;
Y = f(X);
plot(X,Y);
hold on;
%plot(1,1,'--.r');

n = 0;
for i = 1:N
   x = randInRange(a,b);
   y = randInRange(0,c); 
   if (y <= f(x))
       plot(x,y,'--.b');
       n = n + 1;
   else
       plot(x,y,'--.r');
   end
end

I = (b-a) * c * n / N;

end

