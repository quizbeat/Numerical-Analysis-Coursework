classdef MonteCarlo
    % Integration by Monte Carlo method
    % 	tasks to do:
    %	1. find beta value
    %	2. write function to find a and b limits
    
    methods(Static)
        
        function [x] = randInRange(a,b)
        	% returns random value in range [a,b]
            n = length(a);
            x = zeros(1,n);
            for i = 1:n
                x(i) = a(i) + (b(i) - a(i)) .* rand();
            end
        end
        
        function [inArea] = checkPoint(x,G)
            n = length(x);
            inArea = 1;
            for i = 1:n
                inArea = inArea && G{i}(x);
            end
        end
        
        function [I] = nDimIntegral(f,a,b,G,N)
            % computing definite integral at G area
            % which no more than n-dimensional parallelepiped
            % with properties a and b
            n = length(a); % dimension of random vector x
            fSum = 0; % sum of computed values f(x)
            count = 0; % amount of points found in G
            % generating N random vectors x
            for i = 1:N
                x = MonteCarlo.randInRange(a,b);
                % check conditions, x must be in [a,b]
                % inArea = isequal(a<=x,x<=b,ones(1,n)); wrong check
                inArea = MonteCarlo.checkPoint(x,G);
                % adding f(x)
                if (inArea)
                    fSum = fSum + f(x);
                    count = count + 1;
                end
            end
            % computing n-dimensional volume of figure
            V = prod(b - a);
            % computing integral value
            I = V * fSum / N;
        end
        
        function [I] = test3d(N)
            % 3-dimensional integral
            % SSS x^2 dxdydz at G area
            n = 3;
            f = @(x) x(1)^2;
            
            x1Cond = @(x) (0<=x(1) && x(1)<=1);
            x2Cond = @(x) (0<=x(2) && x(2)<=1-x(1));
            x3Cond = @(x) (0<=x(3) && x(3)<=10*(x(1)+3*x(2)));
            
            G = {x1Cond,x2Cond,x3Cond};
            
            a = zeros(1,n);
            b = zeros(1,n);
            a(1) = 0; b(1) = 1;
            a(2) = 0; b(2) = 1;
            a(3) = 0; b(3) = 40;
            
            I = MonteCarlo.nDimIntegral(f,a,b,G,N);
        end
        
    end
    
end

