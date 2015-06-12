classdef MonteCarlo
    % Integration by Monte Carlo method
    
    methods(Static)
        
        %% randInRange: Returns random value in range [a,b].
        function [x] = randInRange(a,b)
            n = length(a);
            x = zeros(1,n);
            for i = 1:n
                x(i) = a(i) + (b(i) - a(i)) .* rand();
            end
        end
        
        %% checkPoint: Checks, is point x in G area, or not.
        function [inArea] = checkPoint(x,G)
            n = length(x);
            inArea = 1;
            for i = 1:n
                inArea = inArea && G{i}(x);
            end
        end
        
        %% ndIntegral: Computing n-dimensional definite integral 
         % at G area which no more than 
         % n-dimensional parallelepiped with properties a and b.
        function [I] = ndIntegral(f,a,b,G,N)
            t_beta = 3; % beta = 0.997
            fSum = 0; % sum of computed values f(x)
            fSumSquared = 0; % squared sum of f(x)
            n = 0; % amount of points found in G
            
            % generating N random vectors x
            for i = 1:N
                x = MonteCarlo.randInRange(a,b);
                % check conditions, x must be in [a,b]
                inArea = MonteCarlo.checkPoint(x,G); % bool value
                % adding f(x)
                if (inArea)
                    fSum = fSum + f(x);
                    fSumSquared = fSumSquared + f(x)^2;
                    n = n + 1;
                end
            end
           
            % computing n-dimensional volume of figure
            V = prod(b - a);
            
            % computing integral value
            I = V * fSum / N;
            
            % computing Sigma
            sigma = n / N
            
            % computing standard deviation
            sd1 = fSumSquared/n - (fSum/n)^2
            sd2 = sigma * (1 - sigma)
            
            % computing error
            error = V * t_beta * (sigma*sqrt(sd1)/sqrt(n) + ...
                                  abs(fSum/n)*sqrt(sd2)/sqrt(N))
        end
        
        %% TESTS
        
        %% test1d: computing 1-dimensional definite integral
        function [I] = test1d(N)
            n = 1;
            f = @(x) 1/sqrt(2*pi) * exp(-(x^2)/2);
            
            % conditions of G area
            x1Cond = @(x) (0<=x(1) && x(1)<=3);
            G = {x1Cond};
            
            % limitations for every element in x
            a = zeros(1,n);
            b = zeros(1,n);
            a(1) = 0; b(1) = 3;
            
            I = MonteCarlo.ndIntegral(f,a,b,G,N);
        end
        
        %% test2d: computing 2-dimensional definite integral
        function [I] = test2d(N)
            n = 2;
            f = @(x) x(1) + x(2);
            
            % conditions of G area
            x1Cond = @(x) (0<=x(1) && x(1)<=2);
            x2Cond = @(x) (x(1)^2<=x(2) && x(2)<=2*x(1));
            G = {x1Cond,x2Cond};
            
            % limitations for every element in x
            a = zeros(1,n);
            b = zeros(1,n);
            a(1) = 0; b(1) = 2;
            a(2) = 0; b(2) = 4;
            
            I = MonteCarlo.ndIntegral(f,a,b,G,N);
        end
        
        %% test3d: computing 3-dimensional definite integral
        function [I] = test3d(N)
            n = 3;
            f = @(x) 10*x(1);
            
            % conditions of G area
            x1Cond = @(x) (0<=x(1) && x(1)<=1);
            x2Cond = @(x) (0<=x(2) && x(2)<=sqrt(1-x(1)^2));
            x3Cond = @(x) (0<=x(3) && x(3)<=((x(1)^2+x(2)^2)/2));
            G = {x1Cond,x2Cond,x3Cond};
            
            % limitations for every element in x
            a = zeros(1,n);
            b = zeros(1,n);
            a(1) = 0; b(1) = 1;
            a(2) = 0; b(2) = 1;
            a(3) = 0; b(3) = 1;
            
            I = MonteCarlo.ndIntegral(f,a,b,G,N);
        end

        %% test6d: computing 6-dimensional definite integral
        function [I] = test6d(N)

        end
        
    end
    
end

