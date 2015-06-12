classdef MonteCarlo
    % Integration by Monte Carlo method
    %   
    
    methods(Static)
        
        function [x] = randInRange(a,b)
          % returns random value in range [a,b]
            n = length(a);
            x = zeros(1,n);
            for i = 1:n
                x(i) = a(i) + (b(i) - a(i)) .* rand();
            end
        end
        
        function [min,max] = findMinMax1Arg(lhs,rhs,x_min,x_max)
            h = 0.05;
            min = inf;
            max = -inf;
            for x = x_min:h:x_max
               if (lhs(x) < min)
                   min = lhs(x);
               end
               if (rhs(x) > max)
                   max = rhs(x);
               end
            end
        end
        
        function [min,max] = findMinMax2Arg(lhs,rhs,x_min,x_max,y_min,y_max)
            h = 0.05;
            min = inf;
            max = -inf;
            for x = x_min:h:x_max
                for y = y_min:h:y_max
                    if (lhs(x,y) < min)
                        min = lhs(x,y);
                    end
                    if (rhs(x,y) > max)
                        max = rhs(x,y);
                    end
                end
            end
        end
        
        function [I] = int1d(f,x_lhs,x_rhs,N)
          % Sf(x)dx = I
            % single definite integral with conditions:
            % x_lhs <= x <= x_rhs
            
            sum = 0;
            for i = 1:N
               x = MonteCarlo.randInRange(x_lhs,x_rhs);
               sum = sum + f(x);
            end
            
            I = (x_rhs-x_lhs) * sum / N;
        end

        function [I] = int2d(f,x_lhs,x_rhs,y_lhs,y_rhs,N)
          % SSf(x,y)dxdy = I
            % double definite integral with conditions:
            % x_lhs <= x <= x_rhs
            % y_lhs <= y <= y_rhs
            
          x_min = x_lhs;
          x_max = x_rhs;

          [y_min,y_max] = MonteCarlo.findMinMax1Arg(y_lhs,y_rhs,x_min,x_max);

          sum = 0;

          for i = 1:N
            x = MonteCarlo.randInRange(x_min,x_max);
            y = MonteCarlo.randInRange(y_min,y_max);
            if (x_lhs <= x && x <= x_rhs && ...
              y_lhs(x) <= y && y <= y_rhs(x))
              sum = sum + f(x,y);
            end
          end

          S = (x_max-x_min) * (y_max-y_min);
          I = S * sum / N;
        end
        
        function [I] = int3d(f,x_lhs,x_rhs,y_lhs,y_rhs,z_lhs,z_rhs,N)
            % SSSf(x,y,z)dxdydz = I
            % triple definite integral with conditions:
            % x_lhs <= x <= x_rhs
            % y_lhs <= y <= y_rhs
            % z_lhs <= z <= z_rhs
            
            x_min = x_lhs; 
            x_max = x_rhs;
            
            [y_min,y_max] = MonteCarlo.findMinMax1Arg(y_lhs,y_rhs,x_min,x_max);
            
            [z_min,z_max] = MonteCarlo.findMinMax2Arg(z_lhs,z_rhs,x_min,x_max,y_min,y_max);
                
            sum = 0;
            
            for i = 1:N
              % generating vector of random values
                x = MonteCarlo.randInRange(x_min,x_max);
                y = MonteCarlo.randInRange(y_min,y_max);
                z = MonteCarlo.randInRange(z_min,z_max);
                % check if point (x,y,z) in G area
                if (x_lhs <= x && x <= x_rhs && ...
                    y_lhs(x) <= y && y <= y_rhs(x) && ...
                    z_lhs(x,y) <= z && z <= z_rhs(x,y))
                    sum = sum + f(x,y,z);
                end
            end
            
            V = (x_max-x_min) * (y_max-y_min) * (z_max-z_min);
            I = V * sum / N;
        end

        function [I] = int6d(f,x1_lhs,x1_rhs,y1_lhs,y1_rhs,z1_lhs,z1_rhs,...
                     x2_lhs,x2_rhs,y2_lhs,y2_rhs,z2_lhs,z2_rhs,N)
            a = zeros(1,6);
            b = zeros(1,6);
            
            a(1) = x1_lhs;
            b(1) = x1_rhs;
            [a(2),b(2)] = MonteCarlo.findMinMax1Arg(y1_lhs,y1_rhs,a(1),b(1));
            [a(3),b(3)] = MonteCarlo.findMinMax2Arg(z1_lhs,z1_rhs,a(1),b(1),a(2),b(2));
            
            a(4) = x2_lhs;
            b(4) = x2_rhs;
            [a(5),b(5)] = MonteCarlo.findMinMax1Arg(y2_lhs,y2_rhs,a(4),b(4));
            [a(6),b(6)] = MonteCarlo.findMinMax2Arg(z2_lhs,z2_rhs,a(4),b(4),a(5),b(5));
            
            f_sum = 0;
            
            for i = 1:N
                % generating random vector
                x = zeros(1,6);
                for j = 1:6
                    x(j) = MonteCarlo.randInRange(a(j),b(j));
                end
                % check conditions
                inArea = true;
                for j = 1:6
                    if (x < a(j) || x > b(j))
                        inArea = false;
                        break;
                    end
                end
                % adding f(x)
                if (inArea)
                    f_sum = f_sum + f(x);
                end
            end
            
            % what?
            V = 0;
            I = -1 * V;
            
      end
        
        function [I] = test1d(N)
          f = @(x) sqrt(7-3.*(sin(x)).^2);
          x_lhs = 0;
          x_rhs = 8;
          I = MonteCarlo.int1d(f,x_lhs,x_rhs,N);
        end

        function [I] = test2d(N)
          f = @(x,y) 1;
          x_lhs = 0;
          x_rhs = 1;
          y_lhs = @(x) x;
          y_rhs = @(x) x + 1;
          I = MonteCarlo.int2d(f,x_lhs,x_rhs,y_lhs,y_rhs,N);
        end
        
        function [I] = test3d(N)
            f = @(x,y,z) x^2;
            x_lhs = 0;
            x_rhs = 1;
            y_lhs = @(x) 0;
            y_rhs = @(x) 1 - x;
            z_lhs = @(x,y) 0;
            z_rhs = @(x,y) 10*(x+3*y);
            I = MonteCarlo.int3d(f,x_lhs,x_rhs,y_lhs,y_rhs,z_lhs,z_rhs,N);
        end
        
    end
    
end

