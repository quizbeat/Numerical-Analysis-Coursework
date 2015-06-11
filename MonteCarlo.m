classdef MonteCarlo
    % ?????????? ??????? ?????????? ??????? ????? ?????
    %   ??? ?????-?? ????
    
    methods(Static)
        
        function [x] = randInRange(a,b)
        % returns random value in range [a,b]
            x = a + (b - a)*rand();
        end
        
        function [I] = int1d(f,a,b,N)
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
               if (y <= f(x)) % point under function plot
                   plot(x,y,'--.b');
                   n = n + 1;
               else
                   plot(x,y,'--.r');
               end
            end

            I = (b-a) * c * n / N;
        end
        
        function [I] = int3d(f,x_lhs,x_rhs,y_lhs,y_rhs,z_lhs,z_rhs,N)
            % ???f(x,y,z)dxdydz = I
            % triple definite integral with conditions:
            % x_lhs ? x ? x_rhs
            % y_lhs ? y ? y_rhs
            % z_lhs ? z ? z_rhs
            step = 0.5;
            
            min_x = x_lhs();
            max_x = x_rhs();
            
            min_y = inf;
            max_y = -inf;
            % find min_y, max_y
            for x = min_x:step:max_x
               if (y_lhs(x) < min_y)
                   min_y = y_lhs(x);
               end
               if (y_rhs(x) > max_y)
                   max_y = y_rhs(x);
               end
            end
            
            min_z = inf;
            max_z = -inf;
            % find min_z, max_z
            for x = min_x:step:max_x
                for y = min_y:step:max_y
                    if (z_lhs(x,y) < min_z)
                        min_z = z_lhs(x,y);
                    end
                    if (z_rhs(x,y) > max_z)
                        max_z = z_rhs(x,y);
                    end
                end
            end
            
            max_z = 30;
                
            n = 0;
            sum = 0;
            
            for i = 1:N
                x = MonteCarlo.randInRange(min_x,max_x);
                y = MonteCarlo.randInRange(min_y,max_y);
                z = MonteCarlo.randInRange(min_z,max_z);
                if (x_lhs() <= x & x <= x_rhs() & ...
                    y_lhs(x) <= y & y <= y_rhs(x) & ...
                    z_lhs(x,y) <= z & z <= z_rhs(x,y))
                    n = n + 1;
                    sum = sum + f(x,y,z);
                end
            end
            
            V = (max_x-min_x) * (max_y-min_y) * (max_z-min_z);
            I = V / N * sum;
        end
        
    end
    
end

