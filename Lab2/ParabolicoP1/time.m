function [ t ] = time( M, t_0, t_max)

   t = zeros(1,M);
   t(1) = t_0;
   t_step = 1/M;
   
   for i=2:M
       t(i) = t(i-1) + t_step;
   end

end

