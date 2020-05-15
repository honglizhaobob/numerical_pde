hmax = 0.15;
pv_square = [0,0; 1,0; 1,1; 0,1; 0,0];
pv_polygon = [0,0; 1,0; .5,.5; 1,1; 0,1; 0,0];

for pv = {pv_square, pv_polygon} 
  errors = poiconv(pv{1}, hmax, 3)
  hold on
  grid on
  figure(1)
  loglog(hmax ./ 2.^(0:2), errors, 'LineWidth', 2)
  title('error convergence rate')
  xlabel('refinement'); ylabel('log-log error'); 
  disp('displaying the convergence rate ...')
  disp('====================')
  
  rate = log2(errors(end-1)) - log2(errors(end))
end
legend('square mesh', 'polygon mesh')
