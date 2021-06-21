function varproj_figgen_obj(varargin)

if mod(nargin,2)~=0
    disp('Function Use: varproj_figgen_obj(data, name)')
end
cv = {'-b','-k', '-r', '--b', '--k', '--r', ':b', ':k', ':r', '-.b', '-.k', '-.r'};
leg = cell(nargin/2,1); 
figure;    
xlabel('Iteration (log-scale)', 'Fontsize', 22, 'Fontweight', 'Bold')
ylabel('||X - X_{true}||_2 (log-scale)', 'Fontsize', 22, 'Fontweight', 'Bold')
set(gca, 'Fontsize', 22, 'Fontweight', 'Bold')
title('L_2 Norm: True - Interpolated', 'Fontsize', 22, 'Fontweight', 'Bold')

hold on

for i = 1:nargin/2
    obj_value = varargin{2*i-1};
    plot(log(1:numel(obj_value)),log(obj_value/max(obj_value)), cv{i}, 'LineWidth', 2); 
    leg{i} = varargin{2*i}; 
end
legend(leg)




end