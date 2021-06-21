function varproj_figgen_svd(varargin)

if mod(nargin,2)~=0
    disp('Function Use: varproj_figgen_svd(data, name)')
end
cv = {'-b','-k', '-r', '--b', '--k', '--r', ':b', ':k', ':r', '-.b', '-.k', '-.r'};
leg = cell(nargin/2,1); 
figure;    
ylabel('Relative Singular Value \sigma_i/\sigma_{max}', 'Fontsize', 22, 'Fontweight', 'Bold')
xlabel('ith index', 'Fontsize', 22, 'Fontweight', 'Bold')
set(gca, 'Fontsize', 22, 'Fontweight', 'Bold')
title('SV Decay: Residuals', 'Fontsize', 22, 'Fontweight', 'Bold')

hold on

for i = 1:nargin/2
    temp = svd(varargin{2*i-1});
    plot(1:numel(temp),temp/max(temp), cv{i}, 'LineWidth', 2); 
    leg{i} = varargin{2*i}; 
end
legend(leg)

end