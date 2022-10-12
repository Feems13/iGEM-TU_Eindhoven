clear all;close all;
Defaults
p_names = labels_parameters([12,13, 20,21, 31,32]); %labels_parameters([12, 13, 20, 21, 31, 32]);
units = units_parameters;
feature_label = {'IL-10 concentration at t = 36 h'};
%%
load('results_MPSA_n10000_p3_28-Sep-2022.mat')
load('init_MPSA_n10000_p3_28-Sep-2022.mat')
rng(10,'twister') 
%p_tweak_values_dummies(:,1:3) = p_tweak_values;

%%
V = V(:,1);

% Classification of acceptable vs unacceptable simulations
flag  = zeros(n_MPSA, 1);
Sa    = zeros(n_MPSA, Np);
Su    = zeros(n_MPSA, Np);
value = zeros(n_MPSA, Np);

threshold = mean(V);    %mean as threshold or median 
acc       = find(V <= threshold);
unacc     = find(V  > threshold);
flag(acc) = 1;

% Cumulative distributions (for model parameters and dummies)
for i = 1 : Np+Nd
    temp = [p_tweak_values_dummies(:,i), flag];     %associate 1 to acceptable cases and 0 to unacceptable parameter values
    temp = sortrows(temp,1);        %sorts temp based on column 1

    value(:,i) = temp(:,1);
    Sa(:,i) = cumsum(temp(:,end));
    Su(:,i) = cumsum(-1*temp(:,end)+ones(n_MPSA,1));
    Sa(:,i) = Sa(:,i)/max(Sa(:,i));
    Su(:,i) = Su(:,i)/max(Su(:,i)); 
end

figure()
for i = 1:Np
    subplot(2,3,i)
    plot(value(:,i),Sa(:,i))
    hold on
    plot(value(:,i),Su(:,i))
    title(p_names(i))
end
 %sgtitle(sprintf('MPSA %s',feature_labels{Vi}))
%     saveas(gcf,sprintf('MPSA_all_p_%s.png',t))
 hold off
 
% Kolmogorov-Smirnov statistic
K_S = max(abs(Sa-Su)); %Higher K_S means more sensitive
ratio_groups = length(unacc)/length(acc);

for ksi = 1:Np
    if K_S(ksi) <= max(K_S(Np+1:Np+Nd))
        K_S_scaled(ksi) = 0;
    else
        K_S_scaled(ksi) = K_S(ksi)/max(K_S(Np+1:Np+Nd));
    end
end

%%Plots
figure(10)
plot(K_S(1:Np),'o')
title('K_S scores')

figure(20)
plot(K_S(Np+1:Np+Nd),'o')
title('K_S scores dummies')

figure(4)
d_mean = mean(K_S(Np+1:Np+Nd),2);
d_std = std(K_S(Np+1:Np+Nd),0,2);
errorbar(d_mean,d_std,'o')


%Final plot for report
figure(5)
set(gcf, 'Color', 'w');
h = heatmap(K_S(1:Np));
%h.XData = p_names;
%h.CellLabelFormat = '%1.2f';
%h.YData = feature_label;
%h.title('Multiple parameter sensitivity analysis')
h.Colormap = parula;
s = struct(h);
s.XAxis.TickLabelRotation = 90;
s.YAxis.TickLabelRotation = 90;
set(gca, 'FontSize', 18);