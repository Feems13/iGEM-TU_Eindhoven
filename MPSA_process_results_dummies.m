
%p_names = p_names_tex();
%units = p_units_tex();
%BAs = [{'CA'}, {'CDCA'}, {'DCA'}, {'UDCA'}, {'LCA'},{'LCAs'},{'Other'}];
%states = [{'t'}, {'g'}, {'u'}];

%l = [strcat(BAs,' Peripheral'),strcat(states,' Peripheral'),strcat(BAs,' Portal'),strcat(states,' Portal'),strcat(BAs,' Fecal'),strcat(states,' Fecal')];

%%
clear all;close all;

load('results_MPSA_n10000_p3_11-Oct-2022.mat')
load('init_MPSA_n10000_p3_11-Oct-2022.mat')
Vinit = V;
clear V
N= linspace(log(500),log(10000),40);
N= round(exp(N));
s = RandStream('mlfg6331_64');
rng(10,'twister') 

p_tweak_values_dummies(:,1:3) = p_tweak_values;

figure()
for Ni = 1:length(N)
        %%
    [Vtot,i_random] = datasample(s,Vinit,N(Ni),1,'Replace',false);
    par_MC_i = p_tweak_values_dummies(i_random,:);
    n_MPSA = N(Ni);
    K_S = zeros(size(Vtot,2),Np+Nd);
    K_S_scaled = NaN(size(Vtot,2),Np);
    ratio_groups = zeros(size(Vtot,2),1);

    V = Vtot(:,1);
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
        temp = [par_MC_i(:,i), flag];     %associate 1 to acceptable cases and 0 to unacceptable parameter values
        temp = sortrows(temp,1);        %sorts temp based on column 1

        value(:,i) = temp(:,1);
        Sa(:,i) = cumsum(temp(:,end));
        Su(:,i) = cumsum(-1*temp(:,end)+ones(n_MPSA,1));
        Sa(:,i) = Sa(:,i)/max(Sa(:,i));
        Su(:,i) = Su(:,i)/max(Su(:,i)); 
    end
        
        %% Plot in report of different sensitivities
%         if Vi == 28
%             subplot(1,3,1)
%             plot(value(:,19),Sa(:,19),'LineWidth',2)
%             hold on
%             plot(value(:,19),Su(:,19),'LineWidth',2)
%             title([l(Vi)])
%             xlabel(sprintf('%s [%s]',p_names{19},units{19}))
%             hold off
%         elseif Vi == 11
%             subplot(1,3,2)
%             plot(value(:,24),Sa(:,24),'LineWidth',2)
%             hold on
%             plot(value(:,24),Su(:,24),'LineWidth',2)
%             title([l(Vi)])
%             xlabel(sprintf('%s [%s]',p_names{24},units{24}))
%             hold off
%         elseif Vi == 2
%             subplot(1,3,3)
%             plot(value(:,25),Sa(:,25),'LineWidth',2)
%             hold on
%             plot(value(:,25),Su(:,25),'LineWidth',2)
%             title([l(Vi)])
%             xlabel(sprintf('%s [%s]',p_names{25},units{25}))
%             hold off
%         end
        %%
%         for i = 1:Np
%             subplot(6,6,i)
%             plot(value(:,i),Sa(:,i))
%             hold on
%             plot(value(:,i),Su(:,i))
%             title(p_names(i))
%         end
    %     sgtitle(sprintf('MPSA %s',t))
    %     saveas(gcf,sprintf('MPSA_all_p_%s.png',t))
    %     hold off
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

    % figure()
    % plot(K_S(Vp,1:Np),'o')
    % title(t)
    % figure()
    % plot(K_S(Vp,Np+1:Np+Nd),'o')
    % title(t)

    d_mean(:,Ni) = mean(K_S(:,Np+1:Np+Nd),2);
    d_std(:,Ni) = std(K_S(:,Np+1:Np+Nd),0,2);
    %errorbar(d_mean,d_std,'o')
end

%% Dummy plot in report of DCA feature
for i=1%1:size(Vinit,2)
    figure()
    set(gcf, 'color', 'w');
    %subplot(1,1,i)
    e = errorbar(N([1,3:end]),d_mean(i,[1,3:end]),d_std(i,[1,3:end]),'s',...
        'LineWidth',2,'MarkerSize',6, 'Color', [89/255 185/255 150/255],...
        'MarkerEdgeColor', [19/255 155/255 105/255], ...
        'MarkerFaceColor', [19/255 155/255 105/255]);
    %hold on;
    ylabel('Sensitivity ', 'FontSize', 16)
    xlabel('Sample size (N)', 'FontSize', 16)
    %title(['MPSA dummies';l(i)])
    h = legend('Mean + std');
    set(h,'FontSize',12);
    ylim([0,0.15]);
    xlim([0 10100]);
end

%% MPSA dummy parameters plot appendix
% figure()
% 
% for i=1:size(Vinit,2)
%     subplot(5,6,i)
%     errorbar(N([1,3:end]),d_mean(i,[1,3:end]),d_std(i,[1,3:end]),'ko','LineWidth',1,'MarkerSize',1)
%     %hold on;
%     if ismember(i,[1,7,13,19,25])
%         ylabel('Sensitivity')
%     end
%     if ismember(i,[25:30])  
%         xlabel('Sample size')
%     end
%     title([l(i)])
% end
% legend('Mean \pm standard deviation')
% %suptitle('MPSA dummies')

%% Final plot for report
feature_l = [1:size(K_S,1)]; %[1,3,9,11,13,19,21,23,29];

figure()
h = heatmap(K_S(feature_l,1:Np));
h.XData = p_names;
%h.CellLabelFormat = '%1.2f';
h.YData = l(feature_l);
%h.title('Multiple parameter sensitivity analysis')
h.Colormap = parula;
s = struct(h);
s.XAxis.TickLabelRotation = 90;
