close all; clear all;
Defaults_deg;

options = odeset('RelTol', 1e-9, 'nonnegative',1);

Kd1 = 0.741; % [nM]
KDarray = 10.^[-3:0.1:3]*Kd1;
P0 = paramIL10_deg;
IL10 = zeros(9,length(KDarray));
IL10_echt = zeros(4,length(KDarray));

for i = 1:length(KDarray)
    ii = [1005, i]; 
    [T,Y] = ode15s( @(t,y) ODEs_IL10_deg(t,y,P0,ii), T, Y0, options);
    [IL10_1, IL10_2, IL10_5, IL10_10, IL10_15, IL10_20, IL10_25, IL10_30, IL10_36]= features_kd(T, Y(:,17));
    
    IL10(:,i) =  [IL10_1, IL10_2, IL10_5, IL10_10, IL10_15, IL10_20, IL10_25, IL10_30, IL10_36];
    IL10_echt(:,i) = [IL10_1, IL10_2, IL10_5, IL10_10];
end

% figure(1)
% semilogx(KDarray, IL10(1,:), 'o')
% figure(5)
% semilogx(KDarray, IL10(2,:), 'o')
% set(gcf, 'color', 'w');

%% Plotjes
figure(1)
names = {'t = 1h', 't = 2h', 't = 5h', 't = 10h'};


for j = 1:size(IL10_echt,1)
    subplot(2,2,j)
    set(gcf, 'Color', 'w')
    semilogx(KDarray, IL10_echt(j,:), '*', 'MarkerEdgeColor', [0.969 0.702 0.169])  
    title(names(j))
    xlabel('K_D1 [nM]', 'FontSize', 14)
    ylabel('Concentration IL-10 [nM]', 'FontSize', 14)
    xlim([10^(-3), 10^3])
end 
subplot(2,2,1)
ylim([0,0.9])
subplot(2,2,2)
ylim([0,6])
subplot(2,2,3)
ylim([0,17])
subplot(2,2,4)
ylim([0,17])

figure(2);
set(gcf, 'Color', 'w');
semilogx(KDarray, IL10(9,:), '*', 'MarkerEdgeColor', [0.969 0.702 0.169]);
xlabel('K_D1 [nM]', 'FontSize', 14);
ylabel('Concentration IL-10 [nM]', 'FontSize', 14);
xlim([10^(-3), 10^3])
ylim([0,17])