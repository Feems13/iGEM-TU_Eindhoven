clear all; close all;

Defaults_deg;

P0 = paramIL10_deg;

options = odeset('RelTol', 1e-9, 'nonnegative', 1);

RJ = 24;
Rarray = 10.^[-2:0.1:1.5]*RJ;
IL10 = zeros(5, length(Rarray));
RJ2Lp = zeros(5, length(Rarray));

for i = 1:length(Rarray)
    ii = [1008, i];
    [T,Y] = ode15s( @(t,y) ODEs_IL10_deg(t,y,P0,ii), T, Y0, options);
    [IL10_1, IL10_10, IL10_20, IL10_30, IL10_36, RJ2Lp_1, RJ2Lp_10, RJ2Lp_20, RJ2Lp_30, RJ2Lp_36] = features_RJ(T, Y(:,17), Y(:,22));
    
    IL10(:,i) = [IL10_1, IL10_10, IL10_20, IL10_30, IL10_36];
    RJ2Lp(:,i) = [RJ2Lp_1, RJ2Lp_10, RJ2Lp_20, RJ2Lp_30, RJ2Lp_36];
end

L = 0.24; %[nM]
ratio = Rarray/L;

names = {'t = 1h', 't = 10h', 't = 20h', 't = 30h', 't = 36h'};

figure()
set(gcf, 'color', 'w');
semilogx(ratio, RJ2Lp(5,:), '*', 'MarkerSize', 7, 'MarkerEdgeColor', [0.969 0.702 0.169]);
hold on;
h = stem([5.0119, 501.1872],[0.1111,11.4939], '--*', 'MarkerSize', 7, 'LineWidth', 1.5); % ]
h.Color =  [0.0681 0.5556 0.3763];
%stem(794.3282, 5.6012, ':*');
%title(names(5));
xlabel('Ratio Receptor: Ligand [-]', 'FontSize', 14);
ylabel('Dimerized and phosphorylated receptor concentration [nM]', 'FontSize', 14);
xlim([0.9, 10^3.5]);