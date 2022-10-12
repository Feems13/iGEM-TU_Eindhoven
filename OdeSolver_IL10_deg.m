clear all;close all;
Defaults_deg;
P0 = paramIL10_deg();


options = odeset('RelTol',1e-9,'nonnegative',1);

%% Mainloop
% note:
% the upper for loop is parfor capable, if we need it. just know that
% parfor and plotting does not go well together. However, feel free to
% analyse the results from parfor, or save them for later plotting 

for i = 1:1 % vary one constant
    for j = 1:1 % vary another
        
        %ii = [i,j];
        ii=[]; % Keep ii empty to disable constant variation
        
        % solve the model
        [T,Y]=ode15s( @(t,y) ODEs_IL10_deg(t,y,P0,ii) ,T,Y0,options);
        
        [IL10_36] = features(T, Y(:,17));

        % any analasis here:
        % or define a function later to analyze the results with.
        % out = DoAnalasis(Y,); % or something like it.
        colors = [0.969 0.702 0.169; 0.0681 0.5556 0.3763; 0.2201 0.3042 0.4757; 0 0 0; 0.3588 0.2669 0.3744]; %STAT visualisation
        lines = {'-', '--', ':'};
        %Title = "\fontsize{14}\color{black}\bfIL-10 secretion";
        PlotYx(T,Y, [14:17],labels, colors, lines);
        xlim([0, 50]);
        legend('mRNA of IL-10 in cytoplasm', 'IL-10 protein of ER', 'IL-10 protein of Golgi', 'IL-10 protein extracellular');
        % ylim([0, 5000]);
%         titel1 = sprintf('IL-10 waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.jpg', P0(34), P0(3), P0(4));
%         titel2 = sprintf('IL-10 waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.fig', P0(34), P0(3), P0(4));
%         titel3 = sprintf('IL-10 waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.svg', P0(34), P0(3), P0(4));
%         saveas(gcf, titel1)
%         saveas(gcf, titel2)
%         saveas(gcf, titel3)
        
        colors = [0.969 0.702 0.169; 0.0681 0.5556 0.3763; 0.2201 0.3042 0.4757; 0 0 0; 0.3588 0.2669 0.3744]; %STAT visualisation
        %Title = "\fontsize{14}\color{black}\bfSTAT3 circle";
        PlotYx(T,Y, [1,7:9],labels, colors, lines); %,3,12,11,6,2,4
        xlim([0, 50]);
        legend('STAT3 in cytoplasm', 'Dimerized and phosphorylated STAT3 in cytoplasm', 'Dimerized and phosphorylated STAT3 in nucleus', 'Phosphorylated STAT3 in nucleus');
        % ylim([0, 1000]);
%         titel1 = sprintf('STAT3 waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.jpg', P0(34), P0(3), P0(4));
%         titel2 = sprintf('STAT3 waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.fig', P0(34), P0(3), P0(4));
%         titel3 = sprintf('STAT3 waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.svg', P0(34), P0(3), P0(4));
%         saveas(gcf, titel1)
%         saveas(gcf, titel2)
%         saveas(gcf, titel3)

        
        colors = [0.969 0.702 0.169; 0.0681 0.5556 0.3763; 0.2201 0.3042 0.4757; 0 0 0; 0.3588 0.2669 0.3744]; %STAT visualisation %receptor and ligand
        %Title = "\fontsize{14}\color{black}\bfReceptor species";
        PlotYx(T,Y, [19:22,2],labels, colors, lines);
        xlim([0 50]);
        ylim([0 18]);
        legend('Receptor', 'Receptor bound to ligand', 'Dimerized receptor', 'Dimerized and phosphorylated receptor', 'STAT3 bound to receptor');
%         titel1 = sprintf('Receptor waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.jpg', P0(34), P0(3), P0(4));
%         titel2 = sprintf('Receptor waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.fig', P0(34), P0(3), P0(4));
%         titel3 = sprintf('Receptor waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.svg', P0(34), P0(3), P0(4));
%         saveas(gcf, titel1)
%         saveas(gcf, titel2)
%         saveas(gcf, titel3)
        
        colors = [0.969 0.702 0.169; 0.0681 0.5556 0.3763; 0.2201 0.3042 0.4757; 0 0 0; 0.3588 0.2669 0.3744]; %STAT visualisation
        %Title = "\fontsize{14}\color{black}\bfConcentration of STAT3cp over time";
        PlotYx(T,Y, [3],labels, colors, lines);
        xlim([0, 50]);
        legend('Phosphorylated STAT3 in the cytoplasm');
%         titel1 = sprintf('STAT3cp waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.jpg', P0(34), P0(3), P0(4));
%         titel2 = sprintf('STAT3cp waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.fig', P0(34), P0(3), P0(4));
%         titel3 = sprintf('STAT3cp waarden IL-10, KDc1 %.4f KM3 %.3f KM4 %.0f.svg', P0(34), P0(3), P0(4));
%         saveas(gcf, titel1)
%         saveas(gcf, titel2)
%         saveas(gcf, titel3)
    end
end

%% Plot
function PlotYx(t,Y,x,labels,colors,lines)
figure; hold on;
colors = [0.969 0.702 0.169; 0.0681 0.5556 0.3763; 0.2201 0.3042 0.4757; 207/255 154/255 216/255; 0 0 0]; %STAT visualisation
set(gcf, 'color', 'w'); 
set(gca, 'ColorOrder', colors, 'LineStyleOrder', lines');
%reset(gca);
plot(t/3600, Y(:,x), 'Linewidth', 1.5);
%ylim([0,200]);
xlabel('Time [hour]')
ylabel('Concentration [nM]')
%legend(labels(x))

% Create a title and a subtitle with the initial conditions
subtitle = '\fontsize{10}\color{gray}\rmInitial [nM]: ';
Y0 = Y(1,:);
for n = 1:length(Y0)
    if Y0(n) ~= 0
        subtitle = strcat(subtitle, labels(n), '=', num2str(Y0(n)), ',');
    end
end
subtitle = string(subtitle); 
Title = '';
title({Title; subtitle});
end