close all;
clear all;
clc;

%% PLOT CONVERGENCE ANALYSIS 
N_el = [46 81 92 138 184]; % Number of elements
nat_el = [6.7116 6.7054 6.7054 6.7050 6.7050
          8.2194 8.2141 8.2135 8.2132 8.2131 
          14.2797 14.2581 14.2564 14.2551 14.2549];% First 3 Natural Frequencies

figure(1) 
plot(N_el, nat_el(3,:),'o-','MarkerEdgeColor','r','MarkerSize',7,'LineWidth',1.5); hold on; grid on
xlabel('Number of Elements [/]');
ylabel('Natural Frequency [Hz]');
%legend('Monotonic \sigma - \epsilon curve','Cyclic curve','location','SouthEast');
%xlim([0 1200]);
ylim([14.25 14.285]);