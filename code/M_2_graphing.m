 %% Set up and load results 
cd('/Volumes/GoogleDrive/My Drive/Github/Urban-Frictions-Sim/code');
addpath(genpath('depends'));
friction_results= load('../output/friction_results.mat').friction_results;
no_friction_results= load('../output/no_friction_results.mat').no_friction_results;

%% Now create the graphs (without frictions)
dimensions=size(no_friction_results);
length = dimensions(2)-2;
red = [1, 0, 0];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),length)', linspace(red(2),pink(2),length)', linspace(red(3),pink(3),length)'];


%distance graph
 scatter(no_friction_results(:,1), no_friction_results(:,3),10, colors_p(1,:))
  xlabel('Log Distance')
  ylabel('Log Welfare (normalized to zero)')
  hold on
 
  for i =2:length
     scatter(no_friction_results(:,1),no_friction_results(:,i+2),10, colors_p(i,:));
     
  end
  legend('10%','9%', '8%', '7%', '6%', '5%', '4%', '3%', '2%','1%');
    hold off 
  f = gcf;
  exportgraphics(f,'../figures/link_distance_no_frictions.png','Resolution',300)
  
 %popularity graph
  scatter(no_friction_results(:,2), no_friction_results(:,3),10, colors_p(1,:))
  xlabel('Negative Log popularity')
  ylabel('Log Welfare (normalized to zero)')
  hold on
 
  for i =2:length
     scatter(no_friction_results(:,2),no_friction_results(:,i+2),10, colors_p(i,:));
   
  end
    hold off 
  f = gcf;
  exportgraphics(f,'../figures/popularity_no_frictions.png','Resolution',300)
 
 %% Create the Friction graphs
dimensions=size(friction_results);
length = dimensions(2)-2;
red = [1, 0, 0];
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),length)', linspace(red(2),pink(2),length)', linspace(red(3),pink(3),length)'];


%distance graph
 scatter(friction_results(:,1), friction_results(:,3),10, colors_p(1,:))
  xlabel('Log Distance')
  ylabel('Log Welfare (normalized to zero)')
  hold on
 
 % for i =2:length
  %   scatter(friction_results(:,1),friction_results(:,i+2),10, colors_p(i,:));
 % end
  %legend('0% real ','25% real', '50% real', '75 real', '100% real');
    hold off 
  f = gcf;
  exportgraphics(f,'../figures/link_distance_frictions.png','Resolution',300)
  %%
 %popularity graph
  scatter(friction_results(:,2), friction_results(:,3),10, colors_p(1,:))
  xlabel('Negative Log popularity')
  ylabel('Log Welfare (normalized to zero)')
  hold on
 
  for i =2:length
     scatter(friction_results(:,2),friction_results(:,i+2),10, colors_p(i,:));
   
  end
    hold off 
  f = gcf;
  exportgraphics(f,'../figures/popularity_frictions.png','Resolution',300)
