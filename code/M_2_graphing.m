 %% Set up and load results 
cd('/Volumes/GoogleDrive/My Drive/ grad school/Research/RA work/Gabriel /Nick stuff /old_model_simulation/code');
addpath(genpath('depends'));

 results= load('../output/results.mat').results;
 results_graphable= [results (results(:,3)-results_init_eqbm.U_bar)/results_init_eqbm.U_bar*100 (results(:,4)-results_init_eqbm.U_bar)/results_init_eqbm.U_bar (results(:,5)-results_init_eqbm.U_bar)/results_init_eqbm.U_bar];
 
  %% Now create the graphs 
 scatter(results_graphable(:,1), results_graphable(:,6))
 title('Welfare Variation based on original link distance')
  xlabel('Original Link Distance')
  ylabel('percentage change in welfare')
  hold on
  scatter(results_graphable(:,1),results_graphable(:,7));
  scatter(results_graphable(:,1),results_graphable(:,8));
  legend('10% schock','5% shock','1% shock');
  hold off 
  f = gcf;
  exportgraphics(f,'../output/figures/link_distance_no_agglomeration.png','Resolution',300)

 
 scatter(results_graphable(:,2),results_graphable(:,6))
 title('Welfare Variation based on original link popularity')
  xlabel('Original Link Popularity')
  ylabel('percentage change in welfare')
  hold on
  scatter(results_graphable(:,2),results_graphable(:,7));
  scatter(results_graphable(:,2),results_graphable(:,8));
  legend('10% schock','5% shock','1% shock');
  hold off 
  f = gcf;
  exportgraphics(f,'../output/figures/link_popularity_no_agglomeration.png','Resolution',300)

  