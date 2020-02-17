function [xtot,Measurement,  Prediction] = fitLDS_simulateddata(drmat0,  rmat0, Ntrials, rngtot)

   Nstims = 2;

    rtot = horzcat(rmat0{:});
    stimblocks  = [ones([1, length(rngtot) * Ntrials{1}, 1, 1]),zeros([1, length(rngtot) * Ntrials{2}, 1, 1]); ...
                                             zeros([1, length(rngtot) *  Ntrials{1}, 1, 1]), ones([1, length(rngtot) * Ntrials{2}, 1, 1])];
    
     Lambdatot_full = [rtot;stimblocks];  % only allow weights to change                                     
    


drmat0tot = horzcat(drmat0{:});
           
              

      
xtot = drmat0tot / Lambdatot_full;
Prediction_tot = xtot * Lambdatot_full;


Measurement = horzcat(drmat0{:});
Prediction = Prediction_tot(:, 1:size(Measurement,2));

    
end