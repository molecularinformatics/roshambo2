from roshambo2 import Roshambo2



roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True)

for mode in [0,1,2]:
    scores = roshambo2_calculator.compute(backend='cuda', optim_mode='combination', reduce_over_conformers=True, start_mode=mode, write_scores=False) 

    print(scores)

    roshambo2_calculator.write_best_fit_structures(hits_sdf_prefix=f'hits_start_mode_{mode}')
