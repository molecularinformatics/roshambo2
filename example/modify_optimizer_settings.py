from roshambo2 import Roshambo2
 

roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=False, verbosity=2)


# to change the optimizer settings create a dict like below
# these are the default settings I found to work well
optimizer_settings={'lr_q':0.1, 'lr_t':0.1, 'steps':100}

scores = roshambo2_calculator.compute(backend='cuda', optim_mode='shape', reduce_over_conformers=False, start_mode=1, optimizer_settings=optimizer_settings) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)

# we can run with more steps

optimizer_settings={'lr_q':0.1, 'lr_t':0.1, 'steps':200}

scores = roshambo2_calculator.compute(backend='cuda', optim_mode='shape', reduce_over_conformers=False, start_mode=1, optimizer_settings=optimizer_settings) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)

# we can run with less steps
optimizer_settings={'lr_q':0.1, 'lr_t':0.1, 'steps':50}

scores = roshambo2_calculator.compute(backend='cuda', optim_mode='shape', reduce_over_conformers=False, start_mode=1, optimizer_settings=optimizer_settings) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)

# the lr's can be changed but I have not had much success changing them
