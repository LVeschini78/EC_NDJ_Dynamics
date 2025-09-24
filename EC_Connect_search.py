import EC_Connect_core
import nevergrad as ng
import numpy as np
import concurrent.futures
from time import time
import os

curr_wd = os.getcwd()
store = os.path.join(curr_wd, "Testing", "Scan_1.zarr")
"""
Define the parameter space for optimization
Nevergrad parameter space definition uses Instrumentation to group parameters
Each parameter is defined with its initial value and bounds
Nevergrad documentation: https://facebookresearch.github.io/nevergrad/
"""

param_space = ng.p.Instrumentation(
    cr = ng.p.Scalar(init=10).set_bounds(0,15),
    ps_bd4 = ng.p.Scalar(init=100).set_bounds(50,150),
    ps_bn1 = ng.p.Scalar(init=100).set_bounds(50,150),
    ps_bj1 = ng.p.Scalar(init=100).set_bounds(0,150),
    ps_ind = ng.p.Scalar(init=60).set_bounds(0,100), 
    ps_ja = 5000,
    ps_Kdni = 100,
    ps_ni = 200000,
    ps_kconv = ng.p.Scalar(init=10).set_bounds(0,15),
    ps_K_nui = 100,
)

"""
Define the optimizer. We choose Differential Evolution (DE)
https://facebookresearch.github.io/nevergrad/optimizers_ref.html#nevergrad.families.DifferentialEvolution
"""
optimizer = ng.optimizers.DE(param_space, budget=5000, num_workers=10)

"""
Objective function to minimize (or maximize)
This function runs the EC_Connect_core simulation with given parameters and returns 
the ks_stat value to be minimized
"""
def objective_function(**params):
    local_store = store
    # Generate a unique experiment name for each run to avoid overwriting
    local_expname = f"{np.random.randint(0, 100000)}"  # Randomly generated experiment name for each run
    # Initialize and run the simulation
    sim_core = EC_Connect_core.ECConnectCore(store=local_store,
                                             exp_name=local_expname,
                                             params=params)
    # Step the simulation for 1000 MCS and write results
    for mcs in range(1001):
        if mcs != 1000:
            sim_core.step()
        else:
            results = sim_core.step_write()
    ks_stat = results[0]
    # print(f"Finished simulation with params: {params}, ks_stat: {ks_stat}")
    return ks_stat

# Function to run optimization in parallel
def run_parallel_optimization():
    # Use ProcessPoolExecutor to run simulations in parallel 
    # https://docs.python.org/3/library/concurrent.futures.html
    with concurrent.futures.ProcessPoolExecutor(max_workers=optimizer.num_workers) as executor:
        futures = {}
        # Submit initial jobs
        for _ in range(optimizer.num_workers):
            param = optimizer.ask()
            # print(f"Submitting initial job with params: {param.kwargs}")
            future = executor.submit(objective_function, **param.kwargs)
            futures[future] = param
        # As jobs complete, submit new ones until budget is exhausted
        for i in range(optimizer.budget // optimizer.num_workers):
            done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                value = future.result()
                param = futures[future]
                # print(f"Job completed for params: {param.kwargs}, ks_stat: {value}")
                optimizer.tell(param, value)
                # Launch a new job
                new_param = optimizer.ask()
                # print(f"Submitting new job with params: {new_param.kwargs}")
                new_future = executor.submit(objective_function, **new_param.kwargs)
                futures[new_future] = new_param
                del futures[future]

        # Wait for any remaining jobs
        for future in futures:
            value = future.result()
            param = futures[future]
            # print(f"Final job completed for params: {param.kwargs}, ks_stat: {value}")
            optimizer.tell(param, value)

    recommendation = optimizer.provide_recommendation()
    print("Best parameters found:", recommendation.kwargs)
    print("Best objective value:", objective_function(**recommendation.kwargs))

# Define default parameters for simple parameter search
default_params = dict(cr = 10,
                      ps_bd4 = 100,
                      ps_bn1 = 100,
                      ps_bj1 = 0.0, 
                      ps_ind = 0.0, 
                      ps_ja = 5000, 
                      ps_Kdni = 100, 
                      ps_ni = 200000,
                      ps_kconv = 8,
                      ps_K_nui = 0.5)

# Simple function to run a single simulation with given parameters
def simple_sim_run(param_set:dict, expname:str):
    # We use the same store for simplicity, but in typical scenario we should consider unique stores
    # to avoid write conflicts
    store = store    
    # Initialize and run the simulation
    sim = EC_Connect_core.ECConnectCore(store=store,
                                        exp_name=expname,
                                        params=param_set)
    # Step the simulation for 1000 MCS and write results
    for mcs in range(1001):
        if mcs != 1000:
            sim.step()
        else:
            results = sim.step_write()
    print(f"Simulation {expname} completed. Ks Stat:", results[0])

# Function to a parameter search over two parameters with parallel execution
def run_two_params_search_parallel(param_1:str, range_1:tuple, param_2:str, range_2:tuple, n_points:int = 10, max_workers=10):
    # Create a grid of parameter values
    values_1 = np.round(np.linspace(range_1[0], range_1[1], n_points), 3)
    values_2 = np.round(np.linspace(range_2[0], range_2[1], n_points), 3)
    param_grid = [(v1, v2) for v1 in values_1 for v2 in values_2]

    # Run simulations in parallel using ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for idx, (val1, val2) in enumerate(param_grid):
            params = default_params.copy()
            params[param_1] = val1
            params[param_2] = val2
            expname = idx+1
            print(f"Submitting simulation {idx+1}/{len(param_grid)}: {expname}")
            futures.append(executor.submit(simple_sim_run, params, expname))
        # Wait for all to finish
        concurrent.futures.wait(futures)

"""
Main function, uncomment the desired function to run parameter optimization, simple 2-parameter searches, or
or individual simulations
"""
if __name__ == "__main__":
    start_time = time()
    # run_parallel_optimization()
    # results = simple_sim_run(default_params, 'test_default')
    run_two_params_search_parallel(param_1='ps_bj1', 
                                   range_1=(0.0, 150), 
                                   param_2='ps_ind', 
                                   range_2=(0.0, 100),
                                   n_points=30,
                                   max_workers=10)
    end_time = time()
    print(f"Optimization completed in {end_time - start_time:.2f} seconds.")