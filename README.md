# jfactors

Code associated with the two papers: Evans, Sanders & Geringer-Sameth (2016) (in [spherical](spherical/)) and Sanders, Evans,  Geringer-Sameth & Dehnen (2016) (in [flattened](flattened/)).

## Spherical (Evans, Sanders & Geringer-Sameth 2016)

Associated scripts in [spherical](spherical/). The figures in the paper are generated by the following scripts:

1. [against_distance.py](spherical/against_distance.py) -- Figure 1
2. [J_D_table.py](spherical/J_D_table.py) -- Table 1
3. [J_D_profiles.py](spherical/J_D_profiles.py) -- Figures 2, 3, 4, 5, 6, 7
4. [sweet_spot.py](spherical/sweet_spot.py) -- Figure 8

## Flattened (Sanders, Evans,  Geringer-Sameth & Dehnen 2016)

Associated scripts in [flattened](flattened/) and also [C++ code](c++_code/) . The following scripts generate the figures and tables in the paper

1. [fig1_flattened.py](flattened/fig1_flattened.py) -- Figure 1
2. [generate_ret2_table.py](flattened/generate_ret2_table.py) -- Table 1
3. [compute_J_ret2_sims.py](flattened/compute_J_ret2_sims.py) -> [sim_profiles.py](flattened/sim_profiles.py) -- Figures 2 & 3
4. [flattened.py](flattened/flattened.py) -> Figures 5,6 & 7 & Table 3 (need Wyn's data and output from [C++ code](c++_code/))
5. [compute_corrections.py](flattened/compute_corrections.py) -- Table 3 & 4
6. [corr_plot.py](flattened/corr_plot.py) -- Figure 4
7. [unit_sphere.py](flattened/unit_sphere.py) -- Figure 8
8. [sampler.py](flattened/sampler.py) -> [process_data.py](flattened/process_data.py) -- Table 5
9. [sampler.py](flattened/sampler.py) -> [ret2_hists.py](flattened/ret2_hists.py) and [all_plot.py](flattened/all_plot.py) -- Figures 9 & 10

The manuscript is in [paper](paper/) along with the data for the tables and plots.

## AUTHORS

Jason Sanders - jls at ast dot cam dot ac dot uk

Wyn Evans

Alex Geringer-Sameth


