using Revise
using ElectroPhysiology
using PhysiologyModeling

import PhysiologyModeling.calculate_dendrogram_distance
#%%=[Run branch generation]__________________________________________________________________________________#
radial = 5
branches = 2
layers = 5

xs1, ys1, connection_list1 = create_dendrogram_map(radial, branches, layers; origin = (0.0, 0.0))
xs2, ys2, connection_list1 = create_dendrogram_map(radial, branches, layers; origin = (0.0, 2.0))
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))