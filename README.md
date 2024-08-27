# Event-based-flood-risk
These scripts are used to generate the event-based probalistic framework to capture flood spatial dependence for assessing coastal flood risk at the global scale. For details, see the manuscript _Li, H., Eilander, D., Haer, T., & Ward, P. J. (2024). Improving global-scale coastal risk assessments by considering spatial dependence. In review._ A preprint doi will be soon provided.

**Input dataset**
1. The improved version of the global dataset of spatially dependent extreme sea level (esl) events (Li et al., 2024), extended from Li et al. (2023). The new dataset is freely accessed at: https://doi.org/10.5281/zenodo.12722313
2. The GLOFRIS coastal inundation maps from Mortensen et al. (2024), which is publicly available: https://doi.org/10.5281/zenodo.10637089

**Scripts**
1. Hazard

inundation.py-->constructing inundation maps for synthetic events using the esl dataset and the GLOFRIS inundation maps
source_station.py-->developing a specific coastal segment for each station (or assigning the inundated cells to their sourcing station)
source_station_basin.py-->developing a source station map for the entire basin
merge_source_station_globe.py-->merging basin source station maps into a global source station map

3. Impact
linking_damage_to_location.py-->making damage maps for each location

dmg_individual_station.py-->calculating the damage cells caused by a given station (i.e. for its coastal segment)

damage_percentage.py-->calculating the damage percentage for segments which are in multiple subnational basins

add_damage_perc.py-->spliting the damages for segments which are in mupltiple subnational basins based on the damage percentages

agg_dmg_annual_damage.py-->aggregating annual damages across different basins

5. Plotting (for producing figures in the manuscript Li et al., 2024)
risk_curve_plotting_bootstrap.py-->plotting national aggregated risk curves

ead_map.py-->plotting a global map showing the expected annual damage (EAD) esimate difference

dmg_rp200_map.py-->plotting a global map showing the RP200 damage esimate difference

coastline_length_diff_plot.py-->plotting the risk estimate differences with different coastline lengths

worst_year_plotting.py-->plotting continetal maps showing the flood damages in the year with the highest combined annual damages


**References**
Li, H., Eilander, D., Haer, T., & Ward, P. J. (2024). Improving global-scale coastal risk assessments by considering spatial dependence. In review.
Li, H., Haer, T., Couasnon, A., Enríquez, A. R., Muis, S., & Ward, P. J. (2023). A spatially-dependent synthetic global dataset of extreme sea level events. Weather and Climate Extremes, 41, 100596. https://doi.org/10.1016/j.wace.2023.100596 
Mortensen, E., Tiggeloven, T., Haer, T., van Bemmel, B., Le Bars, D., Muis, S., Eilander, D., Sperna Weiland, F., Bouwman, A., Ligtvoet, W., & Ward, P. J. (2024). The potential of global coastal flood risk reduction using 690 various DRR measures. Natural Hazards and Earth System Sciences, 24(4), 1381–1400. https://doi.org/10.5194/nhess-24-1381-2024 
