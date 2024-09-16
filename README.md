# Event-based-flood-risk
These scripts are used to generate the event-based probalistic framework to capture flood spatial dependence for assessing coastal flood risk at the global scale. For details, see the manuscript _Li, H., Eilander, D., Haer, T., & Ward, P. J. (2024). Improving global-scale coastal risk assessments by considering spatial dependence. In review. The preprint can be accessed at: 10.22541/essoar.172641608.83190937/v1. 

**Input dataset**
1. The improved version of the global dataset of spatially dependent extreme sea level (esl) events (Li et al., 2024), extended from Li et al. (2023). The new dataset is freely accessed at: https://doi.org/10.5281/zenodo.12722313
2. The GLOFRIS coastal inundation maps from Mortensen et al. (2024), which is publicly available: https://doi.org/10.5281/zenodo.10637089

**Scripts**

Hazard
1. inundation.py-->constructing inundation maps for synthetic events using the esl dataset and the GLOFRIS inundation maps
2. source_station.py-->developing a specific coastal segment for each station (or assigning the inundated cells to their sourcing station)
3. source_station_basin.py-->developing a source station map for the entire basin
4. merge_source_station_globe.py-->merging basin source station maps into a global source station map

Impact
1. linking_damage_to_location.py-->making damage maps for each location
2. dmg_individual_station.py-->calculating the damage cells caused by a given station (i.e. for its coastal segment)
3. damage_percentage.py-->calculating the damage percentage for segments which are in multiple subnational basins
4. add_damage_perc.py-->spliting the damages for segments which are in mupltiple subnational basins based on the damage percentages
5. agg_dmg_annual_damage.py-->aggregating annual damages across different basins

Plotting (for producing figures in the manuscript Li et al., 2024)
1. risk_curve_plotting_bootstrap.py-->plotting national aggregated risk curves
2. ead_map.py-->plotting a global map showing the expected annual damage (EAD) esimate difference
3. dmg_rp200_map.py-->plotting a global map showing the RP200 damage esimate difference
4. coastline_length_diff_plot.py-->plotting the risk estimate differences with different coastline lengths
5. worst_year_plotting.py-->plotting continetal maps showing the flood damages in the year with the highest combined annual damages

**References**
1. Li, H., Eilander, D., Haer, T., & Ward, P. J. (2024). Improving global-scale coastal risk assessments by considering spatial dependence. In review.
2. Li, H., Haer, T., Couasnon, A., Enríquez, A. R., Muis, S., & Ward, P. J. (2023). A spatially-dependent synthetic global dataset of extreme sea level events. Weather and Climate Extremes, 41, 100596. https://doi.org/10.1016/j.wace.2023.100596
3. Mortensen, E., Tiggeloven, T., Haer, T., van Bemmel, B., Le Bars, D., Muis, S., Eilander, D., Sperna Weiland, F., Bouwman, A., Ligtvoet, W., & Ward, P. J. (2024). The potential of global coastal flood risk reduction using 690 various DRR measures. Natural Hazards and Earth System Sciences, 24(4), 1381–1400. https://doi.org/10.5194/nhess-24-1381-2024 
