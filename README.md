# CATS-H0

Based on "sn_analysis_new.csv", there is a two-step detection to give TRGB of each field.

We have different spatial clipping extents, color band widths, smoothing scales, minimum relative thresholds of peaks and internal extinction corrections. For each of the combinations, there is a detection file providing the TRGB magnitudes, and an information file providing essential statistics of all the fields.

We also provide the initial csv files in "phot_csv" folder which pass through EDD photometry selection (Tully+09, Anand+21b), clipped csv files in "clipped_csv" folder, files of stars removed after spatial clipping in "cut_csv" folder, and spatial clipping maps in "clipping_map" folder.

To reproduce our results, just run "Comprehensive_part1.py" and "Comprehensive_part2.py" in turn.

see EDD Phot in
1. Tully, R. B., Rizzi, L., Shaya, E. J., et al. 2009, The Astronomical Journal, 138, 323, doi: 10.1088/0004-6256/138/2/323
2. Anand, G. S., Lee, J. C., Van Dyk, S. D., et al. 2021b, MNRAS, 501, 3621, doi: 10.1093/mnras/staa3668
