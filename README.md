# CATS-H0

Based on "sn_analysis_new.csv", there is a two-step detection to give TRGB of each field.

We have different spatial clipping extents, color band widths, smoothing scales, minimum relative thresholds of peaks and internal extinction corrections. For each of the combinations, there is a detection file providing the TRGB magnitudes, and an information file providing essential statistics of all the fields.

We also provide the initial csv files in "phot_csv" folder which pass through EDD photometry selection, clipped csv files in "clipped_csv" folder, files of stars removed after spatial clipping in "cut_csv" folder, and spatial clipping maps in "clipping_map" folder.

To reproduce our results, just run "Comprehensive_part1.py" and "Comprehensive_part2.py" in turn.
