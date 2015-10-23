# LC-ICPMS-ESIMS-feature-detection
R scripts for aligning LC-ICPMS and LC-ESIMS data and identifying chromatographically coherent isotope features

Requires XCMS package for LC-ESIMS data analysis.

LC-ICPMS data should be imported as a csv file with separate columns for the time and intensity of each isotope X. The headers must be in the format "Time X" and "X" (e.g. "Time 56Fe",  "56Fe").

Alignment feature requires a cyanocobalamin peak in both data sets.
