# $${\color{lightgreen}Methyl\\_arrays.r}$$

The arguments should be placed in the following order:

                                       1 - a directory where microarray data are stored,
                                       2 - a csv file with microarray data description,
                                       3 - a directory where the analysis results should be saved,
                                       4 - a Boolean value (TRUE/FALSE) determining if the data binarization step should be performed,
                                       5 - a number of CPU threads to be used,
                                       6 - a variable in the aforementioned csv file containing sample names,
                                       7 - a variable in the aforementioned csv file containing sample sources,
                                       8 - a variable in the aforementioned csv file containing sample genders,
                                       9 - a variable in the aforementioned csv file containing sample descriptions,
                                      10 - independent factor variables in the aforementioned csv file to be used, separated with '+', starting with a main categorical variable,
                                      11 - a variable containing the Sentrix IDs,
                                      12 - a variable containing the Sentrix Positions,
                                      13 - a path to the optional txt file containing the list of CpG sites (one per a line) for which the methylation status visualization is to be performed (NA if none).
