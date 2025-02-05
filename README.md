# Predictive ability of hop (Humulus lupulus L.) grown in single hills on plot environments
</p>
Crop Science, 2025 <br>
DOI: 
</p>

#### About the study
In this study, we evaluated seven hop genotypes in five different spacing and trellis configurations common to hop breeding programs. The study took place in 2022 and 2023 in Prosser, WA, U.S. A number of agronomic and cone-related traits were
evaluated. The objective of the study was to evaluate the predictve ability of early stage selection environments on later stage, plot environments to determine if selection earlier in the pipeline could be sufficiently
predictive and possibly justify the elimination of the time, space, and labor intensive "single-hill" stages. 

#### What's inlcluded
The R script, "Data Analyses.R," utilzies "historical_blups.csv", and "Spearman Rank Correlations Table.csv", in addition to the dataset "spacing_trial_dat.csv" available on the dryad repository associated with the paper: DOI: 10.5061/dryad.dbrv15fcg. 

#### About the script
The script is long, but includes all the analyses in the paper. The analyses are separated by #### headers #### to allow for easy navigation in R Studio. 

First, it evaluates correlations between cone traits (to determine which to report), prepares the data for running the linear mixed effects models, and calcuates the means of the 3.5' spacings. 
Second, it evaluates the genotype's performance relative to their historical yield data to confirm the validity of our "gold standard" as described. 
Third, it evaluates the effect of position within the 3.5' spacing plot. 

It then iterates through each of the traits, running the linear models and ANVOAs to evaluate the differences between spacings, spearman rank correlations between genotypes/cultivars, and creates the bar plots. 

Finally, it compiles all the means and ANOVAs for the supplemental, and explores a few more questions: what happens when you extrapolate yield from single hills; creates a figure to represent the correlations in an easier to digest format, 
and summarizes the frequency of the significance of the relationships. 
