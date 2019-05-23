# ResCap_analysis

ResCap_analysis_1.0

Hello world! 
This script allows you to analyze the results from the ResCap capture platform. It performs an on-target selection, then a
subsampling to normalize the samples and erase outliers, and then a selection of the resistance genes with KMA.

Let's hope it works! ¯\_(ツ)_/¯

Software needed to work:
  Python3
  Pandas
  Bowtie2

The input can be added by having all the files with the labels R1 and R2 in a folder, or a csv file with two columns named R1 and R2 with the names or path of the files. (all of them have to be contained in the same folder).

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        Number of clusters available to use
  -t TABLE, --table TABLE
                        CSV files (comma separated) with two columns, called
                        R1 and R2
