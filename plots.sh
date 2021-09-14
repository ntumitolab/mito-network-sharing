# visualise bio graphs
Rscript bingo-biograph-script.R Output/orig-1.xml-amlist.csv-0-0-graphs.txt
Rscript bingo-biograph-script.R Output/orig-2.xml-amlist.csv-0-0-graphs.txt
Rscript bingo-biograph-script.R Output/orig-3.xml-amlist.csv-0-0-graphs.txt
Rscript bingo-biograph-script.R Output/orig-4.xml-amlist.csv-0-0-graphs.txt
Rscript bingo-biograph-script.R Output/orig-5.xml-amlist.csv-0-0-graphs.txt
Rscript bingo-biograph-script.R Output/orig-6.xml-amlist.csv-0-0-graphs.txt

# visualise synthetic graphs for reference case
Rscript bingo-allgraphs-script.R Output/orig-1.xml-amlist.csv-0-0-graphs.txt

# visualise graphs for variants on reference case: singletons removed, trajectories pruned
Rscript bingo-allgraphs-script.R Output/orig-1.xml-amlist.csv-0-1-graphs.txt
Rscript bingo-allgraphs-script.R Output/orig-1.xml-prune-10-amlist.csv-0-0-graphs.txt

# final performance of different networks
Rscript bingo-results-script.R Output/orig-1.xml-amlist.csv-0-0-results-overall.txt 0
Rscript bingo-results-script.R Output/orig-2.xml-amlist.csv-0-0-results-overall.txt 0
Rscript bingo-results-script.R Output/orig-3.xml-amlist.csv-0-0-results-overall.txt 0
Rscript bingo-results-script.R Output/orig-4.xml-amlist.csv-0-0-results-overall.txt 0
Rscript bingo-results-script.R Output/orig-5.xml-amlist.csv-0-0-results-overall.txt 0
Rscript bingo-results-script.R Output/orig-6.xml-amlist.csv-0-0-results-overall.txt 0
# test cases: singletons removed, trajectories pruned
Rscript bingo-results-script.R Output/orig-1.xml-amlist.csv-0-1-results-overall.txt 0
Rscript bingo-results-script.R Output/orig-1.xml-prune-10-amlist.csv-0-0-results-overall.txt 0

# dynamics for different L in bio networks
Rscript bingo-biotrace-script.R Output/orig-1.xml-amlist.csv-0-0-results-frame.txt 0

# dynamics for different L across all networks (also with singletons removed)
Rscript bingo-alltraces-script.R Output/orig-1.xml-amlist.csv-0-0-results-frame.txt 0
Rscript bingo-alltraces-script.R Output/orig-1.xml-amlist.csv-0-1-results-frame.txt 0

# dynamics for m > 0
Rscript bingo-alltraces-script.R Output/orig-1.xml-amlist.csv-0-0-results-frame.txt 0.02
Rscript bingo-biotrace-script.R Output/orig-1.xml-amlist.csv-0-0-results-frame.txt 0.01
Rscript bingo-results-script.R Output/orig-1.xml-amlist.csv-0-0-results-overall.txt 0.01
Rscript bingo-biotrace-script.R Output/orig-1.xml-amlist.csv-0-0-results-frame.txt 0.02
Rscript bingo-results-script.R Output/orig-1.xml-amlist.csv-0-0-results-overall.txt 0.02

# predictions of bingo score
Rscript bingo-prediction-script.R Output/orig-1.xml-amlist.csv-0-0-graphs.txt Output/orig-1.xml-amlist.csv-0-0-results-frame.txt

# correlations of network statistics with bingo performance
Rscript bingo-correlations-script.R Output/orig-1.xml-amlist.csv-0-0-results-overall.txt
