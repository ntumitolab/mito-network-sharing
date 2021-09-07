# process XML output of analysed video data to get encounter network adjacency matrices
# 1.6 is threshold distance for an encounter in microns; -1 denotes no cap on trajectory lengths
Rscript trajectory-analysis.R Data/orig-1.xml 1.6 -1
Rscript trajectory-analysis.R Data/orig-2.xml 1.6 -1
Rscript trajectory-analysis.R Data/orig-3.xml 1.6 -1
Rscript trajectory-analysis.R Data/orig-4.xml 1.6 -1
Rscript trajectory-analysis.R Data/orig-5.xml 1.6 -1
Rscript trajectory-analysis.R Data/orig-6.xml 1.6 -1
# prune trajectories over 10 frames in length
Rscript trajectory-analysis.R Data/orig-1.xml 1.6 10

# housekeeping
mkdir Output
cp Data/*amlist.csv .

# compile simulation code
g++ bingo.c -I /usr/local/include/igraph -ligraph -lm -g -o ./bingo.ce

# run simulations:
# reference case with simulation and example trajectory output 
./bingo.ce orig-1.xml-amlist.csv --output-simulation-traces --output-trajectories 10 > tmpo1 &
# reference case with singletons removed
./bingo.ce orig-1.xml-amlist.csv --output-simulation-traces --remove-singletons --output-trajectories 10 > tmpo1a &
# reference case allowing simulated mitochondria to be inactive
./bingo.ce orig-1.xml-amlist.csv --inactivity > tmpo1b &
# reference case with pruned trajectories
./bingo.ce orig-1.xml-prune-10-amlist.csv > tmpo1c &

# other networks
./bingo.ce orig-2.xml-amlist.csv > tmpo2 &
./bingo.ce orig-3.xml-amlist.csv > tmpo3 &
./bingo.ce orig-4.xml-amlist.csv > tmpo4 &
./bingo.ce orig-5.xml-amlist.csv > tmpo5 &
./bingo.ce orig-6.xml-amlist.csv > tmpo6 &

