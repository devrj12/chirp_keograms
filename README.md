# chirp_keograms

The algorithms in this repository can be used to make Range-Time-Intensity plots or keograms from the HF observations made by the ionospheric sounding mode of the HamSCI Personal Space Weather Station (PSWS) in Scranton, Pennsylvania currently implemented on an Ettus N200 Universal Software Radio Peripheral (USRP) using the open source GNU Chirpsounder data collection and analysis code chirpsounder2. 

1. filter_ionograms.py : It runs through the lfm files created during data collection using chirpsounder2 and saves the data in a dictionary. 

2. plot_ionograms.py : It reads the data saved by filter_ionograms.py to make the stacked subplots of the RTI images (keograms). 

For both of the python codes, there are folders named to read/save data from/to. These folder names need to be changed as according to the names of the folders to be used by other users of these python codes.

Example RTI Plot:![RTI-2021-02-18k](https://user-images.githubusercontent.com/15792043/131446616-0458d7ac-6893-4502-a50c-f5c1ad10feda.png)
