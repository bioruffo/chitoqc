# chitoqc
Chip Topology QC for Ion Proton chips  

This script creates an .mp4 movie of active/inactive Ion Torrent chip wells, to perform a visual analysis similar to the "Per tile sequence quality" of FastQC, but for aligned Ion Torrent data (e.g. exomes, gene panels).  

It very simply uses CIGAR data from alignments, so it doesn't have a fine flow-by-flow resolution; on the plus side, it can be ran even with an input BAM file that has no flow information.  


## Prerequisites:  
The script uses [FFmpeg](https://ffmpeg.org/).


## Usage:  
1) If needed, alter any of the settings in config.txt.   
2) Run:  
  `python ChiToQC.py -i <BAM file> -o <movie name>`  
  The mp4 movie will be saved in the output directory indicated in config.txt.  
  
## Example output

https://user-images.githubusercontent.com/24945128/160299016-5f817e81-a5de-480b-9da3-b020cd06ad48.mp4

The color scale (blue to red) represents the number of events for each position in the chip.  
The first row of images shows, from left to right: Insertions/deletions, Skips/Clips, Aligned wells, and Expired wells (that have reached the end of the sequence).  
The second row shows the same information as a ratio to the total amount of active wells in each position.  
The third row shows the difference between the current frame and the last one.  
The fourth row shows the % of wells as graph, through time. 
