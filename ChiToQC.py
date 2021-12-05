'''This script performs a topology quality control on an Ion Proton PI chip, by
visualizing sequence quality metrics mapped to chip locations.


'''

import sys, os, getopt
import pysam
import re
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from collections import OrderedDict
import time
import bisect


DEFAULTS = OrderedDict([("LOCATION", None),
                        ("START", None),
                        ("STOP", None),
                        ("BRACKMIN", 0),
                        ("BRACKMAX", None),
                        ("BRACKSTEP", 5),
                        ("QUALITY", 4),
                        ("IMGDIR", "ChiToQC_img")])


# Code to be returned by ReadTopology.retrieve_pos_cat() when sequence ended
NUC_OVER_LEN = 64

LOGFILE = "ChiToQC_log.txt"   

class ReadTopology(object):
    '''The ReadTopology class includes methods to store and retrieve topology
    data from each read in a BAM file.
    
    '''
    
    # Read locations on Ion PI BAMs are currently stored in the qname property
    # in the form: 'LABEL:X:Y' and can be retrieved by the pattern:
    pattern = re.compile('\w+:(\d+):(\d+)')
    
    
    def __init__(self, read):
        '''Extract data from 'read', a pysam.AlignedSegment object.
        Data collected is:
         - X, Y position in chip, as tuple
         - cigartuples, as numpy array
         - cigar "blocks" - i.e. nucleotides 0-2 are similar, then 3-10, then...
         - sequence length
        
        '''
        self.position = tuple(int(i) \
                              for i in self.pattern.search(read.qname).groups())
        self.cigartuples = np.array(read.cigartuples)
        self.blocks = np.cumsum(self.cigartuples[:,1])
        self.length = self.blocks[-1]
    
    def retrieve_pos_cat(self, nucleotide):
        '''Iterate over the cigartuples and return the cigar numeric code
        for the current position.
        Code 64 (NUC_OVER_LEN) reserved for sequence ended.
        
        '''
        if nucleotide >= self.length:
            return NUC_OVER_LEN
        else:
            return self.cigartuples[bisect.bisect(self.blocks, nucleotide)][0]

class TopologyData(object):
    '''The TopologyData class stores methods to analyze and plot the
    topology data of a run, based on locations saved in a BAM file.
    
    '''
    
    #Placeholder for no data, helps with max, np.amax and if
    nodata = [0]

  
    # Categories to be plotted; change the order here to switch plot order
    cigarcats = OrderedDict([("In, Del, Pad", [1, 2, 6]), # Insertions, deletions, padding: category 1
                             ("Skip, Clip", [3, 4, 5]), # Skipped, soft-clipped, hard-clipped: category 2
                             ("Aligned", [0, 7, 8]), # Alignment match/mismatch: category 0
                             ("Expired", [NUC_OVER_LEN])]) # sequence ended
    
    # Sequential keys
    plotcats = dict([(i, title) for i, title in enumerate(cigarcats.keys())])
    
    # Categorized cigar codes
    cigarcatkeys = dict([(cigar, category) for category, values in enumerate(cigarcats.values()) for cigar in values])
    
    # other plots
    variations = ["", " / load data", " (diff)", " / %alive, by cycle"]

    # Plot boundaries
    PI_extent = [0, 15456, 0, 10656] # The PI chip has 15456x10656 wells



    def __init__(self, bamfile, location=None, start=None, stop=None, mag=2):
        '''Read a BAM file and extract data from the reads within the desired region.
        'location' is the chromosome to be returned
        'start', 'stop' are nucleotide boundaries for pysam.fetch()
        'mag' is the magnification from the minimal plot size (see self.bins)
        
        The data extracted (and saved as numpy array) is:
         - a list of reads as ReadTopology objects, as 'reads';
         - a numpy array of each read's positions, as 'all_positions';
         - the total number of reads, as 'numreads';
         - the maximum length of a read, as 'maxlen'.
        A "loading" heatmap will also be created at initialization time.
               
        '''
        self.bins = [29*mag, 20*mag] #number of x and y bins for the histogram2d plot; ratio equals the PI x:y ratio
        self.empty_heatmap = np.zeros(self.bins)
        self.filename = bamfile.split(r"/")[-1]
        linelog("Reading BAM file: "+bamfile)
        sys.stdout.flush()
        if os.path.isfile(bamfile):
            bamfile = pysam.AlignmentFile(bamfile, "rb")
        else:
            linelog("\nFile {0:s} not found.".format(bamfile))
            sys.exit(2)
        
        # Data retrieval starts with a fetch, pointed to a region if needed
        if not location: #Retrieve the whole data
            selected_region = bamfile.fetch()
        elif location and not (start or stop): # Retrieve one chromosome
            selected_region = bamfile.fetch(location)
        else:
            selected_region = bamfile.fetch(location, start, stop)
        
        # Create the reads list
        self.reads = [ReadTopology(read) for read in selected_region]
        linelog("Storing read locations...")
        self.all_positions = np.array([read.position for read in self.reads])
        self.numreads=len(self.reads)
        self.maxlen = max([read.length for read in self.reads])
        linelog(" Done ("+str(len(self.reads)) + " reads).")
        # Create loading heatmap
        self.make_loadmap()
        

    def make_loadmap(self):
        '''Creates the "loading" heatmap, that includes all available reads.
        Return the heatmap as 'load_heat'.
        
        '''
        linelog("Creating the load map...")
        self.lastbins = self.bins
        x = self.all_positions[:,1]
        y = self.all_positions[:,0]
        self.load_heat, xedges, yedges = np.histogram2d(x, y, bins=self.bins)
        linelog(" Done.")


    def categorize_reads(self, start, step=1):
        '''Separate reads by category, in order to be plotted.
        If 'step' is provided, merge several points into one.
        Return categories as sets of (x,y) coordinates, in the form
        dict(category: np.array([(point_tuples), ...]))        
        
        '''
        # first create an empty dictionary populated with the right categories
        categorized_data = dict([(category, set()) for category in self.plotcats])
        # For each step point, pass through the reads and categorize them
        this_point = np.zeros((self.numreads, step+2), dtype=np.uint16)
        this_point[:,0:2] = self.all_positions
        for point in range(step):
            # Create an array containing the category to which each read belongs
            this_point[:,point+2] = [self.cigarcatkeys[read.retrieve_pos_cat(start+step)] for read in self.reads]
        for category in categorized_data.keys():
            categorized_data[category] = this_point[np.array([x for x in map(any, zip(*[this_point[:,i] == category for i in range(2,2+step)]))], dtype=np.bool)]
            categorized_data[category] = categorized_data[category][:,:2]
        self.test = categorized_data
        return categorized_data


    def default_plot(self, ax, heatmap, cmap, vmax, vmin):
        # "jet" or "rainbow"
        ax.pcolor(heatmap.transpose(), cmap=cmap, vmax=vmax, vmin=vmin)
    
    def plot_data(self, heatmaps, scattrange, scatters, maxes, mins, title="", savepic=None):
        # TODO this has more arguments than a divorce hearing
        #from http://stackoverflow.com/questions/2369492/generate-a-heatmap-in-matplotlib-using-a-scatter-data-set
        # altered from: http://stackoverflow.com/a/11367847

        fig, ax = plt.subplots(4, 4, figsize=(12,8))
        fig.suptitle("{0:s} cigarstring position(s) {1:s}".format(self.filename, title), fontsize=14)
        hmp_pencils = ["jet", "jet", "coolwarm"]
        # Plotting load data at 0,0
        #self.default_plot(ax[0, 0], self.load_heat, vmax=None)
        # Plotting other categories
        for i in self.plotcats:
            for j in range(3):
                print(i, j)
                self.default_plot(ax[j, i], heatmaps[j][i], cmap=hmp_pencils[j], vmax=maxes[j][i], vmin=mins[j][i])
            ax[3, i].plot(scattrange, scatters[i])
        # remove axes and set titles
        for i in range(4):
            for j in range(4):
                if i < 3:
                    ax[i,j].axis('off')
                    #ax[i,j].set_aspect(1)
                else:
                    ax[i,j].tick_params(axis='x', labelsize=8)
                    ax[i,j].tick_params(axis='y', labelsize=6)
                    #ax.set_aspect('auto')
                    
                ax[i,j].set_title(self.plotcats[j] + self.variations[i], fontsize=10)
        if savepic:
            plt.savefig(savepic)
        plt.close(fig)


    def make_heatmaps(self, point_data):
        curr_point_hmps = []
        for i in self.plotcats.keys(): # this will break if you mess wiht plotcats
            if point_data[i].any():
                # Plotting total values at 2,1..2,4
                x = point_data[i][:,1]
                y = point_data[i][:,0]
                heatmap, xedges, yedges = np.histogram2d(x, y, bins=self.bins)
                curr_point_hmps.append(heatmap)
            else:
                # Need to leave a placeholder
                curr_point_hmps.append(self.empty_heatmap)
        return curr_point_hmps
    
    def make_ratios_heatmaps(self, point_heatmaps):
        base_heatmap=self.load_heat
        curr_point_ratio_hmps = []
        for heatmap in point_heatmaps:
            if heatmap is not self.nodata:
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.true_divide(heatmap, base_heatmap)
                    ratio[ratio == np.inf] = 0
                    ratio = np.nan_to_num(ratio)
                curr_point_ratio_hmps.append(ratio)
            else:
                curr_point_ratio_hmps.append(self.empty_heatmap)
        return curr_point_ratio_hmps
        
    def make_diff_heatmaps(self, last, curr):
        # TODO ok this is so repeated
        curr_point_diff_hmps = []
        for i in range(len(curr)):
            if last[i] is not self.nodata:
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.subtract(curr[i], last[i])
                    ratio = np.nan_to_num(ratio)
                curr_point_diff_hmps.append(ratio)
            else:
                curr_point_diff_hmps.append(self.nodata)
        return curr_point_diff_hmps


    def make_movie(self, directory, moviename, brackstep=1, brackmin=0, brackmax=None):
        brackmin = 0
        img_dir = directory+"/"
        call(["rm", "-r", img_dir])
        call(["mkdir", img_dir])
        if not brackmax or brackmax < brackmin:
            brackmax = self.maxlen
        else:
            brackmax = min([brackmax, self.maxlen])
       # Create heatmaps first, to check maxes
        linelog("Creating heatmaps:")
        titles=[]
        # TODO it's a prototype, bear with me I have few time
        heatmaps=[[],[],[]] # hmp, ratio, diff
        scatters = [[],[],[],[]] # one per category
        # TODO horrifying [[0 for i in self.plotcats], [0 for i in self.plotcats], [0 for i in self.plotcats]]
        heatmaps_maxes = [[0 for i in self.plotcats], [0 for i in self.plotcats], [0 for i in self.plotcats]]
        heatmaps_mins = [[0 for i in self.plotcats], [0 for i in self.plotcats], [0 for i in self.plotcats]]
        # Recreate the Load data, in case the bins have been altered
        if self.lastbins != self.bins:
            self.make_loadmap()
        past_ratios = None
        for i, min_limit in enumerate(range(brackmin, brackmax, brackstep)):
            times=[0, 0, 0, 0, 0, 0, 0] #TODO even this is ugly
            starttime = time.time()
            self.filtered = self.categorize_reads(min_limit, brackstep)
            times[0]=int(time.time() - starttime)
            # Creating current point heatmaps
            inittime = time.time()
            curr_heatmaps=self.make_heatmaps(self.filtered)
            heatmaps[0].append(curr_heatmaps)
            times[1] = int(time.time()-inittime)
            # Creating ratios heatmaps
            inittime = time.time()
            curr_ratios = self.make_ratios_heatmaps(curr_heatmaps)
            heatmaps[1].append(curr_ratios)
            times[2] = int(time.time()-inittime)
            # Creating differential heatmaps
            inittime = time.time()
            if i == 0:
                curr_diff = [self.empty_heatmap for i in range(4)]
            else:
                curr_diff = self.make_diff_heatmaps(past_ratios, curr_ratios)
            heatmaps[2].append(curr_diff)
            past_ratios = curr_ratios
            times[3] = int(time.time()-inittime)
            # Making scatterplots
            inittime = time.time()
            sum_alive = sum([np.sum(this_array) for this_array in curr_heatmaps[:3]])
            for j in self.plotcats:
                if j < 3:
                    scatters[j].append(np.sum(curr_heatmaps[j])/sum_alive)
                else:
                    scatters[j].append(np.sum(curr_heatmaps[j]))
            times[4] = int(time.time()-inittime)
            # Adding titles
            titles.append(str(min_limit) + " to " + str(min_limit+brackstep-1) + " sequence bases")
            # Finding mins and maxes
            inittime = time.time()
            for j in self.plotcats:
                for k in range(3):
                    sort_9x=np.sort(heatmaps[k][i][j].flat)
                    sort_9x=sort_9x[int(len(sort_9x)*0.02):int(len(sort_9x)*0.98)]
                    heatmaps_maxes[k][j] = max([heatmaps_maxes[k][j], np.amax(sort_9x)])
                    heatmaps_mins[k][j] = min([heatmaps_mins[k][j], np.amin(sort_9x)])
            times[5] = int(time.time()-inittime)
            times[6] = int(time.time()-starttime)
            linelog("Heatmap: {0:s}, filter {1:d}s, maps {2:d}s, ratios {3:d}s, diff {4:d}s, scat {4:d}s, max {5:d}s, tot {6:d}s".format(titles[-1], *times))
        linelog("Done.")
        linelog("Plotting and saving figures...")
        scatters=np.array(scatters)
        scattrange=[i+brackstep-1 for i in range(brackmin, brackmax, brackstep)]
        for i, min_limit in enumerate(range(brackmin, brackmax, brackstep)):
            savepic = img_dir+str(i).rjust(3,"0")+".png"
            linelog("Plotting fig " + savepic + ", " + titles[i])
            # TODO this is awful: [heatmaps[0][i], heatmaps[1][i], heatmaps[2][i]]
            self.plot_data([heatmaps[0][i], heatmaps[1][i], heatmaps[2][i]], scattrange[:i+1], scatters[:,:i+1], heatmaps_maxes, heatmaps_mins, title=titles[i], savepic=savepic)
        linelog("Done.")
        linelog("Creating movie...")
        call(["ffmpeg", "-framerate", "4", "-i", directory+"/%03d.png", moviename, "-y"])
        linelog("Done!")


def linelog(message, mode="a"):
    print(message)
    message = message.replace("\n", "\n                    | ")
    logfile = open(LOGFILE, mode)
    logfile.write("{1:d}_{2:0>2d}_{3:0>2d}_{4:0>2d}:{5:0>2d}:{6:0>2d} | {0:s}\n".format(message, *time.localtime()[0:6]))
    logfile.close()


def load_configs(filename):
    configs = DEFAULTS
    numconfigs = ["START", "STOP", "QUALITY", "BRACKSTEP", "BRACKMIN", "BRACKMAX"]
    linelog("Loading configurations...")
    try:
        with open(filename, "r") as configs_file:
            lines = [line.split('#')[0] for line in configs_file.readlines()]
        configs_file.close()
        for line in lines:
            line_data = [item.strip(' \t\n\r') for item in line.split("=") if "=" in line]
            if len(line_data) == 2:
                if line_data[0] in numconfigs:
                    configs[line_data[0]] = int(line_data[1])
                else:
                    configs[line_data[0]] = line_data[1]
    except FileNotFoundError:
        message = "### WARNING ###\nConfiguration file {0:s} not found; using built-in defaults.".format(filename)
        linelog(message)
    message = ""
    for option, value in configs.items():
        if value is not None:
            message = message + option + ": "+str(value)+"\n"
        else:
            message = message + option + ": None\n"
    linelog(message.strip("\n"))
    return configs

def main(argv, configs):
    # BAM name must always be provided
    bamname = None
    # Default output file name:
    moviename = 'ChiToQC_output.mp4'
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('ChiToQC.py -i <BAM file> -o <movie name>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            bamname = arg
        elif opt in ("-o", "--ofile"):
            moviename = arg
    if not bamname:
        print('Please provide a name for the BAM file.')
        sys.exit(2)

    my_data = TopologyData(bamfile=bamname,
                           location=configs["LOCATION"],
                           start=configs["START"],
                           stop=configs["STOP"],
                           mag=configs["QUALITY"])
                           
    my_data.make_movie(directory=configs["IMGDIR"],
                       moviename=moviename,
                       brackstep=configs["BRACKSTEP"],
                       brackmin=configs["BRACKMIN"],
                       brackmax=configs["BRACKMAX"])
    linelog("ChiToQC ended.")
    

if __name__ == "__main__":
    # initializing logfile
    linelog("ChiToQC started", mode="w")
    # Load configuration files
    configs=load_configs("config.txt")
    # Run the main script
    main(sys.argv[1:], configs)
