/*------------------------------------------------------------------------------
 * timepix3cam version 3
 * Copyright (C) 2016-2022
 * Liisa Hirvonen, CMCA, The University of Western Australia
 * lmhirvonen@gmail.com
 *----------------------------------------------------------------------------*/

/* This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*------------------------------------------------------------------------------
 * USAGE:
 * <last file> <input file(s)>
 *
 * <last file> is the number of last file to be processed (in each subdirectory),
 *	use zero to process all files
 * <input file(s)> is an absolute (full) path to the input file. For multiple files,
 *  enter the full path for each file.
 * The file extension should be .csv (ASCII files)
 *----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Compiling in Ubuntu 12.04 / 16.04 / 18.04 with Code::Blocks
 *  - Requires libtiff-dev (libtiff5-dev) and g++ packages
 *  - Add libm.so and libtiff.so in build bath
 *    (in Build options -> Linker settings -> Link libraries)
 *----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 * WHAT IT DOES:
 * Reads ASCII files with columns of photon event data.
 * This program was specifically made for processing data from Timepix3Cam.
 * Organises the data into an xyt cube.
 * The output is an .ics file that can be loaded in Tri2 for fluorescence
 * lifetime analysis.
 * Also outputs an ASCII file containing a histogram of the photon arrival times
 * and an intensity image.
 * Besides sum (all pixels), also outputs centroided data, where only the centroid
 * pixel of each cluster (photon event) is counted.
 *----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 * Changelog
 *
 *  Version 0.0.0: (April 2, 2019) * - Copy of timepixcam-4.0.1.c with old data format
 *  Version 1.0.0: (April 2, 2019) * - Changed input data format to timepix3cam
 *  Version 1.0.1: (April 11, 2019) * - New time calibration - the time is not in ps
 *  Version 1.0.2: (April 26, 2019) * - Corrected image flip & rotation * - Changes to outparam file
 *  Version 1.0.3: (April 26, 2019) * - Sub-pixel centroiding
 *  Version 1.0.4: (June 13, 2019) * - Fixed ICS header, now opens with SLIM Curve
 *  Version 1.0.4e: (Feb 29, 2020) * - Change output to print total photon numbers (for Rajannya's paper)
 *  Version 2.0: (Jul 23, 2020) * - New detector in Cork - copy of 1.0.4e * - Changed p_flimhistp from short to int (fix overflow)
 *  Version 2.1: (Sept 23, 2020) * - Add histogram for photons / frame
 *  Version 2.2: (Oct 2, 2020) * - Remove bright pixels from image
 *  Version 2.3: (Oct 2, 2020) * - Option to process more than 1 input file * - Reorganised input arguments
 *  Version 2.4: (Oct 6, 2020) * - Add histogram for photons / frame
 *  Version 2.4a: (Oct 7, 2020) * - Testing to find source of double peak/bump in decay
 *  Version 2.5: (Oct 9, 2020) * - Add count rate to outparam file
 *  Version 3.0: (Nov 2, 2020) * - Cleaned up version for compile on Windows
 *  Version 3.1: (Oct 13, 2022) * - Cleaned up version for general distribution
 *
 */

#define VERSION  "3.1"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <unistd.h>
#include <stdbool.h>
#include <signal.h>
#include <fcntl.h>
#include <time.h>

#ifndef O_BINARY
# define O_BINARY 0
#endif

// Input parameters
#define TIMEBINS_IN 200000000 // Time points in input file
#define TIMEBIN_IN_US 0.000006103516 // Time bin width =  t/4096*25 ns
// ***Note: ToT is already in ns
#define DIMX 256    // Input file dimension
#define DIMY 256    // Input file dimension

// For processing part of the image
// NOTE: if changing image size, remember to check sum image pixel removal
#define X_START	0
#define Y_START 0
#define XY_WIDTH_IN 256
#define T_START  0
#define BINNING_TIME 100000  // Binning for time axis
#define TIMEBINS_OUT_NOBIN 32800000 // Without binning

// Error codes
#define ERR_USAGE 1
#define ERR_IMG 2
#define ERR_BADDIR 2
#define ERR_ALLOC 3
#define ERR_INFILE 4
#define ERR_OUTFILE 5
#define ERR_BADCOUNT 6
#define ERR_NOEVENT 10

int abort_received = 0;

// Output parameters
int TIMEBINS_OUT = TIMEBINS_OUT_NOBIN/BINNING_TIME;
double TIMEBIN_OUT_US = TIMEBIN_IN_US*BINNING_TIME;
int XY_WIDTH_OUT = XY_WIDTH_IN + 2;
int X_END = X_START + XY_WIDTH_IN;
int Y_END = Y_START + XY_WIDTH_IN;
int T_END = T_START + TIMEBINS_OUT_NOBIN;
double timebin_out_s = TIMEBIN_IN_US*BINNING_TIME/1000000; // convert from us to s

void sigproc(int n);

void printUsage();

void make_ics_header(FILE *fp, char *name, int channels, float length, int dimY, int dimX);

void FindSurroundingPixelsAmplitude(unsigned short x, unsigned short y,
    unsigned short dimx, unsigned short *p_area, unsigned int *p_xc,
    unsigned int *p_yc, unsigned short *p_xp, unsigned short *p_yp,
    unsigned int *p_intensity, unsigned int *p_maxint, bool p_donetab_a[],
    const unsigned int p_tab_a[], unsigned int *maxarea);

void FindSurroundingPixelsTime(unsigned short x, unsigned short y,
    unsigned short dimx, unsigned short *p_area, unsigned int *p_xc,
    unsigned int *p_yc, unsigned short *p_xp, unsigned short *p_yp,
    unsigned int *p_intensity, unsigned int *p_maxint, bool p_donetab_t[],
    const unsigned int p_tab_t[], unsigned int *maxarea);

void init(char *argv[], int argd, int zoom,
    FILE **p_eventlist, FILE **p_outparam, // Text files
    FILE **p_imageavraw, FILE **p_imagecraw, FILE **p_imagetraw, // Image files
    FILE **p_flimimageavraw, FILE **p_flimimagepraw, // FLIM image files
    FILE **p_flimhist, FILE **p_counthist, // Histogram files
    unsigned int **p_tab_a, unsigned int **p_tab_t,    // Arrays for current frame
    bool **p_donetab_a, bool **p_donetab_t,   // Arrays for current frame
    unsigned long int **p_imageav, unsigned int **p_imagec, unsigned int **p_imaget, // Arrays for intensity images
    unsigned int **p_flimimageav, unsigned int **p_flimimagep, // Arrays for FLIM images
    unsigned int **p_flimhistav, unsigned int **p_flimhistp); // Arrays for FLIM images

int main(int argc, char *argv[]) {

	int x, y, time, t, argd, zoom, numFiles;
	unsigned int xca, yca, xct, yct;    // Weighted average centroid coordinates
	unsigned short xpa, ypa, xpt, ypt;  // Earliest timecode pixel coodrinates
	unsigned short no_blank=0, area=0, minarea=0;
	unsigned int intensity_a, intensity_t, maxint_a, maxint_t, maxarea=0;
	unsigned long int event=0, event_frame=0, blank_frame=0,
               discarded_under=0, discarded_under_frame=0,
               discarded_over=0, discarded_over_frame=0,
               discarded_outOfBounds=0;
	long long int counter=0, frame=0, emptyframe=0;
	unsigned long long int imagec_max=0, imaget_max=0, imageav_max=0, index_xy=0;
    long long int trigId=0, tot=0;
    unsigned long long int trigTime=0, toa=0;
    bool endOfFile = false;

	unsigned int *tab_a, *tab_t;  // Current frame (amplitude and time)
	bool *donetab_a, *donetab_t;      // Keeps track of processed pixels in current frame
	unsigned long int *imageav; // Array for sum intensity image. Size: (XY_WIDTH_OUT)^2
	unsigned int *imaget; // Array for timecode centroided intensity image. Size: (XY_WIDTH_OUT)^2
	unsigned int *imagec; // Array for amplitude centroided intensity image. Size: (XY_WIDTH_IN*zoom+2)^2
	unsigned int *flimimageav, *flimimagep; // Data cubes for flim images. Size: (XY_WIDTH_OUT)^2*TIMEBINS_OUT
	unsigned int *flimhistav;   // Arrays for flim histograms. Size: TIMEBINS_OUT
	unsigned int *flimhistp;   // Arrays for flim histograms. Size: TIMEBINS_OUT

	FILE *in, *eventlist, *outparam;
	FILE *imageavraw, *imagecraw, *imagetraw, *flimimageavraw, *flimimagepraw, *flimhist, *counthist;

    signal(SIGINT, sigproc);

	// Check number of input arguments
	//if (argc < 6) { printUsage(argc, argv); }
	//else {argd = 5;} // 1st argument that contains the file path
	if (argc < 3) { printUsage(argc, argv); }
	else {argd = 2;} // 1st argument that contains the file path

	// Read input arguments
	// Changed many of these for hard-coded ones for Rajannya compile test
	minarea = 3; //atoi(argv[1]);
	maxarea = 15; //atoi(argv[2]);
	zoom = 1; //atoi(argv[3]);
	long long int lastframe = atoi(argv[1]);//atoi(argv[4]);
	numFiles = argc - 2;//5;
	printf("\nNumber of input files = %d\n",numFiles);

    // Check that input variables area are sensible
	if (minarea < 1 || minarea > maxarea) {
		printf("\nERROR: Minimum area must be bigger than 0 and smaller than maxarea.\n");
		printf("Exiting...\n\n");
		exit(ERR_USAGE);
	}
	if (zoom < 1) {
		printf("\nERROR: Zoom must be bigger or equal to 1.\n");
		printf("Exiting...\n\n");
		exit(ERR_USAGE);
	}
	if (lastframe<1) lastframe=100000000000000;

    // Initialise: allocate memory for arrays, open files
    init(argv, argd, zoom,
        &eventlist, &outparam, // Text files
        &imageavraw, &imagecraw, &imagetraw, // Intensity image files
        &flimimageavraw, &flimimagepraw, // FLIM image files
        &flimhist, &counthist, // Histogram files
        &tab_a, &tab_t, &donetab_a, &donetab_t, // Current frame arrays
        &imageav, &imagec, &imaget, // Intensity image arrays
        &flimimageav, &flimimagep, // FLIM image arrays
        &flimhistav, &flimhistp); // FLIM histogram arrays

    int currentTrigId = 0;
    int avgEventCounter = 0;
    long long int avgTimeCounter = 0; // For average arrival time time calculation
    long long unsigned int fileStartTime = 0; // Calculation of duration of experiment in each file
    double totalTimeCounter = 0; // Duration of experiment across multiple files
    int eventsAvg = 0;
    int curFileFrame = 0;

    int startFrame = 0;     // For processing only part of the data
    int skippedFrame = 0;   // For processing only part of the data
    lastframe += startFrame;

    // Loop over input files
    for (int curFile=0; curFile<numFiles; curFile++) {
         // Open input file for reading
        in = fopen(argv[argd+curFile],"r");
        if (in==NULL) {
            printf("error opening %s\n\n", argv[argd+curFile]);
            exit(ERR_BADDIR);
        }
        printf("\nOpened %s for reading", argv[argd+curFile]);

        // Read the header
        char header_row[1024];
        fgets(header_row, sizeof header_row, in);

        // Read the first line
        if (fscanf(in, "%lld,%llu,%d,%d,%llu,%lld,%d,\n", &trigId, &trigTime, &x, &y, &toa, &tot, &time) != 7) {
            printf("\n\nUnable to read first line of data. Exiting...\n\n");
            exit(ERR_INFILE);
        }
        currentTrigId = trigId;
        avgEventCounter = 0;
        avgTimeCounter = 0;
        fileStartTime = trigTime;
        curFileFrame = 0;

        // Go through the file line by line
        while (curFileFrame<lastframe && !endOfFile) {

            frame++;
            curFileFrame++;

            // Check if user has pressed Ctrl+C, if so quit the loop
            if (abort_received == 1 ) {
                if ( frame != 0 || abort_received == 2) {
                    printf("Abort received - exiting...\n");
                    abort_received = 2;
                }
                else {
                    printf("No file processed\n");
                    exit(1);
                }
                break;
            }

            // Reset tabs and donetab
            bzero(tab_a, (XY_WIDTH_OUT)*(XY_WIDTH_OUT)*sizeof(unsigned int));
            bzero(tab_t, (XY_WIDTH_OUT)*(XY_WIDTH_OUT)*sizeof(unsigned int));
            bzero(donetab_a, (XY_WIDTH_IN+3)*(XY_WIDTH_IN+3)*sizeof(bool));
            bzero(donetab_t, (XY_WIDTH_IN+3)*(XY_WIDTH_IN+3)*sizeof(bool));

            // Loop for frame
            while (trigId==currentTrigId) {
                // Process the line of data and put it in correct pixel in time stacks
                if (frame > startFrame && x>X_START && x<X_END && y>Y_START && y<Y_END && time>T_START && time<T_END) {

                    // Add pixel to sum result arrays
                    flimimageav[(x-X_START) + XY_WIDTH_IN*(y-Y_START) + XY_WIDTH_IN*XY_WIDTH_IN*((time-T_START)/BINNING_TIME)] += 1;
                    flimhistav[(time-T_START)/BINNING_TIME] += 1;

                    // Increase counters
                    counter++;

                    // Add to current frame events
                    index_xy = (x-X_START) + XY_WIDTH_IN*(y-Y_START);
                    no_blank = 1;

                    // Add the pixel in the image pixel arrays
                    tab_a[index_xy] = tot; // For amplitude centroiding
                    tab_t[index_xy] = time; // For time centroiding and FLIM image

                    // Add to sum image, 
                    // Uncomment if/else to discard hot pixels
                    // coordinates from "_sum.raw" image: x=xIn, y=257-yIn
                    //if ((x==131) & (y==254)) {}
                    //else if ((x==254) & (y==254)) {}
                    //else if ((x==20) & (y==149)) {}
                    //else if ((x==7) & (y==1)) {}
                    //else if ((x==6) & (y==138)) {}
                    //else {
                        imageav[index_xy]+=tot; // Amplitude sum image
                        if (imageav[index_xy]>imageav_max && index_xy!=65278) {
                            imageav_max=imageav[index_xy]; // Update max
                        }
                    //}


                    if (fscanf(in, "%lld,%llu,%d,%d,%llu,%lld,%d,\n", &trigId, &trigTime, &x, &y, &toa, &tot, &time) != 7) {
                        printf("\nEnd of file.\n");
                        endOfFile = true;
                        trigId++;
                        //break;
                    }
                }
                else if (frame < startFrame) { // If not at start frame, do nothing
                   if (fscanf(in, "%lld,%llu,%d,%d,%llu,%lld,%d,\n", &trigId, &trigTime, &x, &y, &toa, &tot, &time) != 7) {
                        printf("\nEnd of file.\n");
                        endOfFile = true;
                        trigId++;
                        //break;
                    }
                }
                else { // End of file
                    discarded_outOfBounds++;
                    if (fscanf(in, "%lld,%llu,%d,%d,%llu,%lld,%d,\n", &trigId, &trigTime, &x, &y, &toa, &tot, &time) != 7) {
                        printf("\nEnd of file.\n");
                        endOfFile = true;
                        trigId++;
                        //break;
                    }
                }
            }	// end of while loop

            if (!endOfFile) {
                currentTrigId=trigId;

                // If there are no events in this frame, skip forward to next frame
                if ( no_blank == 0 || frame < startFrame) {
                    if (frame < startFrame) skippedFrame++; else blank_frame++;
                    event_frame=0;
                    discarded_under_frame=0;
                    discarded_over_frame=0;
                    fflush(eventlist);
                    continue;
                }

                // Scan the frame for events
                for (x=1; x<XY_WIDTH_IN-2; x+=1) {
                    for (y=1; y<XY_WIDTH_IN-2; y+=1) {

                        area=0;
                        intensity_a=0;
                        maxint_a=0;
                        intensity_t=0;
                        maxint_t=1000000000;
                        xca=0; yca=0;
                        xpa=0; ypa=0;
                        xct=0; yct=0;
                        xpt=0; ypt=0;

                        // If the pixel is over the threshold and hasn't already been checked,
                        // finds the event area, intensity and centre pixel.
                        if ( tab_a[x+y*XY_WIDTH_IN] !=0 && !donetab_a[x+y*XY_WIDTH_IN] ) {
                            FindSurroundingPixelsAmplitude(x, y, XY_WIDTH_IN, &area, &xca, &yca, &xpa, &ypa,
                                                  &intensity_a, &maxint_a, donetab_a, tab_a, &maxarea) ;
                        }

                        if ( tab_t[x+y*XY_WIDTH_IN] !=0 && !donetab_t[x+y*XY_WIDTH_IN] ) {
                            FindSurroundingPixelsTime(x, y, XY_WIDTH_IN, &area, &xct, &yct, &xpt, &ypt,
                                                  &intensity_t, &maxint_t, donetab_t, tab_t, &maxarea) ;
                        }

                        // Only take this event into account if its size is correct
                        if (area >= minarea && area <= maxarea) {

                            // Increase counters
                            event++;
                            event_frame++;

                            // Centre of mass calculation
                            xca = (int) zoom*xca/intensity_a;
                            yca = (int) zoom*yca/intensity_a;

                            // Write event to eventlist file
                            fprintf(eventlist,"%lld\t%d\t%d\t%d\t%d\t%hd\t%hd\n", frame, xpa, ypa, xca, yca, area, maxint_a);

                            // Add event to amplitude centroided image
                            imagec[xca+yca*XY_WIDTH_IN*zoom]++;
                            if (imagec[xca+yca*XY_WIDTH_IN*zoom]>imagec_max) {
                                imagec_max=imagec[xca+yca*XY_WIDTH_IN*zoom];
                            }

                            // Add event to timecode centroided image
                            imaget[xpt+ypt*XY_WIDTH_IN]++;
                            if (imaget[xpt+ypt*XY_WIDTH_IN]>imaget_max) {
                                imaget_max=imaget[xpt+ypt*XY_WIDTH_IN];
                            }

                            // Add to FLIM image
                            // NOTE: Now adds amplitude centroided pixel
                            flimimagep[xpa+XY_WIDTH_IN*ypa+XY_WIDTH_IN*XY_WIDTH_IN*((maxint_t-T_START)/BINNING_TIME)] += 1;
                            // Add to intensity histogram list
                            flimhistp[(maxint_t-T_START)/BINNING_TIME] += 1;

                            // Add to average arrival time counter
                            avgTimeCounter += (maxint_t-T_START);

                        }
                        // Discard if the event is too big or too small
                        else if ( area > maxarea ) {
                            discarded_over++; discarded_over_frame++;
                        }
                        else if ( area < minarea && area !=0 ) {
                            discarded_under++; discarded_under_frame++;
                        }
                    }
                } // end of scan for events

                if (frame%100!=0) {
                    avgEventCounter += event_frame;
                }
                else {
                    eventsAvg = avgEventCounter/100;
                    fprintf(counthist, "%lld\t%d\t%.3f\n", frame, eventsAvg, (avgTimeCounter/avgEventCounter)*TIMEBIN_IN_US);
                    avgEventCounter = 0;
                    avgTimeCounter = 0;
                }

                event_frame=0;
                discarded_under_frame=0;
                discarded_over_frame=0;

                fflush(eventlist);
            }
        }	// end of file loop
        endOfFile = false;
        totalTimeCounter += (trigTime-fileStartTime);
    } // end of loop over files

    // Calculate total experiment duration in s
	double totalTimeS = totalTimeCounter * TIMEBIN_IN_US / 1000000;

    fprintf(stderr,"\nWriting image files... \n");

    // Display the progress (percentage of image files written)
    int num_imgs = 7;
    fprintf(stderr,"\b\b\b%d%%", 0*100/num_imgs); // Update progress %

	/* Write the output data matrices into the output files */

    // FLIM histogram, i.e. binned time decay for the whole image
	double timeinterval = TIMEBIN_OUT_US; // Keep the time in microseconds
	for (t=0; t<TIMEBINS_OUT; t++) {
		fprintf(flimhist, "%.9f\t%u\t%u\n", t*timeinterval, flimhistav[t], flimhistp[t]);
	}
	fprintf(flimhist, "*END\n");
    fprintf(stderr,"\b\b\b%d%%", 1*100/num_imgs); // Update progress %

	// Sum intensity image
	for (y=XY_WIDTH_OUT-1; y>=0; y--) {
		for (x=0; x<XY_WIDTH_OUT; x++) {
			fwrite( &imageav[x+y*XY_WIDTH_IN], 4, 1, imageavraw);
		}
	}
    fprintf(stderr,"\b\b\b%d%%", 2*100/num_imgs); // Update progress %

	// Amplitude centroided intensity image
	for (y=XY_WIDTH_IN*zoom+1; y>=0; y--) {
		for (x=0; x<XY_WIDTH_IN*zoom+2; x++) {
			fwrite( &imagec[x+y*XY_WIDTH_IN*zoom], 4, 1, imagecraw);
		}
	}
    fprintf(stderr,"\b\b\b%d%%", 3*100/num_imgs); // Update progress %

	// Timecode centroided intensity image
	for (y=XY_WIDTH_OUT-1; y>=0; y--) {
		for (x=0; x<XY_WIDTH_OUT; x++) {
			fwrite( &imaget[x+y*XY_WIDTH_IN], 4, 1, imagetraw);
		}
	}
    fprintf(stderr,"\b\b\b%d%%", 4*100/num_imgs); // Update progress %

    // FLIM image without centroiding
	for (y=XY_WIDTH_IN-1; y>=0; y--) {
		for (x=0;x<XY_WIDTH_IN;x++) {
			for (t=0; t<TIMEBINS_OUT; t++) {
				fwrite(&flimimageav[x+XY_WIDTH_IN*y+XY_WIDTH_IN*XY_WIDTH_IN*t], 2, 1, flimimageavraw);
			}
		}
	}
    fprintf(stderr,"\b\b\b%d%%", 5*100/num_imgs); // Update progress %

    // Centroided FLIM image
	for (y=XY_WIDTH_IN-1; y>=0; y--) {
		for (x=0;x<XY_WIDTH_IN;x++) {
			for (t=0; t<TIMEBINS_OUT; t++) {
				fwrite(&flimimagep[x+XY_WIDTH_IN*y+XY_WIDTH_IN*XY_WIDTH_IN*t], 2, 1, flimimagepraw);
			}
		}
	}
    fprintf(stderr,"\b\b\b%d%%", 6*100/num_imgs); // Update progress %

	printf("\n\n..done writing images.");

	// Print results
	printf("\n");
	printf("\nExperiment duration: %.2f s", totalTimeS);
	printf("\nCount rate: %.0f kHz", event/totalTimeS/1000);
	printf("\nHit rate: %.2f MHz", counter/totalTimeS/1000000);
	printf("\n\nTotal frames: %lld", frame-startFrame);
	printf("\nSkipped frames: %d", skippedFrame);
	printf("\nEmpty frames: %lld", emptyframe);
	printf("\nTotal triggered pixels: %lld", counter);
	printf("\nTotal photon events: %ld", event);
	printf("\nDiscarded undersize: %lu", discarded_under);
	printf("\nDiscarded oversize: %lu", discarded_over);

	printf("\n\nNon-centroided image size = %d x %d.", XY_WIDTH_OUT, XY_WIDTH_OUT);
	printf("\nCentroided image size = %d x %d.", XY_WIDTH_IN*zoom+2, XY_WIDTH_IN*zoom+2);
	printf("\nFLIM image size = %d x %d x %d.", XY_WIDTH_OUT, XY_WIDTH_OUT, TIMEBINS_OUT);
	printf("\n\n");

	// Print to outparam file
    fprintf(outparam, "\nInput file(s):\n");
    for (int curFile=0; curFile<numFiles; curFile++) {
        fprintf(outparam, "%s\n", argv[argd+curFile]);
    }
	fprintf(outparam, "\nPhoton event size: min = %d, max = %d pixels.", minarea, maxarea);
	fprintf(outparam, "\nMagnification for centroiding = %d.", zoom);
	fprintf(outparam, "\nOutput time bin width = %.0f ns.", TIMEBIN_OUT_US*1000);
	fprintf(outparam, "\n");
	fprintf(outparam, "\nFLIM image size = %d x %d x %d.", XY_WIDTH_OUT, XY_WIDTH_OUT, TIMEBINS_OUT);
	fprintf(outparam, "\nCentroided image size = %d x %d.", XY_WIDTH_IN*zoom+2, XY_WIDTH_IN*zoom+2);
	fprintf(outparam, "\n");
	fprintf(outparam, "\nNumber of input files: %d", numFiles);
	fprintf(outparam, "\nTotal processed frames: %lld", frame-startFrame);
	fprintf(outparam, "\nSkipped frames: %d", skippedFrame);
	fprintf(outparam, "\nEmpty frames: %lld", emptyframe);
	fprintf(outparam, "\nTotal triggered pixels: %lld", counter);
	fprintf(outparam, "\nTotal photon events: %ld", event);
	fprintf(outparam, "\nDiscarded undersize: %lu", discarded_under);
	fprintf(outparam, "\nDiscarded oversize: %lu", discarded_over);
	fprintf(outparam, "\n");
	fprintf(outparam, "\nExperiment duration: %.2f s", totalTimeS);
	fprintf(outparam, "\nHit rate: %.0f kHz", counter/totalTimeS/1000);
	fprintf(outparam, "\nCount rate: %.0f kHz", event/totalTimeS/1000);
	fprintf(outparam, "\n\n");

    // Close the output files
	fclose(flimhist);
	fclose(counthist);
	fclose(imageavraw);
	fclose(imagecraw);
	fclose(imagetraw);
	fclose(flimimageavraw);
	fclose(flimimagepraw);
    fclose(eventlist);
    fclose(outparam);

    // Free memory
    free(tab_a);
    free(tab_t);
    free(donetab_a);
    free(donetab_t);
	free(imagec);
	free(imaget);
	free(imageav);
	free(flimhistav);
	free(flimhistp);
	free(flimimageav);
	free(flimimagep);

	return 0;

} // end of main()


/*****************************************************************************************
* make_ics_header
* - Makes a header for the .ics file
*****************************************************************************************/
void make_ics_header(FILE *fp, char *name, int channels, float timebin_out_s, int dimY, int dimX) {
  fprintf(fp, \
"\t\nics_version\t2.0\n\
filename\t%s\n\
layout\tparameters\t4\n\
layout\torder\tbits\tz\tx\ty\n\
layout\tsizes\t16\t%d\t%d\t%d\n\
layout\tcoordinates\tvideo\n\
layout\tsignificant_bits\t16\n\
representation\tformat\tinteger\n\
representation\tsign\tunsigned\n\
representation\tcompression\tuncompressed\n\
representation\tbyte_order\t1\t2\n\
parameter\torigin\t0.000000\t0.000000\t0.000000\t0.000000\n\
parameter\tscale\t1.000000\t%.5e\t1.000000\t1.000000\n\
parameter\tunits\trelative\tns\tpixels\tpixels\n\
parameter\tlabels\tintensity\tmicro-time\tx-position\ty-position\n",\
  name, channels, dimY, dimX, timebin_out_s*1000000000);

  fprintf(fp,\
"history\tsoftware\ttimepix3cam version %s\n\
history\texperimenter\tLiisa Hirvonen\n\
history\tmicroscope\tTimepix3cam @ UCC\n\
history\tcreation date\t%s %s\n\
history\ttype\tLifetime\n\
history\tlabels\tt x y\n\
history\tdimensions\t%d %d %d\n\
history\toffsets\t0 0 0\n\
history\tunits\ts m m\n\
history\textents\t%.5e %d %d\n\
end\t\n", \
  VERSION, __TIME__, __DATE__, channels, dimY, dimX, timebin_out_s*channels, dimY, dimX);
}


/*****************************************************************************************
* sigproc
* - For termination of running program
*****************************************************************************************/
void sigproc(int n) {

    char res;
    printf("\n\tCtrl-C pressed, are you sure you want to quit the program?\n"\
           "\tProcessing done so far will be saved.\n"\
           "\tPress Y to abort, anything else to continue.");

    scanf("%c", &res);

    if (res == 'Y' ) {
        printf("Aborting...\n");
        abort_received = 1;
    }
}


/*****************************************************************************************
* FindSurroundingPixels
* - Recursive function to find all triggered pixels in a cluster
*****************************************************************************************/
void FindSurroundingPixelsAmplitude(unsigned short x, unsigned short y,
                           unsigned short dimx, unsigned short *p_area, unsigned int *p_xc,
                           unsigned int *p_yc, unsigned short *p_xp, unsigned short *p_yp,
                           unsigned int *p_intensity, unsigned int *p_maxint, bool p_donetab[],
                           const unsigned int p_tab[], unsigned int *maxarea) {

    // Mark this pixel as checked
    *(p_donetab + x+y*dimx) = 1;

    // If the event is not too big
    if (*p_area <= *maxarea ) {
        *p_intensity += *(p_tab+x+y*dimx); // Add the pixel intensity to event intensity
        *p_area += 1;   // Add 1 to event area
        *p_xc += x* *(p_tab+x+y*dimx);  // For the centre position, multiply the pixel
        *p_yc += y* *(p_tab+x+y*dimx);  // coordinates by intensity and add to the position
    }

    // If it is brighter than other pixels in this event area
    if ( *(p_tab+x+y*dimx) > *p_maxint ) {
        *p_maxint = *(p_tab+x+y*dimx);  // Update maximum intensity
        *p_xp = x;  // Mark this as the brightest pixel
        *p_yp = y;
    }

    // Check every surrounding pixel. If they have not been checked and are not zero, call this function again.
    if ( *(p_donetab + x+1+y*dimx)==0 && *(p_tab+x+1+y*dimx) !=0 )
        FindSurroundingPixelsAmplitude(x+1, y, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
    if ( *(p_donetab + x-1+y*dimx)==0 && *(p_tab+x-1+y*dimx) !=0 )
        FindSurroundingPixelsAmplitude(x-1, y, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
    if ( *(p_donetab + x+(y+1)*dimx)==0 && *(p_tab+x+(y+1)*dimx) !=0 )
        FindSurroundingPixelsAmplitude(x, y+1, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
    if ( *(p_donetab + x+(y-1)*dimx)==0 && *(p_tab+x+(y-1)*dimx) !=0 )
        FindSurroundingPixelsAmplitude(x, y-1, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
}

void FindSurroundingPixelsTime(unsigned short x, unsigned short y,
                           unsigned short dimx, unsigned short *p_area, unsigned int *p_xc,
                           unsigned int *p_yc, unsigned short *p_xp, unsigned short *p_yp,
                           unsigned int *p_intensity, unsigned int *p_maxint, bool p_donetab[],
                           const unsigned int p_tab[], unsigned int *maxarea) {

    // Mark this pixel as checked
    *(p_donetab + x+y*dimx) = 1;

    // If the event is not too big
    if (*p_area <= *maxarea ) {
        *p_intensity += *(p_tab+x+y*dimx); // Add the pixel intensity to event intensity
        *p_area += 1;   // Add 1 to event area
        *p_xc += x* *(p_tab+x+y*dimx);  // For the centre position, multiply the pixel
        *p_yc += y* *(p_tab+x+y*dimx);  // coordinates by intensiy and add to the position
    }

    // If it is brighter than other pixels in this event area
    // NOTE: here we want the EARLIEST timecode, it is now finding the highest?
    if ( *(p_tab+x+y*dimx) < *p_maxint ) {
        *p_maxint = *(p_tab+x+y*dimx);  // Update maximum intensity
        *p_xp = x;  // Mark this as the brightest pixel
        *p_yp = y;
    }

    // Check every surrounding pixel. If they have not been checked and are not zero, call this function again.
    if ( *(p_donetab + x+1+y*dimx)==0 && *(p_tab+x+1+y*dimx) !=0 )
        FindSurroundingPixelsTime(x+1, y, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
    if ( *(p_donetab + x-1+y*dimx)==0 && *(p_tab+x-1+y*dimx) !=0 )
        FindSurroundingPixelsTime(x-1, y, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
    if ( *(p_donetab + x+(y+1)*dimx)==0 && *(p_tab+x+(y+1)*dimx) !=0 )
        FindSurroundingPixelsTime(x, y+1, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
    if ( *(p_donetab + x+(y-1)*dimx)==0 && *(p_tab+x+(y-1)*dimx) !=0 )
        FindSurroundingPixelsTime(x, y-1, dimx, p_area, p_xc, p_yc, p_xp, p_yp,
                              p_intensity, p_maxint, p_donetab, p_tab, maxarea);
}


/*****************************************************************************************
* init
* - Initialise: allocate memory, open result files for writing
*****************************************************************************************/
void init(char *argv[], int argd, int zoom,
          FILE **p_eventlist, FILE **p_outparam, // Text files
          FILE **p_imageavraw, FILE **p_imagecraw, FILE **p_imagetraw, // Intensity image files
          FILE **p_flimimageavraw, FILE **p_flimimagepraw, // FLIM image files
          FILE **p_flimhist, FILE **p_counthist, // Histogram files for counts & time decay
          unsigned int **p_tab_a, unsigned int **p_tab_t, // Arrays for current frame
          bool **p_donetab_a, bool **p_donetab_t,   // Arrays for current frame
          unsigned long int **p_imageav, unsigned int **p_imagec, unsigned int **p_imaget, // Arrays for intensity images
          unsigned int **p_flimimageav, unsigned int **p_flimimagep, // Arrays for FLIM images
          unsigned int **p_flimhistav, unsigned int **p_flimhistp) {

    fprintf(stderr,"\nInitialisation... ");

    char name[512], fname[512];

    /*  tab_a is the amplitude image in a 1 dimensional array,
        after the frame events have been read in. index tab[x+y*dimx] */
    *p_tab_a=(unsigned int*)calloc ((XY_WIDTH_OUT)*(XY_WIDTH_OUT),sizeof(unsigned int));
    if (*p_tab_a==NULL) {
        printf("Allocation error: tab_a\n");
        exit(ERR_ALLOC);
    }

    /*  tab_t is the timecode image in a 1 dimensional array,
        after the frame events have been read in. index tab[x+y*dimx] */
    *p_tab_t=(unsigned int*)calloc ((XY_WIDTH_OUT)*(XY_WIDTH_OUT),sizeof(unsigned int));
    if (*p_tab_t==NULL) {
        printf("Allocation error: tab_t\n");
        exit(ERR_ALLOC);
    }

    /*  donetab_a keeps track of amplitude image pixels that have been checked,
        donetab[x+y*dimx]=0 if the pixel hasn't be counted yet, and 1 otherwise */
    *p_donetab_a=(bool*)calloc ((XY_WIDTH_OUT+1)*(XY_WIDTH_OUT+1), sizeof(bool));
    if (*p_donetab_a==NULL) {
        printf("Allocation error: donetab_a\n");
        exit(ERR_ALLOC);
    }

    /*  donetab_t keeps track of timecode image pixels that have been checked,
        donetab[x+y*dimx]=0 if the pixel hasn't be counted yet, and 1 otherwise */
    *p_donetab_t=(bool*)calloc ((XY_WIDTH_OUT+1)*(XY_WIDTH_OUT+1), sizeof(bool));
    if (*p_donetab_t==NULL) {
        printf("Allocation error: donetab_t\n");
        exit(ERR_ALLOC);
    }

	// Allocate memory for sum image array
	*p_imageav = (unsigned long int*) calloc ((XY_WIDTH_OUT)*(XY_WIDTH_OUT), sizeof(unsigned long int));
	if (*p_imageav==NULL) {
		printf("Allocation error: imageav\n");
		exit(ERR_ALLOC);
	}

	// Allocate memory for amplitude centroided image array
	*p_imagec = (unsigned int*) calloc ((XY_WIDTH_IN*zoom+2)*(XY_WIDTH_IN*zoom+2), sizeof(unsigned int));
	if (*p_imagec==NULL) {
		printf("Allocation error: imagec\n");
		exit(ERR_ALLOC);
	}

	// Allocate memory for timecode centroided image array
	*p_imaget = (unsigned int*) calloc ((XY_WIDTH_OUT)*(XY_WIDTH_OUT), sizeof(unsigned int));
	if (*p_imaget==NULL) {
		printf("Allocation error: imaget\n");
		exit(ERR_ALLOC);
	}

	// Allocate memory for FLIM sum histogram
	*p_flimhistav = (unsigned int*) calloc (TIMEBINS_OUT, sizeof(unsigned int));
	if (*p_flimhistav==NULL) {
		printf("Allocation error: flimhistav\n");
		exit(ERR_ALLOC);
	}

	// Allocate memory for FLIM centroided histogram
	*p_flimhistp = (unsigned int*) calloc (TIMEBINS_OUT, sizeof(unsigned int));
	if (*p_flimhistp==NULL) {
		printf("Allocation error: flimhistp\n");
		exit(ERR_ALLOC);
	}

	// Allocate memory for FLIM sum data cube
	*p_flimimageav = (unsigned int*) calloc ((XY_WIDTH_OUT)*(XY_WIDTH_OUT)*TIMEBINS_OUT, sizeof(unsigned int));
	if (*p_flimimageav==NULL) {
		printf("Allocation error: flimimageav\n");
		exit(ERR_ALLOC);
	}

	// Allocate memory for FLIM centroided data cube
	*p_flimimagep = (unsigned int*) calloc ((XY_WIDTH_OUT)*(XY_WIDTH_OUT)*TIMEBINS_OUT, sizeof(unsigned int));
	if (*p_flimimagep==NULL) {
		printf("Allocation error: flimimagep\n");
		exit(ERR_ALLOC);
	}

	// Open sum image file
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
    strcat(name, "_image-sum.raw");
	if((*p_imageavraw = fopen(name, "w")) == NULL){
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
	}

	// Open amplitude-centroided image file
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_image-centroided-ampl.raw");
	if((*p_imagecraw = fopen(name, "w")) == NULL){
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
	}

	// Open timecode-centroided image file
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_image-centroided-time.raw");
	if((*p_imagetraw = fopen(name, "w")) == NULL){
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
	}

	// Open FLIM sum image file
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_flim-image-sum.ics");
	if((*p_flimimageavraw = fopen(name, "w")) == NULL){
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
	}
	// Make header for .ics file, total time in seconds
	make_ics_header(*p_flimimageavraw, name, TIMEBINS_OUT, timebin_out_s, XY_WIDTH_IN, XY_WIDTH_IN);

	// Open FLIM centroided image file
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_flim-image-centr.ics");
	if((*p_flimimagepraw = fopen(name, "w")) == NULL){
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
	}
	// Make header for .ics file, total time in seconds
	make_ics_header(*p_flimimagepraw, name, TIMEBINS_OUT, timebin_out_s, XY_WIDTH_IN, XY_WIDTH_IN);

	// Open decay histogram file & print header
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_decay-histogram.dat");
    if ((*p_flimhist=fopen(name,"w"))==NULL) {
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
    }
    fprintf(*p_flimhist,"Time (us)\tPhotons sum\tPhotons centr\n");

	// Open count histogram file & print header
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_count-histogram.dat");
    if ((*p_counthist=fopen(name,"w"))==NULL) {
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
    }
    fprintf(*p_counthist,"Frame\tAvg n\tAvg t\n");

	// Open eventlist file & print header
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_res-eventlist.dat");
    if ((*p_eventlist=fopen(name,"w"))==NULL) {
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
    }
    fprintf(*p_eventlist,"#frame\txp\typ\txc\tyc\tarea\tintensity\n");

    // Open outparam file & print header
	strcpy(fname, argv[argd]);
    memset(name,0,sizeof(name));
    strncpy(name,fname,strlen(fname)-4);
	strcat(name, "_res-outparam.txt");
    if ((*p_outparam=fopen(name, "w"))==NULL){
        printf("\nCould not open %s for writing\n", name);
        exit(ERR_OUTFILE);
    }
    fprintf(*p_outparam,"Processed with timepix3cam version %s,\ncompiled on %s at %s,\nfrom %s \n",
           VERSION, __DATE__, __TIME__, __FILE__);

    fprintf(stderr,"done.");

}


/*****************************************************************************************
* printUsage
* - Prints info about required input parameters
*****************************************************************************************/
void printUsage() {
//	printf("\nInput arguments: <minarea> <maxarea> <zoom> <last frame> <input file 1> <input file 2> ...");
	printf("\nInput arguments: <last frame> <input file 1> <input file 2> ...");
	//printf("\n<minarea> = Minimum area in pixels for a photon event.");
	//printf("\n<maxarea> = Maximum area in pixels for a photon event.");
	//printf("\n<zoom> = Magnification for sub-pixel centroiding.");
	printf("\n<last frame> = Number of frames to be processed, enter 0 to process all.");
	printf("\n<input file> = Full path to the input .csv file. Enter multiple file names if input is split.");
    printf("\n\n");
	exit(ERR_USAGE);
}
