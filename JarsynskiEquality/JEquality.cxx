/*
I used this script to calculate the free energy change of a biophysical simulation 
of the nanoparticle going through the bilayer, using the Jarsynski equality and the steered 
molecular dynamics simulation. The simulation protocol is following the LAMMPS simulation package, 
specifically the page:

http://lammps.sandia.gov/doc/fix_smd.html

The approach I used was to the following: 


 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

int main(int argv, char** argc) {
  
  // When compiling, it is useful to type in the following commmand to prevenet
  // segmentation faults
  
  /*ulimit -s unlimited*/

  const int numberoftrajectories = 2; /* The number of pull files we have. I have set it here as 2 but let's say if you have 50 files, then you can set it to 50. You need to index each pull file at the end with the appropriate number - pull.1, pull.2, pull.3 .. etc  */
  const int numbins = 191; // The approximate number of bins I divided the z coordinate of the nanoparticle 

  double binmin[numbins];  // minimum value in each bin
  double binmax[numbins]; // maximum value in each bin

  double variance[numbins]; // bin to append the variance
  int bin_counter[numbins]; // bin to append the counter for averaging purposes  
  
  double minz = -95.5;
  double binwidth = 1.0;
  double averageforce[numbins]; /*The average force*/
  double energysum[numbins];
  double energydif[numbins];
  double work1[numbins];
  double newav[numbins];

  //    double num[nlines][numberoftrajectories];
  //    double cz[nlines][numberoftrajectories]; /*ditto*/
  //    double Un[nlines][numberoftrajectories];
  //    double force[nlines][numberoftrajectories];
  //    double work[nlines][numberoftrajectories];

  double boltzmann = 0.0019872041; /* Boltzmann constant units are kcal */
  char pull[40]; /* buffer which cocatenates the string of filename and index */
  char filename[40] = "pull";
  char line[100];

  // -- random integers, mostly used for for loops and other random stuff

  int i;
  long int j;
  int k;
  
  int index; 
  char test_char;
  float deltaG;
  double beta = 1 / (298.0 * boltzmann); /* boltzmann factor */


  /*--- input parameters to read from each force/work file produced from the SMD run ---*/

  double f; /*force*/
  double w; /*work*/
  double z; /*z coordinate*/
  int number; /*number*/
  double U; /*second column*/
  FILE *ipf; /* input file */
  FILE *new_file;

  /*Set bin_counter to 0*/
  
  for (i = 0; i < numbins; i++) {
    bin_counter[i] = 0;
  }

  // memory allocation is required here as the work/force file will be of different lengths - Hence, the
  // total number of line need to be read, allocated to the heap and THEN processed. 
   
  int **num;
  double **cz;
  double **Un;
  double **force;
  double **work;

  // The number of trajectories (or the number of lines, poor wording I must admit from my part..) now is known,
  // and hence the stack allocation can be commenced. 
  
  num = malloc(sizeof (int *) * numberoftrajectories);
  cz = malloc(sizeof (double *) * numberoftrajectories);
  Un = malloc(sizeof (double *) * numberoftrajectories);
  force = malloc(sizeof (double *) * numberoftrajectories);
  work = malloc(sizeof (double *) * numberoftrajectories);

  // Dangling pointer error
  
  if (num == NULL) {
    printf("Cannot allocate memory\n");
    exit(1);
  }
  
  // Cocatenate the pull file names over a index i (so we will get pull.i) and then we can loop over all the work/force
  // files
  
  for(i=0; i<numberoftrajectories; i++) {
    index = i + 1;
    sprintf(pull, "%s.%d", filename, index);
    if((ipf = fopen(pull, "r")) == NULL) {
      printf("Cannot open files\n");
      exit(1);
    } else {
      printf("Opening file: %d\n",index);
    }
    //how many lines do I have
    int nlines = 0;

    while (fgets(line, sizeof (line), ipf) != NULL) {
      nlines++;
    }
    close(ipf);

    printf("Number of lines: %d\n",nlines);

    // allocating memory

    num[i] = malloc(nlines * sizeof (int));
    cz[i] = malloc(nlines * sizeof (double));
    Un[i] = malloc(nlines * sizeof (double));
    force[i] = malloc(nlines * sizeof (double));
    work[i] = malloc(nlines * sizeof (double));  

    if (num[i] == NULL || cz[i] == NULL || Un[i] == NULL || force[i] == NULL || work[i] == NULL) {
      printf("Unable to assign line memory\n");
      exit(1);
    } 
    //  read new_file
    
    new_file = fopen(pull,"r");

    for (j=0;j<nlines;j++)  {
      
       if (fgets(line,sizeof(line),new_file) != NULL) { 
       	   test_char = line[0]; 
           if((test_char != '#') && (test_char != ';')) { 
               	sscanf(line,"%d %lf %lf %lf %lf",&number,&z,&U,&f,&w);
       	   } 
        }
        
      	num[i][j] = number;
      	cz[i][j] = z;
      	Un[i][j] = U;
      	force[i][j] = f;
      	work[i][j] = w;


       for(k=0;k<numbins;k++) {  
           binmin[k] = k*binwidth + minz; 
           binmax[k] = k*binwidth + minz + binwidth; 


       	   if ((cz[i][j] > binmin[k]) &&  (cz[i][j] <= binmax[k])) { 
	     bin_counter[k] = (bin_counter[k])+1;  
                
       	       if (bin_counter[k] == 2) {                                       
		 work1[k] = work[i][j];
		 //if ((cz[i][j] > -0.500000) &&  (cz[i][j] <= 0.500000)) {
		 //    printf("%f %f %d\n",cz[i][j],work[i][j],bin_counter[k]);
		 // }    
       	       } 
       	       energydif[k] = (-beta*(work[i][j]-work1[k]));
	       energysum[k] += exp(-beta*(work[i][j]-work1[k]));
	       //printf("z:%lf work:%lf work1:%lf beta*(work-work1):%lf sum of exponential:%lf \n",cz[i][j],work[i][j],work1[k],energydif[k],energysum[k],bin_counter[k]); 
           }
       } 
    }
    fclose(new_file);
  }

  //outside trajectory file loop    
  
  for (k=0;k<numbins;k++) {

    // Processing the data using the Jarysnki equality
    
    deltaG = (-1/beta) * ( (-beta*work1[k]) + log(1+energysum[k]) ) - log(bin_counter[k]);

    printf("%lf %lf %lf %d %lf \n",binmax[k],deltaG,energysum[k],bin_counter[k],(exp(-beta*work1[k])));
    
  }

  int v=0;

  // free memory  allocation of each double pointer
  for(v=0; v<numberoftrajectories; v++) {
    free(num[v]);
    free(cz[v]);
    free(Un[v]);
    free(work[v]);
    free(force[v]);
  }
  
  free(num);
  free(cz);
  free(Un);
  free(work);
  free(force);

  return (0);
}

