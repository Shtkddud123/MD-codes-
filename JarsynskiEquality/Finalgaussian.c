#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

int main() {
  /*ulimit -s unlimited*/
  const int numberoftrajectories = 2; /* The number of pull files we have */
  const int numbins = 191;
  double binmin[numbins];
  double binmax[numbins];
  double variance[numbins];
  int bin_counter[numbins];
  double minz = -95.5;
  double binwidth = 1.0;
  double averageforce[numbins]; /*The average force*/
  double energysum[numbins];
  double energydif[numbins];
  double work1[numbins];
  double newav[numbins];
  double blockaverage; 
  //    double num[nlines][numberoftrajectories];
  //    double cz[nlines][numberoftrajectories]; /*ditto*/
  //    double Un[nlines][numberoftrajectories];
  //    double force[nlines][numberoftrajectories];
  //    double work[nlines][numberoftrajectories];
  double boltzmann = 0.0019872041; /*units are kcal*/
  char pull[40]; /*buffer which cocatenates the string of filename and index*/
  char filename[40] = "pull";
  char line[100];
  int i;
  long int j;
  int k;
  int index; /*index for the blah blah blah*/
  char test_char;
  float deltaG;
  double beta = 1 / (298.0 * boltzmann); /*boltzmann factor*/


  /*--------------------------------------------*/

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

  int **num;
  double **cz;
  double **Un;
  double **force;
  double **work;
  num = malloc(sizeof (int *) * numberoftrajectories);
  cz = malloc(sizeof (double *) * numberoftrajectories);
  Un = malloc(sizeof (double *) * numberoftrajectories);
  force = malloc(sizeof (double *) * numberoftrajectories);
  work = malloc(sizeof (double *) * numberoftrajectories);
  
  if (num == NULL) {
    printf("Cannot allocate memory\n");
    exit(1);
  } 

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

    //allocating memory
    num[i] = malloc(nlines * sizeof (int));
    cz[i] = malloc(nlines * sizeof (double));
    Un[i] = malloc(nlines * sizeof (double));
    force[i] = malloc(nlines * sizeof (double));
    work[i] = malloc(nlines * sizeof (double));  

    if (num[i] == NULL || cz[i] == NULL || Un[i] == NULL || force[i] == NULL || work[i] == NULL) {
      printf("Unable to assign line memory\n");
      exit(1);
    } 

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

     
	 /*We are defining the block width*/

	 binmin[k] = k*binwidth + minz; 
	 binmax[k] = k*binwidth + minz + binwidth;

	 /* Creating a method to place each cz value in to the appropriate bin */ 
	 	 
	 if ((cz[i][j] > binmin[k]) &&  (cz[i][j] <= binmax[k])) {
 
	   bin_counter[k] = (bin_counter[k])+1;  
           
	   /*we are starting when the bin counter is 2 because...*/ 
	   
	   if (bin_counter[k] == 2) {
                                       
	     work1[k] = work[i][j];
	     //if ((cz[i][j] > -0.500000) &&  (cz[i][j] <= 0.500000)) {
	     //    printf("%f %f %d\n",cz[i][j],work[i][j],bin_counter[k]);
	     // }    
	   } 

	   /* The sum and energy difference*/

	   energydif[k] = (-beta*(work[i][j]-work1[k]));
	   energysum[k] += exp(-beta*(work[i][j]-work1[k]));
	   //printf("z:%lf work:%lf work1:%lf beta*(work-work1):%lf sum of exponential:%lf \n",cz[i][j],work[i][j],work1[k],energydif[k],energysum[k],bin_counter[k]); 

	 }
       }

       /* We are now describing the data between pull.1 and pull.5*/
    }
    fclose(new_file);
  }
  
  //outside trajectory file loop    
  
  for (k=0;k<numbins;k++) {
     
    deltaG = (-1/beta) * ( (-beta*work1[k]) + log(1+energysum[k]) ) - log(bin_counter[k]);
    printf("%lf %lf %lf %d %lf \n",binmax[k],deltaG,energysum[k],bin_counter[k],(exp(-beta*work1[k])));
    
  }
  
  int v=0;
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

