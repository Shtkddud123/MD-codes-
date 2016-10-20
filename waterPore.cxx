#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>

double LJ (double dist_ab, double eps_ab, double sigma_ab); 
double distance(double ax,double ay,double az, double bx, double by, double bz);
double minDist(double ax,double ay,double az, double bx, double by, double bz,
	       double bxx, double byy, double bzz);
void minimum_Image(double *a,double *b,double *c,double bxx,double byy,double bzz);
double centreOfMass(int Oa, int H1a, int H2a,int Ob, int H1b, int H2b,int Oc, int H1c, int H2c);

int main (int argc, char **argv) 

{ 

  /*             b (H2O)           */
  /*  a (Au)         c (H2O)       */
  /*             d (H2O)          */


  /*Declaration of variables*/ 

  const int numberofatoms = 71313;   /* Declaring the paramters for the atoms in the box */   
  const int numberof7 = 2000;   /* The number of group 7 beads in each trajectory */
  const int numberof6 = 4000;   /* The number of group 6 beads in each trajectory */

  /*Parameters to loop over an entire dump file*/ 
  int firstnumberoftraj = 1; /*The number of trajectories in the first dump file*/
  int numberoftraj = 100; /*The number of trajectories in the dump file*/
  int trajno = 0; /*the nth traj we are at*/
  int inittrajno = 0; 
  
  
  int boxdim = 3; 
  const int nbins = 100;

  int numberofwaters = 1090; /*total number of water triads*/
  double x,y,z; /*coordinates for the atoms in the box*/
  double xco[numberofatoms],yco[numberofatoms],zco[numberofatoms]; 
  int a[numberofatoms],b[numberofatoms];

  /* bottom surface of the bilayer */
  
  double botfirstbilayer[numberoftraj][numberof7];
  int botfirstindex[numberoftraj][numberof7];

  /* top surface of the bilayer */
  double topfirstbilayer[numberoftraj][numberof7];
  int topfirstindex[numberoftraj][numberof7]; 
  
  /* top waters */
  double topfirstwater[numberoftraj][numberofatoms];
  int topfirstwaterindex[numberoftraj][numberofatoms]; 
  
  /* bot waters */
  double botfirstwater[numberoftraj][numberofatoms];
  int botfirstwaterindex[numberoftraj][numberofatoms];

  /*------------------------------------*/

  /* bottom surface of the bilayer */
  double botbilayer[numberoftraj][numberof7];
  int botindex[numberoftraj][numberof7];

 /* top surface of the bilayer */
  double topbilayer[numberoftraj][numberof7];
  int topindex[numberoftraj][numberof7]; 

  /* top waters */
  double topwater[numberoftraj][numberofatoms];
  int topwaterindex[numberoftraj][numberofatoms]; 

/* bot waters */

  double botwater[numberoftraj][numberofatoms];
  int botwaterindex[numberoftraj][numberofatoms];
  int atomtype; /*atom type; 1,2 = water. 3 = gold */ 

  /*box dimensions*/

  double box1;
  double box2; 
  double boxlength[boxdim];


  /*We are now declaring the atom coordinates of each individual atom in the water triad*/ 

  double Obcoordinatex; /*cartesian coordinates for the oxygen on water b */ 
  double Obcoordinatey; /*cartesian coordinates for the oxygen on water b */ 
  double Obcoordinatez; /*cartesian coordinates for the oxygen on water b */ 
  
  double Occoordinatex; /*cartesian coordinates for the oxygen on water c */ 
  double Occoordinatey; /*cartesian coordinates for the oxygen on water c */
  double Occoordinatez; /*cartesian coordinates for the oxygen on water c */
  
  double Odcoordinatex; /*cartersian coordinates for the oxygen on water d */
  double Odcoordinatey; /*cartersian coordinates for the oxygen on water d */
  double Odcoordinatez; /*cartersian coordinates for the oxygen on water d */

  /*Declaring the water triplet*/ 
  int w1[numberofwaters]; /*water under consideration*/ 
  int w2[numberofwaters]; /*index of the closest water*/
  int w3[numberofwaters]; /*index of the second closest water*/ 
  int watercounter[numberofwaters]; /*counting how many times the water molecule in question appears in a blob*/
  // double distance[numberofwaters]; 
  double boxLJ; 
  double weight; 
 
  /*We are now declaring the COM for each water molecule in the triad*/
  double COMbx; /*COM of b*/ 
  double COMby; /*COM of b*/
  double COMbz; /*COM of b*/
  double COMsetx;
  double COMsety; 
  double COMsetz;
  double distgoldoxygenb; /*distance between the centre of mass of each water triad and the gold particle*/ 
  double distgoldoxygenc; /*distance between the centre of mass of each water triad and the gold particle*/ 
  double distgoldoxygend; /*distance between the centre of mass of each water triad and the gold particle*/ 
  double COMdistgoldwater;
  double xmass; /*COM of the triad*/
  double ymass; /*COM of the triad*/ 
  double zmass; /*COM of the triad*/
  double COMdist; /*printing out the distance between each water COM*/ 
  double minimumCOMdist; /*prints out the smallest distances between each COM to determine which water molecules to make into a triad*/ 
  
  /*LJ parameters*/ 
  double Oepsilon = 0.15207;
  double Osigma = 3.1507;
  double Hepsilon = 0.000;
  double Hsigma = 0.0;
  double goldsigma = 2.629; 
  double goldepsilon = 5.29;
  double COMcrosssigma; 
  double COMcrossep; 
  double COMgoldCOMcrosssigma; 
  double COMgoldCOMcrossep; 
  double V_ave[nbins]; /*The average interaction energy*/ 
  double sumofweight[nbins]; /*sum of weights*/ 
  double LJb;
  double LJc; 
  double LJd; 
    
  /*We are now declaring the distance between c -- b and c --d */
  double cbdist; 
  double cddist; 
  
  /*Defining each water molecule in the list*/ 
  int firstoxygen; 
  int firsthydrogen;
  int lasthydrogen; 
  int newoxygen;
  int newfirsthydrogen;
  int newlasthydrogen; 

  /*Keeping track of the closest waters*/ 
  double closestdistance1; 
  int closestwater1; 
  double closestdistance2; 
  int closestwater2; 
  double bignumber = 100000000.0; 

  /*misc*/
  char line[100];
  int n = 0; 
  int m = 0;
  int index;
  int i; 
  int j; 
  int k = 0;
  int l = 0;
  int w;
  int u; 
  double r; 
  double xdiff;
  double ydiff;
  double zdiff; 

  /*Parameters for bin*/
  
  double V_average[nbins]; 
  double weights[nbins]; 
  double bin_min[nbins];
  double bin_max[nbins]; 
  int bin_counter[nbins]; 
  double newLJ[nbins];

  /*Counter for the atom 7 headgroups*/
  int botfirstcounter = 0;
  int topfirstcounter = 0;
  
  /* Counter fot the waters*/
  int topfirstwatercounter = 0;
  int botfirstwatercounter = 0;

  /*Counter for the atom 7 headgroups*/
  int botcounter = 0;
  int topcounter = 0; 

  /* Counter fot the waters*/
  int topwatercounter = 0; 
  int botwatercounter = 0;

  // int ww[numberofwaters];

  /*-----------------------------------------------------------------------------*/
  /* We need to get the coordinates from trajectory 0 from whatever dump file we furst used to record the trajectories over time */

  /*Start of code*/ 
  /*Reading in the dump file*/ 

  printf("test \n"); 
 FILE *newfile; /* input file */
  
  /* open file for reading */
  printf("Opening file %s\n","dump.nano1.5");
  newfile = fopen("dump.nano1.5", "r");
  
  /* check for error opening file */
    
  if(newfile == NULL) {
    
    printf("Error opening file\n");
    exit(1);
  }
  

  // printf("test \n"); 

  /* get a line from the file */
  /* fgets() returns NULL when it reaches an error or end of file */ 


  /*Refreshing the bin!*/ 

  double binwidth = 0.25;  

  for (j=0;j<nbins;j++) { 

    V_ave[j] = 0;
    sumofweight[j] = 0; 
    k = 0;
    
  } 
 
  int nlines = numberofatoms + 9;  

  for(inittrajno=0;inittrajno<firstnumberoftraj;inittrajno++) {
  
    //printf("\n");
    // printf("This is the data for trajectory no %d \n", trajno);

    l = 0;
    n = 0;
 

    for(k=0;k<nlines;k++) {
      
      fgets(line,sizeof(line),newfile);

      //while (fgets(line,sizeof(line),ipf) != NULL) {

	 
      if (l < 5) {
	
	/*We are doing nothing*/
	
	//	printf("%s l= %d", line,l);
	
	l++;
	  
      }
      else if ((l > 4 && l < 8)) {
      
	/*We are scanning the bit with just the box parameters*/
      
	sscanf(line, "%lf %lf", &box1, &box2);
	//	printf("%lf %lf \n",box1,box2 );
	boxlength[l-5] = box2-box1;
	l++;
      }
	
      else if (l == 8) {
      
	//	printf(" **** l= %d \n",l);
	/*We are doing nothing*/
      
	l++;
      }
    
      else {

	//	printf(" ***** l = %d \n",l );
	/* convert the text to numbers */
	sscanf(line,"%d %d %lf %lf %lf",&index,&atomtype,&x,&y,&z);
	a[index-1] = index;
	b[index-1] = atomtype;
	xco[index-1] = x*boxlength[0];
	yco[index-1] = y*boxlength[1];
	zco[index-1] = z*boxlength[2];
	n++;
	l++;
	
	
      }

      //printf("%d %d %lf %lf %lf \n", a[3270],b[3270],xco[3270],yco[3270],zco[3270]);

    }

    /*This for-loop prints out all the trajectory coordinates of the entire file, without any condition whatsoever*/    

    /* For the bottom surface*/

    for(i=0;i<numberofatoms;i++) {
      
      if (b[i] == 7 && zco[i] < 100) {

        botfirstbilayer[inittrajno][botfirstcounter] = zco[i];
	botfirstindex[inittrajno][botfirstcounter] = a[i];
	botfirstcounter++;

      }
      
      /*For the top surface*/
     
      if (b[i] == 7 && zco[i] > 100) {

        topfirstbilayer[inittrajno][topfirstcounter] = zco[i];
	topfirstindex[inittrajno][topfirstcounter] = a[i];
	topfirstcounter++;
	
      }

      /*For the top water*/


      if (b[i] == 1 && zco[i] > 100) {

        topfirstwater[inittrajno][topfirstwatercounter] = zco[i];
	topfirstwaterindex[inittrajno][topfirstwatercounter] = a[i];
	topfirstwatercounter++;
      }

      if (b[i] == 1 && zco[i] < 100) {

        botfirstwater[inittrajno][botfirstwatercounter] = zco[i];
	botfirstwaterindex[inittrajno][botfirstwatercounter] = a[i];
	botfirstwatercounter++;
      }
      
    }
   
  }

  fclose(newfile);


  /*----------------------------------------------------------------------------------------------------------------*/

  /*----------------------------------------------------`------------------------------------------------------------*/
  
  /*----------------------------------------------------------------------------------------------------------------*/
  FILE *ipf; /* input file */
  
  /* open file for reading */
  printf("Opening file %s\n","dump.nano1.5");
  ipf = fopen("dump.nano1.5", "r");
  
  /* check for error opening file */
    
  if(ipf == NULL) {
    
    printf("Error opening file\n");
    exit(1);
  }
  
  /* get a line from the file */
  /* fgets() returns NULL when it reaches an error or end of file */ 


  /*Refreshing the bin!*/ 

  for (j=0;j<nbins;j++) { 

    V_ave[j] = 0;
    sumofweight[j] = 0; 
    k = 0;
    
  } 
 
  for(trajno=0;trajno<numberoftraj;trajno++) {

    //  printf("\n");
    // printf("This is the data for trajectory no %d \n", trajno); 

    l = 0;
    n = 0;
 

    for(k=0;k<nlines;k++) { 
      
      fgets(line,sizeof(line),ipf);

      //while (fgets(line,sizeof(line),ipf) != NULL) {

	 
      if (l < 5) {
	
	/*We are doing nothing*/
	
	//	printf("%s l= %d", line,l);
	
	l++;
	  
      }
      else if ((l > 4 && l < 8)) {
      
	/*We are scanning the bit with just the box parameters*/
      
	sscanf(line, "%lf %lf", &box1, &box2);
	//	printf("%lf %lf \n",box1,box2 );
	boxlength[l-5] = box2-box1; 
	l++;
      }
	
      else if (l == 8) {
      
	//	printf(" **** l= %d \n",l);
	/*We are doing nothing*/
      
	l++;
      }
    
      else {

	//	printf(" ***** l = %d \n",l );
	/* convert the text to numbers */
	sscanf(line,"%d %d %lf %lf %lf",&index,&atomtype,&x,&y,&z);
	a[index-1] = index;
	b[index-1] = atomtype;
	xco[index-1] = x*boxlength[0]; 
	yco[index-1] = y*boxlength[1];
	zco[index-1] = z*boxlength[2];  
	n++;
	l++;
	
	
      }

      //printf("%d %d %lf %lf %lf \n", a[3270],b[3270],xco[3270],yco[3270],zco[3270]); 

    }

    /*This for-loop prints out all the trajectory coordinates of the entire file, without any condition whatsoever*/

 
      
    /* For the bottom surface*/ 

    for(i=0;i<numberofatoms;i++) { 
      
      if (b[i] == 7 && zco[i] < 100) { 

        botbilayer[trajno][botcounter] = zco[i]; 
	botindex[trajno][botcounter] = a[i];

	//printf("%d %d %lf %lf %lf %d \n", a[i],b[i],xco[i],yco[i],zco[i],counter);

	botcounter++;

      }
      
      /*For the top surface*/
     
      if (b[i] == 7 && zco[i] > 100) { 

        topbilayer[trajno][topcounter] = zco[i];
	topindex[trajno][topcounter] = a[i];

	//printf("%d %d %lf %lf %lf %d \n", a[i],b[i],xco[i],yco[i],zco[i],counter);

	topcounter++;
	
      }

      /*For the top water*/


      if (b[i] == 1 && zco[i] > 100) { 

        topwater[trajno][topwatercounter] = zco[i]; 
	topwaterindex[trajno][topwatercounter] = a[i]; 
	
	
	//printf("%d %d %lf %lf %lf %d \n", a[i],b[i],xco[i],yco[i],zco[i],topwatercounter);

	topwatercounter++;
 
	
      }

      if (b[i] == 1 && zco[i] < 100) { 

        botwater[trajno][botwatercounter] = zco[i]; 
	botwaterindex[trajno][botwatercounter] = a[i]; 
	
	//printf("%d %d %lf %lf %lf  \n", a[i],b[i],xco[i],yco[i],zco[i]);
	botwatercounter++;
	
      }

    }
    
    
    for (i=0;i<botcounter;i++) { 

      /* Now array botbilayer has the coordinates for the bottom surface of the bilayer */
      //printf("%d %d %lf \n",i,botindex[i],botbilayer[i]);

    }
    
    
    for (i=0;i<topcounter;i++) { 
      
      /* Now array botbilayer has the coordinates for the bottom surface of the bilayer */
      //printf("%d %d %d %lf \n",i,trajno,topindex[trajno][i],topbilayer[trajno][i]);
      
    }

    for (i=0;i<topwatercounter;i++) { 
      // printf("%d %d %d %lf \n",i,trajno,topwaterindex[trajno][i],topwater[trajno][i]);
    }

    for (i=0;i<botfirstwatercounter;i++) { 

      //  printf("%d %d %d %lf \n",i,trajno,botfirstwaterindex[0][i],botfirstwater[0][i]);
     
    }

    /*Now we can compare all the trajectories to the first one with the original cooridinates for the surfaces of the bilayer*/

    /* First, we test if the top waters have gone down*/

 
    int watercounter; 

    watercounter = 0;

    for(i=0;i<topwatercounter;i++) { 

      // for (u=0;u<botcounter;u++) {

      //  printf("%lf \n",zco[botindex[0][i]]);
	
	/* Check each part in the trajectory with array[1][] */ 

	for (j=0;j<numberofatoms;j++) {
	  
	  if (a[j] == botfirstwaterindex[0][i] && zco[j] > zco[topfirstindex[0][i]]) { 

	    watercounter++; 

	    //printf("%d %d %d %lf \n", botfirstwaterindex[0][i],watercounter,a[j],zco[j]);

	  
	  }
	  
	}
      
	//}

    }
    
    printf("%lf %d %d \n",zco[71312],b[71312], watercounter);

    /* We now have data stored of an entire frame */
     /* Now we have to take into account the number of waters on each side of the bilayer */
     /* Each water is course-grained, so we have to measure the number of crossings of each water bead from  */ 
     /* Now decide if there is a pore forming, by counting the number of crossovers of water beads*/      

    
    //  for(trajno=0;trajno<numberoftraj;trajno++) { */
    
    //  printf("This is the data for trajectory no %d \n", trajno);  */
    
    /*     for (i=0;i<numberofatoms;i++) { */
    
  /*      printf("%d %d %lf %lf %lf \n", a[1],b[1],xco[1],yco[1],zco[1]); */
    
  /*      /\* Checking a value of coordinate[i] as we go down the trajectories*\/ */
    
  /*    } */
    
  /*  } */
  }

 fclose(ipf);

 exit(0);
}



      
