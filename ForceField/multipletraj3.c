/* ------------------------------------------------------------------------------------------------------------------------------

This file processes a nanoparticle simulated inside a box of water and using the distance between the water molecules and the particle, creates a force field for use as tabulated values in simulation packages. 

This *is* work in progress, as the force field one gets from this code does not perfectly fit a neat function (such as the LJ potential and/or More potential). Hence, some optimization of the code is required. However, with time, it should provide a working code to provide tabulated values to use as a nanoparticle in MD simulations. 

--------------------------------------------------------------------------------------------------------------------------------------*/


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>

// Calculating of the Lennard-Jones potential using the sigma and epsilon parameters 
double LJ (double dist_ab, double eps_ab, double sigma_ab); 

// Function to calculate the distance between two points that has been defined by cartesian vectors  
double distance(double ax,double ay,double az, double bx, double by, double bz);

// Function to fix the minimum image convention due to overlapping regions between simulations. 

double mindist(double ax,double ay,double az, double bx, double by, double bz,
	       double bxx, double byy, double bzz);

// WIP 

void minimum_image(double *a,double *b,double *c,double bxx,double byy,double bzz);

// Calculating the Centre of Mass of each water molecule in the simulation 

double centreofmass(int Oa, int H1a, int H2a,int Ob, int H1b, int H2b,int Oc, int H1c, int H2c);

/* void continuity(double bx,double by,double bz,double bxx,double byy,double bzz); */

int main () 
{ 

  // A (rather crude..) image of what is being analysed. 

  /*             b (H2O)           */
  /*  a (Au)         c (H2O)       */
  /*             d (H2O)          */


  /*Declaration of variables*/ 
  /*Declaring the paramters for the atoms in the box */   
  const int numberofatoms = 3271; 
  int boxdim = 3; 
  const int nbins = 100;

  int numberofwaters = 1090; /*total number of water triads*/
  double x,y,z; /*coordinates for the atoms in the box*/

  double xco[numberofatoms],yco[numberofatoms],zco[numberofatoms]; 
  int a[numberofatoms],b[numberofatoms];
  
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

  /*Parameters to loop over an entire dump file*/ 
  int numberoftraj = 1000; /*The number of trajectories in the dump file*/
  int trajno = 0; /*the nth traj we are at*/ 

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
  int NewCounter;
  

  // int ww[numberofwaters];
  
  /*Start of code*/ 
  /*Reading in the data file*/  

  FILE *ipf; /* input file */
  
  /* open file for reading */
  printf("Opening file %s\n","Gold_simulation.dat");
  ipf = fopen("Gold_simulation.dat", "r");
  
  /* check for error opening file */
    
  if(ipf == NULL) {
    
    printf("Error opening file\n");
    exit(1);
  }
 
 
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

  /* Going through each trajectory snapshot - snapshot = trajno */

  for(trajno=0;trajno<numberoftraj;trajno++) {

     printf("\n");
     printf("This is the data for trajectory no %d \n", trajno); 

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

    }
   
    /*We now have data stored of an entire frame */

    // printf("%lf %lf %lf  \n", boxlength[0], boxlength[1], boxlength[2]);

    //  printf("%d %d %lf %lf %lf \n", a[3270],b[3270],xco[3270],yco[3270],zco[3270]);      

    /*We are now printing out the COM of each water molecule*/ 


    for (i=0;i<numberofwaters;i++) {

      firstoxygen = 3*(i);
      firsthydrogen =  3*(i) + 1;
      lasthydrogen = 3*(i) + 2;

      double Dx, Dy, Dz;
      
      Dx = xco[firsthydrogen] - xco[firstoxygen]; 
      Dy = yco[firsthydrogen] - yco[firstoxygen];
      Dz = zco[firsthydrogen] - zco[firstoxygen];
      
      minimum_image(&Dx,&Dy,&Dz,boxlength[0],boxlength[1],boxlength[2]); 

      xco[firsthydrogen] =  xco[firstoxygen] + Dx;
      yco[firsthydrogen] =  yco[firstoxygen] + Dy;
      zco[firsthydrogen] =  zco[firstoxygen] + Dz;
      
      Dx = xco[lasthydrogen] - xco[firstoxygen]; 
      Dy = yco[lasthydrogen] - yco[firstoxygen];
      Dz = zco[lasthydrogen] - zco[firstoxygen];

      minimum_image(&Dx,&Dy,&Dz,boxlength[0],boxlength[1],boxlength[2]); 

      xco[lasthydrogen] =  xco[firstoxygen] + Dx;
      yco[lasthydrogen] =  yco[firstoxygen] + Dy;
      zco[lasthydrogen] =  zco[firstoxygen] + Dz;
      
    }
    
    for (i=0;i<numberofwaters;i++) {
      
      firstoxygen = 3*(i);
      firsthydrogen =  3*(i) + 1;
      lasthydrogen = 3*(i) + 2;
      
      /*Now to calculate the COM of each water molecule */ 
      
      COMsetx  = (xco[firstoxygen]*16.0 + xco[firsthydrogen]*1.0 + xco[lasthydrogen]*1.0)/(16.0 + 1.0 + 1.0);
      COMsety  = (yco[firstoxygen]*16.0 + yco[firsthydrogen]*1.0 + yco[lasthydrogen]*1.0)/(16.0 + 1.0 + 1.0);
      COMsetz  = (zco[firstoxygen]*16.0 + zco[firsthydrogen]*1.0 + zco[lasthydrogen]*1.0)/(16.0 + 1.0 + 1.0);
      
      /*A random big number for the distance which will be truncated later in the code*/
      
      closestdistance1 = bignumber;
      closestdistance2 = bignumber; 
      
      for (j=0;j<numberofwaters;j++) { 
	
	if (i != j) { 
	  
	  newoxygen = 3*(j);
	  newfirsthydrogen =  3*(j) + 1;
	  newlasthydrogen =3*(j) + 2;

	  COMbx  = (xco[newoxygen]*16.0 + xco[newfirsthydrogen]*1.0 + xco[newlasthydrogen]*1.0)/(16.0 + 1.0 + 1.0);
	  COMby  = (yco[newoxygen]*16.0 + yco[newfirsthydrogen]*1.0 + yco[newlasthydrogen]*1.0)/(16.0 + 1.0 + 1.0);
	  COMbz  = (zco[newoxygen]*16.0 + zco[newfirsthydrogen]*1.0 + zco[newlasthydrogen]*1.0)/(16.0 + 1.0 + 1.0);


	  /*After applying the PBC to each triad, we work out the COMdist*/ 
	  

	  /*minimum image convention; Here we are describing the absolute distance between each COM of each water*/
 

	  COMdist = mindist(COMsetx,COMsety,COMsetz,COMbx,COMby,COMbz,boxlength[0],boxlength[1],boxlength[2]);
	  
	  // COMdist = distance(COMsetx,COMsety,COMsetz,COMbx,COMby,COMbz);

	  //  printf("%lf \n", COMdist); 

	  if (COMdist < closestdistance1) {

	    /*'If we find any distance which is not the minimum,reset that distance to be closestdistance2*/ 
	    
	    closestdistance2 = closestdistance1;

	    /*'Hence, set the water thought to be closest to be the second closest'*/ 

	    closestwater2 = closestwater1; 

	    /*'Reset the closestdistance1 to be COMdist to make it smaller*/
	    
	    closestdistance1 = COMdist;

	    /*The closestwater described characteristic of the array*/

	    closestwater1 = j; 
	    
	  }
	  
	  else if (COMdist < closestdistance2) { 
	    
	    closestdistance2 = COMdist;
	    closestwater2 = j; 
	  
	  }
	  
	}
     
      }
   
      w1[i] = i; 
      w2[i] = closestwater1; 
      w3[i] = closestwater2; 

    
   
      //printf("The waters closest to water %d are %d and %d,at distances = %8.3lf and %8.3lf \n\n",w1[i],w2[i],w3[i],closestdistance1,closestdistance2);
   
    }
  
    /*Calculating the weights of each water molecule*/ 

    for (i=0;i<numberofwaters;i++) {
      
      watercounter[i] = 0; 

    }
       
    for (i=0;i<numberofwaters;i++) { 

      watercounter[w1[i]]++; 
      watercounter[w2[i]]++; 
      watercounter[w3[i]]++; 
      

      // printf("%d %d %d \n",  watercounter[w1[i]], watercounter[w2[i]], watercounter[w3[i]]);
    }

     
      
    /* 	/\* 'How many times does w1[i] (the water at the c position) appear in all the water triads?'*\/ */
    /* 	/\*'How many times does w2[i] (the water at the the b position) appear in all the water triads?'*\/ */
     
    /* 	/\*'How many times does w3[i] (the water at the d position) appear in all the water triads?'*\/  */
 
    //printf("The overall weight of the triad of positions %d %d %d is %8.3lf \n", w1[i],w2[i],w3[i],weight);   
  
    /*Calculating the interaction*/ 
    /*Calculating the COM for each triad*/ 
     
      for (i=0;i<numberofwaters;i++) { 
	
	
	/*Checking the entire trajectory*/ 
      
	/*Defining the COM for each axes for the triad*/ 
       
	Obcoordinatex = (xco[3*w1[i]]);
	Obcoordinatey = (yco[3*w1[i]]); 
	Obcoordinatez = (zco[3*w1[i]]);

	Occoordinatex = (xco[3*w2[i]]);
	Occoordinatey = (yco[3*w2[i]]); 
	Occoordinatez = (zco[3*w2[i]]);

	Odcoordinatex = (xco[3*w3[i]]);
	Odcoordinatey = (yco[3*w3[i]]); 
	Odcoordinatez = (zco[3*w3[i]]);
    
      /*Defining the LJ parameters*/
	
	COMgoldCOMcrosssigma = (goldsigma + Osigma)/2;
	COMgoldCOMcrossep = sqrt(goldepsilon*Oepsilon);


	distgoldoxygenb = mindist(xco[3270],yco[3270],zco[3270],
				  Obcoordinatex,Obcoordinatey,Obcoordinatez,
				  boxlength[0],boxlength[1],boxlength[2]);
	distgoldoxygenc = mindist(xco[3270],yco[3270],zco[3270],
				  Occoordinatex,Occoordinatey,Occoordinatez,
				  boxlength[0],boxlength[1],boxlength[2]);
	distgoldoxygend = mindist(xco[3270],yco[3270],zco[3270],
				  Odcoordinatex,Odcoordinatey,Odcoordinatez,
				  boxlength[0],boxlength[1],boxlength[2]);

	/* distgoldoxygenb = distance(xco[3270],yco[3270],zco[3270], */
	/* 			  Obcoordinatex,Obcoordinatey,Obcoordinatez); */
	
	/* distgoldoxygenc = distance(xco[3270],yco[3270],zco[3270], */
	/* 			  Occoordinatex,Occoordinatey,Occoordinatez); */
				 
	/* distgoldoxygend = distance (xco[3270],yco[3270],zco[3270], */
	/* 			  Odcoordinatex,Odcoordinatey,Odcoordinatez); */
				
	/* COM of triad parameters*/ 
	
	/*3w#[i] represents the oxygen in each water molecule*/
	/*3*#w[i]+1 represents the first hydrogen in each water molecule.*/
	/*3*#w[i]+2 represents the last hydrogen in each water molecule.*/ 
	


	double COM1waterx,COM1watery,COM1waterz; 
	double COM2waterx,COM2watery,COM2waterz;
	double COM3waterx,COM3watery,COM3waterz;

	COM1waterx = (xco[3*w1[i]]*16 + xco[3*w1[i]+1]*1.0 + xco[3*w1[i]+2]*1.0)/18.0;
	COM1watery = (yco[3*w1[i]]*16 + yco[3*w1[i]+1]*1.0 + yco[3*w1[i]+2]*1.0)/18.0;
	COM1waterz = (zco[3*w1[i]]*16 + zco[3*w1[i]+1]*1.0 + zco[3*w1[i]+2]*1.0)/18.0;


	COM2waterx = (xco[3*w2[i]]*16 + xco[3*w2[i]+1]*1.0 + xco[3*w2[i]+2]*1.0)/18.0;
	COM2watery = (yco[3*w2[i]]*16 + yco[3*w2[i]+1]*1.0 + yco[3*w2[i]+2]*1.0)/18.0;
	COM2waterz = (zco[3*w2[i]]*16 + zco[3*w2[i]+1]*1.0 + zco[3*w2[i]+2]*1.0)/18.0;


	COM3waterx = (xco[3*w3[i]]*16 + xco[3*w3[i]+1]*1.0 + xco[3*w3[i]+2]*1.0)/18.0;
	COM3watery = (yco[3*w3[i]]*16 + yco[3*w3[i]+1]*1.0 + yco[3*w3[i]+2]*1.0)/18.0;
	COM3waterz = (zco[3*w3[i]]*16 + zco[3*w3[i]+1]*1.0 + zco[3*w3[i]+2]*1.0)/18.0;

	double DxCOM2,DyCOM2,DzCOM2; 
	double DxCOM3,DyCOM3,DzCOM3;
	
	DxCOM2 = COM2waterx - COM1waterx;
	DyCOM2 = COM2watery - COM1watery;
	DzCOM2 = COM2waterz - COM1waterz;

	minimum_image(&DxCOM2,&DyCOM2,&DzCOM2,boxlength[0],boxlength[1],boxlength[2]);

	DxCOM3 = COM3waterx - COM1waterx;
	DyCOM3 = COM3watery - COM1watery;
	DzCOM3 = COM3waterz - COM1waterz;	
  
	minimum_image(&DxCOM3,&DyCOM3,&DzCOM3,boxlength[0],boxlength[1],boxlength[2]);

	/*x component*/

	xmass = (3*COM1waterx + DxCOM2 + DxCOM3)/3;

	/*y component*/

	ymass =  (3*COM1watery + DyCOM2 + DyCOM3)/3;

	/*z component*/
	
	zmass =  (3*COM1waterz + DzCOM2 + DzCOM3)/3;

	/* xmass = (xco[3*w1[i]]*16.0 + xco[3*w1[i]+1]*1.0 + xco[3*w1[i]+2]*1.0 + xco[3*w2[i]]*16.0 + xco[3*w2[i]+1]*1.0 + xco[3*w2[i]+2]*1.0 + xco[3*w3[i]]*16.0 + xco[3*w3[i]+1]*1.0 + xco[3*w3[i]+2]*1.0)/(16.0*3 + 1.0*6);  */
	/* ymass = (yco[3*w1[i]]*16.0 + yco[3*w1[i]+1]*1.0 + yco[3*w1[i]+2]*1.0 + yco[3*w2[i]]*16.0 + yco[3*w2[i]+1]*1.0 + yco[3*w2[i]+2]*1.0 + yco[3*w3[i]]*16.0 + yco[3*w3[i]+1]*1.0 + yco[3*w3[i]+2]*1.0)/(16.0*3 + 1.0*6); */
	/* zmass = (zco[3*w1[i]]*16.0 + zco[3*w1[i]+1]*1.0 + zco[3*w1[i]+2]*1.0 + zco[3*w2[i]]*16.0 + zco[3*w2[i]+1]*1.0 + zco[3*w2[i]+2]*1.0 + zco[3*w3[i]]*16.0 + zco[3*w3[i]+1]*1.0 + zco[3*w3[i]+2]*1.0)/(16.0*3 + 1.0*6); */

	/* printf("water 1 %lf %lf %lf  \n",xco[3*w1[i]],yco[3*w1[i]+1],zco[3*w1[i]+2]); */
	/* printf("water 2 %lf %lf %lf  \n",xco[3*w2[i]],yco[3*w2[i]+1],zco[3*w2[i]+2]); */
	/* printf("water 3 %lf %lf %lf  \n",xco[3*w3[i]],yco[3*w3[i]+1],zco[3*w3[i]+2]); */
	/* printf("\n"); */
	
	/* Distance between each COM and the gold nanoparticle*/
	
	COMdistgoldwater = mindist(xco[3270],yco[3270],zco[3270],
				   xmass,ymass,zmass,
				   boxlength[0],boxlength[1],boxlength[2]);

	 /* COMdistgoldwater = distance(xco[3270],yco[3270],zco[3270], */
      	 /* 			  xmass,ymass,zmass); */
													   
	LJb = LJ(distgoldoxygenb, COMgoldCOMcrossep, COMgoldCOMcrosssigma);
	LJc = LJ(distgoldoxygenc, COMgoldCOMcrossep, COMgoldCOMcrosssigma);
	LJd = LJ(distgoldoxygend, COMgoldCOMcrossep, COMgoldCOMcrosssigma);       
	
	/*defining the weights*/  
	
	weight = (1/(double)watercounter[w1[i]] + 1/(double)watercounter[w2[i]] + 1/(double)watercounter[w3[i]])/3;
	
	//	printf("The overall weight of the triad of positions %d %d %d is %8.3lf \n", w1[i],w2[i],w3[i],weight);  
	
	/*LJ energy*/ 

	boxLJ = LJb + LJc + LJd; 
	
	
	//	printf("%lf %lf %lf \n",boxlength[1]*boxLJ,COMdistgoldwater,weight);

      
      
	
      
	for (j=0;j<nbins;j++) { 
	  
	  bin_min[j] = j*binwidth;
	  bin_max[j] = (j+1)*binwidth; 
     
	  if (COMdistgoldwater > bin_min[j] && COMdistgoldwater <= bin_max[j]) {  
	    
	    V_ave[j] +=boxLJ*weight;  
	    sumofweight[j] += weight;
	    bin_counter[j]++; 
	    newLJ[j] = V_ave[j]/sumofweight[j];
	    
	  }
	   
	   
	}
	 
	 
      }
        for (j=0;j<nbins;j++) {  
	
	  printf("%lf %lf %d \n", bin_min[j],newLJ[j],bin_counter[j]); 
		    
	}
	 
 
      
  }
      
   
   return(0); 
}
  
  

      
double LJ (double dist_ab, double eps_ab, double sigma_ab) { 
  
  double LJ;
  
  LJ = 4*eps_ab*(pow((sigma_ab/dist_ab),12)-pow((sigma_ab/dist_ab),6));
  
  if (LJ < 0.0 && dist_ab < 2.0) {

    printf("WARNING! \n");

    exit(0);

  }

  return LJ;

}



  
double distance(double ax,double ay,double az, double bx, double by, double bz) {

    double dist;

    dist = sqrt((bx-ax)*(bx-ax) + (by-ay)*(by-ay) + (bz-az)*(bz-az));

    return dist;
}

double mindist(double ax,double ay,double az, double bx, double by, 
	       double bz,double bxx, double byy, double bzz) {

    double dist;
    double dx; 
    double dy;
    double dz; 

    dx = bx - ax;
    dy = by - ay;
    dz = bz - az;
    
    
    minimum_image(&dx, &dy, &dz, bxx, byy, bzz);
    
    dist = sqrt(dx*dx + dy*dy + dz*dz);

    return dist;
}

  
void minimum_image(double *a,double *b,double *c,double bxx,double byy,double bzz)
{
  /* This function applies the minimum image convention to a vector */

  double half_bxx = bxx/2.0;
  double half_byy = byy/2.0;
  double half_bzz = bzz/2.0;
  double newa; 
  double newb;
  double newc;

 
  /*changing the distance between points to account for the minimum image criterion*/
     
  if (*a >= half_bxx) {

    *a = *a - bxx;
  }
  else if (*a < -half_bxx) {

      *a = *a + bxx;
  }

   if (*b >= half_byy) {

    *b = *b - byy;
  }
  else if (*b < -half_byy) {

      *b = *b + byy;
  }

  if (*c >= half_bzz) {

    *c = *c - bzz;
  }
  else if (*c < -half_bzz) {

      *c = *c + bzz;
  }
    

}
void continuity(double bx,double by,double bz,double bxx,double byy,double bzz)
{
  /* This function applies the minimum image convention to a vector */

  double half_bxx = bxx/2.0;
  double half_byy = byy/2.0;
  double half_bzz = bzz/2.0;
  double newa; 
  double newb;
  double newc;

 
  /*changing the distance between points to account for the minimum image criterion*/
     
  if (bx < -half_bxx) {
    
    bx = bx + bxx;
  }
  else if (bx >= half_bxx) {

      bx = bx - bxx;
  }   

if (by < -half_byy) {

    by = by + byy;
  }
  else if (by >= half_byy) {

      by = by - byy;
  }


if (bz < -half_bzz) {

    bz = bz + bzz;
  }

  else if (bz >= half_bzz) {

      bz = bz - bzz;
  }


}

double centreofmass(int Oa, int H1a,int H2a,int Ob, int H1b, int H2b,int Oc, int H1c, int H2c)

{

  double Oadum = Oa; 
  double Obdum = Ob;
  double Ocdum = Oc;
  double H1adum = H1a; 
  double H2adum = H2a;
  double H1bdum = H1b;
  double H2bdum = H2b; 
  double H1cdum = H1c;
  double H2cdum = H2c;
  double COM; 

  COM = (Oadum*16 + H1adum*1.0 + H2adum*1.0 + Obdum*16 + H1bdum*1.0 + H2bdum*1.0 + Ocdum*16 + H1cdum*1.0 + H2cdum*1.0)/(16*3 + 6); 

  /*With this COM definition we now know the COM in each cartesian coordinate.*/ 

  /*We now have to calculate the difference in distance between each COM point*/ 

 
  return COM; 
  
}  
