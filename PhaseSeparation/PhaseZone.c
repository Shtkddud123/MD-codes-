#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>


double centreofmass(double H7, double H6_1,double H6_2,double T3_1, double T3_2, double T3_3,double T4);

int main () 
{ 


  /*Declaration of variables*/ 
  /*Declaring the paramters for the atoms in the box */
  /*Declare the number of atoms in each screenshot*/ 
  /*All the polymer atoms and the water atoms too*/
  
  const int numberofatoms = 17828; 
  int boxdim = 3; 
  double x,y,z; /*coordinates for the atoms in the box*/
  int atomtype; 

  /*x,y,z array to store all the values we read in from the dump file in order*/ 
  double xco[numberofatoms],yco[numberofatoms],zco[numberofatoms]; 
  /*ditto for a and b, which represent index and atomtype respectively*/ 
  int a[numberofatoms],b[numberofatoms];

  /*box dimensions*/
  double box1;
  double box2; 
  double boxlength[boxdim];
 

  /*Parameters to loop over an entire dump file*/ 
  int numberofSS = 384; /*The number of screenshots in the dump file*/
  int SSno = 0; /*the nth screenshot we are at*/
 

  /*misc*/
  char line[100];
  int n = 0; 
  int index;
  int i; 
  int j;
  int k = 0;
  int l = 0;

  /*Defining C12E2 atoms*/ 

  int A7; 
  int A6_1;
  int A6_2;
  int A3_1;
  int A3_2; 
  int A3_3; 
  int A4;

  /*Defining mimics*/ 

  int A13; 
  int A12_1;
  int A12_2;
  int A9_1;
  int A9_2; 
  int A9_3; 
  int A10;

  /*Array parameters to store in the COM of each molecule*/ 
  /*Array to store in C12E2 only*/ 

  int numberOfPolymers = 551; //The number of polymers of each type - C12E2 or mimic 

  double headGroupC12_E2xCOM[numberOfPolymers];
  double headGroupC12_E2yCOM[numberOfPolymers];
  double headGroupC12_E2zCOM[numberOfPolymers];
  

  /*Array to store in mimics only!*/

  
  double headGroupMIMICxCOM[numberOfPolymers];
  double headGroupMIMICyCOM[numberOfPolymers];
  double headGroupMIMICzCOM[numberOfPolymers];
 

  // int ww[numberofwaters];
  
  /*Start of code*/ 
  /*Reading in the dump file*/
  /*pointer for the dump file*/

  FILE *ipf; /* input file */
  
  /* open file for reading */
  ipf = fopen("dump.new", "r");
  
  /* check for error opening file */
    
  if(ipf == NULL) {
    
    printf("Error opening file\n");
    exit(1);
  }
  
  /* get a line from the file */
  /* fgets() returns NULL when it reaches an error or end of file */   
  
  int nlines = numberofatoms + 9;  
  int PolymerCounter=0;
  int MimicCounter=0; 


  for(SSno=0;SSno<numberofSS;SSno++) {

  
     //  printf("This is the data for trajectory no %d \n", SSno); 

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

	  //printf("%d %d %lf %lf %lf\n",index,atomtype,x,y,z);

	
      }

    }
    
    //This the counter to count the number of atoms completely surrounded by like atoms. 
    
    

    /*  for(i=0;i<3500;i++) {  */

    /*       printf("%i %i %lf %lf %lf %i\n",a[i],b[i],xco[i],yco[i],zco[i],i); */
    /*     } */
      //We can hold a value in the i array, then check it with the j loop
      //we want to put the r threshold (the zone where if in it, we want to consider it as phase separated)

      //We need to first catagorize polymer 7 6 3 4 as one type of polymer, and set 12 11 10 9 as the second. 
      //make sure that we do not read atoms of the same molecule 

      /*I STILL NEED TO ADD IN THE MIC*/ 


      //Defining the C12E2 and the mimics using the array; 


    /*There are 250 polymers, of C12E2 and mimic type*/
    


    PolymerCounter = 0;
    MimicCounter = 0;
    for (i=0;i<=249;i++) { 


      /*C12E2 atoms*/

      A7 =7*(i);
      A6_1=7*(i)+1;
      A6_2=7*(i)+2;
      A3_1 =7*(i)+3;
      A3_2 =7*(i)+4;
      A3_3 =7*(i)+5;
      A4 =7*(i)+6;

      /*mimic atoms*/

      A13 =7*(i) + 1751;
      A12_1=7*(i)+ 1752;
      A12_2=7*(i)+ 1753;
      A9_1 =7*(i)+ 1754;
      A9_2 =7*(i)+ 1755;
      A9_3 =7*(i)+ 1756;
      A10 =7*(i)+ 1757;

      headGroupC12_E2xCOM[i] = centreofmass(xco[A7],xco[A6_1],xco[A6_2],xco[A3_1],xco[A3_2],xco[A3_3],xco[A4]); 
      headGroupC12_E2yCOM[i] = centreofmass(yco[A7],yco[A6_1],yco[A6_2],yco[A3_1],yco[A3_2],yco[A3_3],yco[A4]); 
      headGroupC12_E2zCOM[i] = centreofmass(zco[A7],zco[A6_1],zco[A6_2],zco[A3_1],zco[A3_2],zco[A3_3],zco[A4]); 

      headGroupMIMICxCOM[i] = centreofmass(xco[A13],xco[A12_1],xco[A12_2],xco[A12_1],xco[A9_2],xco[A9_3],xco[A10]); 
      headGroupMIMICyCOM[i] = centreofmass(yco[A13],yco[A12_1],yco[A12_2],yco[A12_1],yco[A9_2],yco[A9_3],yco[A10]); 
      headGroupMIMICzCOM[i] = centreofmass(zco[A13],zco[A12_1],zco[A12_2],zco[A12_1],zco[A9_2],zco[A9_3],zco[A10]); 

      
     /*  printf("\n"); */
/*       printf("%i %lf %lf %lf  COM  %lf %lf %lf %i\n",b[A7],xco[A7],yco[A7],zco[A7],headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i],i); */
/*       printf("%i %lf %lf %lf  COM  %lf %lf %lf %i\n",b[A6_1],xco[A6_1],yco[A6_1],zco[A6_1],headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i],i); */
/*       printf("%i %lf %lf %lf  COM  %lf %lf %lf %i\n",b[A6_2],xco[A6_2],yco[A6_2],zco[A6_2],headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i],i); */
/*       printf("%i %lf %lf %lf  COM  %lf %lf %lf %i\n",b[A3_1],xco[A3_1],yco[A3_1],zco[A3_1],headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i],i); */
/*       printf("%i %lf %lf %lf  COM  %lf %lf %lf %i\n",b[A3_2],xco[A3_2],yco[A3_2],zco[A3_2],headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i],i); */
/*       printf("%i %lf %lf %lf  COM  %lf %lf %lf %i\n",b[A3_3],xco[A3_3],yco[A3_3],zco[A3_3],headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i],i); */
/*       printf("%i %lf %lf %lf  COM  %lf %lf %lf %i\n",b[A4],xco[A4],yco[A4],zco[A4],headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i],i); */



      //We shall ignore taking into account periodic boundary conditions for now// 

      //MAKE SURE TO ADD LATER!
     
      /*Now to check if the COM are within a distance r of each other*/ 

      /*For checking between C12E2 COMs*/ 

      for(j=0;j<=249;j++) {
	
	//set the distance between COM as 1 angstrom*/ 
	
 
      
	    if(headGroupC12_E2xCOM[i] !=headGroupC12_E2xCOM[j] && headGroupC12_E2yCOM[i] !=headGroupC12_E2yCOM[j] && headGroupC12_E2zCOM[i] !=headGroupC12_E2zCOM[j] ){ 

	      //if the atomtype is the same BUT the index of the atom is different 
	      
	      if(headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j] !=0 && sqrt(pow(headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j],2)) <= 10.000 ) {

		if(headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j] !=0 && sqrt(pow(headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j],2)) <= 10.000) {

		  if(headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j] !=0 && sqrt(pow(headGroupC12_E2zCOM[i] - headGroupC12_E2zCOM[j],2)) <= 10.000) {

		    // printf("%i %i %lf %lf %lf %i\n",b[i],b[j],headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j],headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j],headGroupC12_E2zCOM[i] - headGroupC12_E2zCOM[j],NewCounter);


		    //NewCounter takes into account the total number of signifcant interactions between like particles 

		    PolymerCounter++;
		  
		  }

		
		
	      }
	      
	    }
	    
	    
	    
	  }
	  
	    /*----------------------------------------------------------------------------------------------*/
	    
	    /*Now for checking between mimic COMs*/ 
	    
	    
	    if(headGroupMIMICxCOM[i] !=headGroupMIMICxCOM[j] && headGroupMIMICyCOM[i] !=headGroupMIMICyCOM[j] && headGroupMIMICzCOM[i] !=headGroupMIMICzCOM[j] ){ 
	      
	      //if the atomtype is the same BUT the index of the atom is different 
	      
	      if(headGroupMIMICxCOM[i] - headGroupMIMICxCOM[j] !=0 && sqrt(pow(headGroupMIMICxCOM[i] - headGroupMIMICxCOM[j],2)) <= 10.000 ) {

		if(headGroupMIMICyCOM[i] - headGroupMIMICyCOM[j] !=0 && sqrt(pow(headGroupMIMICyCOM[i] - headGroupMIMICyCOM[j],2)) <= 10.000) {
		  
		  if(headGroupMIMICyCOM[i] - headGroupMIMICyCOM[j] !=0 && sqrt(pow(headGroupMIMICzCOM[i] - headGroupMIMICzCOM[j],2)) <= 10.000) {

		    // printf("%i %i %lf %lf %lf %i\n",b[i],b[j],headGroupC12_E2xCOM[i] - headGroupC12_E2xCOM[j],headGroupC12_E2yCOM[i] - headGroupC12_E2yCOM[j],headGroupC12_E2zCOM[i] - headGroupC12_E2zCOM[j],NewCounter);
		    

		    //NewCounter takes into account the total number of signifcant interactions between like particles 
		    MimicCounter++;
		    
		  }
		  
		  
		
		}
		
	      }
	    
	      
	      
	    }
	    
	    
      }

      
    }
    
     printf("%i %.5f %.5f\n",SSno,PolymerCounter/250.00000,MimicCounter/250.00000);
    


   

   /*  while (i <= 249) { */


/*       A13 =7*(i) + 1751; */
/*       A12_1=7*(i)+ 1752; */
/*       A12_2=7*(i)+ 1753; */
/*       A9_1 =7*(i)+ 1754; */
/*       A9_2 =7*(i)+ 1755; */
/*       A9_3 =7*(i)+ 1756; */
/*       A10 =7*(i)+ 1757; */

/*       //  printf("%i %i %i %i %i %i %i %i \n",b[A13],b[A12_1],b[A12_2],b[A9_1],b[A9_2],b[A9_3],b[A10],i); */
/*       //  printf("%lf %lf %lf %lf %lf %lf %lf %i %i \n",xco[A13],xco[A12_1],xco[A12_2],xco[A9_1],xco[A9_2],xco[A9_3],xco[A10],i,j); */

/*       i++; */
     

/*     } */
  
  }
 return 0;

}
   





/*Function to calculate the centre of mass of each molecule*/

double centreofmass(double H7, double H6_1,double H6_2,double T3_1, double T3_2, double T3_3,double T4)

{

  double H7coord = H7;
  double H6coord_1 = H6_1;
  double H6coord_2 = H6_2;
  double T3coord_1 = T3_1; 
  double T3coord_2 = T3_2; 
  double T3coord_3 = T3_3; 
  double T4coord = T4; 
  double COM;

  COM = ((H7coord*31.03) + (H6coord_1*44.051400) + (H6coord_2*44.051400) + (T3coord_1*42.081000)  + (T3coord_2*42.081000)  +(T3coord_3*42.081000)  + (T4coord*29.062000))/(31.03 + (44.051400*2) + (42.081000*3) + 29.0620); 

  /*With this COM definition we now know the COM in each cartesian coordinate.*/ 

  return COM; 

  
}  


/*Function to take into account the periodic boundary conditions*/ 


/* void minimum_image(double *xcoord7, double *xcoord6_1, double *xcoord6_2, double *xcoord3_1, double *xcoordin3_2, double xcoordin3_3, double xcoordin4,double boxlength)  */

/* { */
  
/*   double box = xboxlength;  */

/*   //with 7, we move the headgroup  */

/*   if (*xcoord7 - *xcoordx6_1 > boxlength/2 ) { */

/*     *xccord7 = *xccord7 - xboxlength;  */
 
/*   } */

/*   else if (*xcoord7 - *xcoordx6_1 > boxlength/2) { */




/*   } */
  


/* } */


