/*

Mobility Patterns Generator for ns-2 simulator
Based on a Community Based Mobility Model
	
Copyright (C) Mirco Musolesi University College London
              m.musolesi@cs.ucl.ac.uk

Version 0.3

ns-2 implementation

December 2006

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

The contact details of the author are the following:

Mirco Musolesi
Department of Computer Science - University College London
Gower Street London WC1E 6BT United Kingdom
Email: m.musolesi@cs.ucl.ac.uk
Phone: +44 20 7679 0391 Fax: +44 20 7387 1397
Web: http://www.cs.ucl.ac.uk/staff/m.musolesi


###
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <istream>
#include <iostream>
#include <cstdlib>
#include <ios>
#include "socialNet.h"
using namespace std;

#define RAND_INT(l,h) (((double)(random() * ((double)(h)-(l)+1.0))) + (l))

/* linear congruential generator.  Generator x[n+1] = a * x[n] mod m */

//default seed 
static unsigned  SEED = 93186752;

double getRandomDouble () {

/* input parameters*/

/* static int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901; */

static unsigned int a = 1223106847, m = 4294967291U, q = 3, r = 625646750;

/* static unsigned int a = 279470273, m = 4294967291U, q = 15, r = 102913196;*/
/* static unsigned int a = 1583458089, m = 2147483647, q = 1, r = 564025558; */
/* static unsigned int a = 784588716, m = 2147483647, q = 2, r = 578306215;  */
/* static unsigned int a = 16807, m = 2147483647, q = 127773, r = 2836;      */
/* static unsigned int a = 950706376, m = 2147483647, q = 2, r = 246070895;  */

   SEED = a*(SEED % q) - r*(SEED / q);
   return ((double)SEED / (double)m);
 }


// Sets the random seed of the random number generator
void setSeed (unsigned int init)   {
	if (init != 0.0) SEED = init;
}
		
		
// Returns a randomly generated double number in a specific range
double generateRandomDouble (double min, double max) {
	//double number=getRandomDouble() * (max-min+1.0) + min;
	double number=getRandomDouble()*(max-min)+min;	
	
	//double number=min+((max-min)*getRandomDouble()/RAND_MAX);
	return number;
}

// Returns a randomly generated integer number in a specific range
int generateRandomInteger(int min, int max) {
	
	int number=(int)floor(generateRandomDouble((double)min,(double)max)+1);
	//printf("%s %d\n","The random generated integer number is",number);

	return number;
}


// Print on stdout the header
void header() {
	printf("%s","\nMobility Patterns Generator for ns-2 simulator\n");
	printf("%s","Based on the Community Based Mobility Model\n");
	printf("%s","Version 0.3 December 2006\n");         
	printf("%s","Copyright (C) Mirco Musolesi University College London\n");
   	printf("%s","              m.musolesi@cs.ucl.ac.uk\n");
}


// Print on stdinput the usage
void usage() {
  	printf("%s","Usage: movGen [-options] \n");
	printf("%s","Input parameters:\n");
	printf("%s","-h                           shows help\n");
	printf("%s","-x                           generates XML\n");
	printf("%s","-n                           sets the number of hosts\n");
	printf("%s","-t <totalSimulationTime>     sets the total simulation time\n");
	printf("%s","-r <reconfigurationInterval> sets the reconfiguration interval\n");
	printf("%s","-s <lowerBoundSpeedHost>     sets lower bound of the speed of the hosts (in m/s)\n");
	printf("%s","-S <upperBoundSpeedHost>     sets upper bound of the speed of the hosts (in m/s)\n");
	printf("%s","-p <connectionThreshold>     sets the connection threshold\n");
	printf("%s","-X <sideLengthXcoordinates>  sets the side length of the the simulation area - x coordinates\n");
	printf("%s","-Y <sideLengthYcoordinates>  sets the side length of the the simulation area - y coordinates\n");
	printf("%s","-R <numberOfRows>            sets the number of rows\n");
	printf("%s","-C <numberOfColumns>         sets the number of columns\n");
	printf("%s","-T <transmissionRange>       sets the transmission range (in m)\n");
	printf("%s","-w <rewiring probability>    sets the rewiring probability\n");
	printf("%s","-G <numberOfGroups>          sets the initial number of groups for the Caveman model\n");
	printf("%s","-g <seedRNG>                 sets the seed of the Random Number Generator\n");
	printf("%s","-c <numberOfTravellers>      sets the number of travellers\n");
	printf("%s","-v <travellersSpeed>         sets the speed of the travellers\n");
	printf("%s","-a <on/off>                  sets the Girvan-Newman algorithm on/off\n");
	printf("%s","-d <on/off>                  sets the deterministic selection of the nodes on/off. If true the selection is deterministic, if false is probabilistic\n"); 
	printf("%s","-A <on/off>                  sets the colocation traces on/off\n");
	printf("%s","-b <on/off>                  sets the communities traces on/off\n");
}


// Generates XML header 
void generateXMLHeader(FILE *outputFileXML, double sideLength_WholeAreaX, double sideLength_WholeAreaY, double transmissionRange, int numberOfNodes) {
	
	fprintf(outputFileXML,"<?xml version=\"1.0\"?>\n");
	fprintf(outputFileXML,"<simulation xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
	fprintf(outputFileXML,"xsi:noNamespaceSchemaLocation=\"http://www.i-u.de/schools/hellbrueck/ansim/xml/simulation.xsd\">\n");
	fprintf(outputFileXML,"\t<parameter>\n");
	fprintf(outputFileXML,"\t\t<area_shape>rectangle</area_shape>\n");
	fprintf(outputFileXML,"\t\t<xsize>%f</xsize>\n",sideLength_WholeAreaX);
	fprintf(outputFileXML,"\t\t<ysize>%f</ysize>\n",sideLength_WholeAreaY);
	fprintf(outputFileXML,"\t\t<extension_factor>1</extension_factor>\n");
	fprintf(outputFileXML,"\t\t<numberOfNodes>%d</numberOfNodes>\n",numberOfNodes);
	fprintf(outputFileXML,"\t\t<range>%f</range>\n",transmissionRange);
	fprintf(outputFileXML,"\t\t<mobility_model>socialFoundedMM</mobility_model>\n");
	fprintf(outputFileXML,"\t</parameter>\n");
			
}


// Generates XML section related to the definition of the node position
void generateXMLSetNodePosition(FILE *outputFileXML, int node, double x, double y,double sideLength_WholeAreaX, double sideLength_WholeAreaY) {
	fprintf(outputFileXML,"\t\t<node>\n");
	fprintf(outputFileXML,"\t\t\t<node_id>%d</node_id>\n",node);
	fprintf(outputFileXML,"\t\t\t<position>\n");
	fprintf(outputFileXML,"\t\t\t\t<xpos>%f</xpos>\n",x-sideLength_WholeAreaX/2.0);
	fprintf(outputFileXML,"\t\t\t\t<ypos>%f</ypos>\n",y-sideLength_WholeAreaY/2.0);
	fprintf(outputFileXML,"\t\t\t</position>\n");
	fprintf(outputFileXML,"\t\t</node>\n");
}


// Generates XML section related to the definition of a new goal
void generateXMLSetNewGoal(FILE *outputFileXML, int node, double time, double newGoalX, double newGoalY, double velocity, double sideLength_WholeAreaX, double sideLength_WholeAreaY) {
	fprintf(outputFileXML,"\t\t<position_change>\n");
	fprintf(outputFileXML,"\t\t\t<node_id>%d</node_id>\n",node);
	fprintf(outputFileXML,"\t\t\t<start_time>%f</start_time>\n",time);
	fprintf(outputFileXML,"\t\t\t<destination>\n");
	fprintf(outputFileXML,"\t\t\t\t<xpos>%f</xpos>\n",newGoalX-sideLength_WholeAreaX/2.0);
	fprintf(outputFileXML,"\t\t\t\t<ypos>%f</ypos>\n",newGoalY-sideLength_WholeAreaY/2.0);
	fprintf(outputFileXML,"\t\t\t</destination>\n");
	fprintf(outputFileXML,"\t\t\t<velocity>%f</velocity>\n",velocity);
	fprintf(outputFileXML,"\t\t</position_change>\n");
}
        

int main (int argc, char *argv[]) {

	//----------------------------
	//Definitions of the constants
	//----------------------------

	//number of hosts
	int numHosts=100;

	//connection threshold
	double connectionThreshold=0.1;

	//length of the sides of the area
	double sideLength_WholeAreaX=2000.0;
	double sideLength_WholeAreaY=2000.0;

	//radius of the transmission area
	double radius=200.0;

	//number of rows
	int numberOfRows=20;

	//number of columns
	int numberOfColumns=20;
	
	//simulationTime
	double totalSimulationTime=1000.0;

	//refreshTime
	double stepInterval=1.0;

	//low bound speed of the host
	double minHostSpeed=1.0;
 	double maxHostSpeed=10.00;

	//reconfiguration interval
	//it defines the interval between two reconfigurations
	double reconfigurationInterval=10000.0;

	//rewiring probability
	double rewiringProb=0.1;
	
	//numberOfGroups
	int numberOfGroups=10;

	//simulation time
	double simTime=0.0;
	
	//number of travellers
	int numberOfTravellers=0;

	//speed of the travellers
	double travellerSpeed=20.0;
    
	//verbose on/off
	bool verbose=false;
	
	//generateXML on/off
	bool generateXML=false;
	
	//girvan Newman clustering algorithm on off
	bool girvanNewmanOn=true;
	
	//colocation traces on/off
	bool colocationTracesOn=false;
	
	//communities traces on/off
	bool communitiesTracesOn=false;
	
	//deterministic on/off
	bool deterministicOn=true;
	
	//drift
	float drift;
	
	//debug variables
	bool firstTime=true;
	double timeFirstReconfiguration;
	
							
	
  	//***********************
	// Retrieving input values
	// ***********************

	char ch;
	int inputSeed;	
		
	while ((ch = getopt (argc, argv, "hxy:n:s:S:p:l:d:X:Y:c:C:R:G:r:T:t:w:v:g:a:d:A:b:"))!=EOF)
		
		 switch(ch) {
		 	 case 'h':
			 			//shows helps
			 			header();
			 			usage();
			 			exit(1);
			 case 'x':
			 			//generate XML file using ANSim format
			 			generateXML=true;			
			 			break;
			 case 'n':
			 			//sets the number of hosts
			 			numHosts=atoi(optarg);
			 			break;
			 case 's':
			 			//sets lower bound of the speed of the hosts
			 			minHostSpeed=(double)atof(optarg);
			 			break;
			 case 'S':										 
						//sets upper bound of the speed of the hosts
						maxHostSpeed=(double)atof(optarg);
						break;
			 case 'p':
			 			//sets the connection threshold
			 			connectionThreshold=(double)atof(optarg);
			 			break;
			 case 'X':
			 			//sets the side length of the whole area (x coordinates)
			 			sideLength_WholeAreaX=(double)atof(optarg);
			 			break;
			 case 'Y':
						//sets the side length of the whole area (y coordinates)
						sideLength_WholeAreaY=(double)atof(optarg);
						break;
			 case 'R':
						//sets the number of rows
						numberOfRows=atoi(optarg);
						break;
			 case 'C':
						//sets the number of columns
					    numberOfColumns=atoi(optarg);
						break;	
			 case 'T':
			 			//sets the transmission range (radius)
			 			radius=(double)atof(optarg);			
			 			break;
			 case 't':
			 			//sets total simulation time
			 			totalSimulationTime=(double)atof(optarg);
			 			break;
			 case 'r':
						//sets the reconfiguration interval
						reconfigurationInterval=(double)atof(optarg);
						break;	
			 case 'w':	
						//sets the rewiring probability
						rewiringProb=(double)atof(optarg);
						break;
			 case 'G':   
						//sets the number of groups of the social network
						numberOfGroups=atoi(optarg);
						break;
			 case 'g':
			 			//sets the seed of the random number generator
			 			inputSeed=atoi(optarg);
			 			setSeed(inputSeed);
			 			break;
			 case 'c':	
						//sets the number of carriers/travellers
						numberOfTravellers=atoi(optarg);
						break;
			 case 'v':  
						//sets the speed of the travellers
						travellerSpeed=(double)atof(optarg);
						break;
			 case 'a': 
						//sets the Girvan/Newman algorithm on/off
						if (strcmp(optarg,"on")==0)
							girvanNewmanOn=true;
						else
						    girvanNewmanOn=false;
						break;	
			 case 'A':  
						//sets colocation traces on/off
						if (strcmp(optarg,"on")==0)
							colocationTracesOn=true;
						else
							colocationTracesOn=false;
						break;
			 case 'b': 
						//sets communities on/off
						if (strcmp(optarg,"on")==0)
							communitiesTracesOn=true;
						else
							communitiesTracesOn=false;
						break;
		     case 'd':
						//sets deterministic on/off
						if (strcmp(optarg,"on")==0)
							deterministicOn=true;
						else
							deterministicOn=false;
						break;
			 case ':':   
                        //fprintf(stderr,"Option -%c requires an operand\n", optopt);
						usage();
                        exit(1);			
             case '?':
					    //
          		     	usage();
          		    	exit(1);        
			 default: 
						break;
		}
	

		
	bool isConnected[numHosts][numHosts];

	struct hostsItem{
		double currentX;
		double currentY;
		double relativeX;
		double relativeY;
		double goalRelativeX;
		double goalRelativeY;
		double goalCurrentX;
		double goalCurrentY;
		int cellIdX;
		int cellIdY;
		double speed;
		double absSpeed;
		double af;
		bool isATraveller;
	};

	hostsItem hosts[numHosts];

	struct cellsItem{
		double minX;
		double minY;
		int numberOfHosts;
	};
	
	cellsItem cells[numberOfRows][numberOfColumns];

	double sideLengthX;
	double sideLengthY;

	//cell attractivity
	float CA[numberOfRows][numberOfColumns];

	//Variables used for the generation of the mobility scenario using the results of complex networks theory and in particular the Girvan-Newman algorithm

	//interaction->adjacency threshold
	double threshold=connectionThreshold;

	int numberOfMembers[numHosts];
	int **adjacency;
	double **interaction;
	int **groups;

	//output stream
	std::ofstream out1; 

	//input fileName
	char resultsFileName[]="results.txt";	
	
	//double last values registered
	double lastValues[numHosts][numHosts];
	
	//array [hostId,communitiesId]
	int communities[numHosts][2];
	
	//probability of moving to the cell [c][r]
	float a[numberOfRows][numberOfColumns];
	
	struct probRange {
		float min;
		float max;
	};
	
	probRange p[numberOfRows][numberOfColumns];
	
	//probability space distribution in the range [0..1]

	// **********************************	
	// Initialisation of the output files	
	// **********************************
	
	//Creation of the ns-2 traces output file
	//if ((outputFile = fopen ("socMov.tr","wt"))==NULL)
	//	fprintf (stderr,"Cannot open %s\n","socMov.tr");

	//Creation of the output file XML
	FILE* outputFileXML;
	generateXML=true;
	if (generateXML==true)
		if ((outputFileXML=fopen("socMov.xml","wt"))==NULL)
			fprintf(stderr,"Cannot open %s\n","socMov.xml");
	
	
	//Creation of the colocation statistics files
	FILE* outputColocFiles[numHosts];
	if (colocationTracesOn==true)
		for (int i=0;i<numHosts;i++) {
			static char fileName[10];
			static char prefix[]=".idat";
			sprintf(fileName,"%d%s",i,prefix);
            if ((outputColocFiles[i]=fopen(fileName,"wt"))==NULL)
				fprintf(stderr,"Cannot open %s\n",fileName);
		}//end for
	
	//Creation of the communities statistics files
	FILE* outputCommunitiesFile;
	if (communitiesTracesOn==true)
		if ((outputCommunitiesFile=fopen("communities.idat","wt"))==NULL)
			fprintf(stderr,"Cannot open communities.idat\n");
			
		
	// ***************************
	// Initialisation of the model
	//****************************
	sideLengthX=sideLength_WholeAreaX/((double)numberOfRows);
	sideLengthY=sideLength_WholeAreaY/((double)numberOfColumns);
		
	for (int i=0;i<numberOfRows;i++)
		for (int j=0;j<numberOfColumns;j++) {
			cells[i][j].minX=((double)i)*sideLengthX;
			cells[i][j].minY=((double)j)*sideLengthY;	
			cells[i][j].numberOfHosts=0;
		}	

	//setup of the links
	for (int i=0;i<numHosts;i++) {

		for (int l=0;l<numHosts;l++) {
			isConnected[i][l]=false;
			lastValues[i][l]=0.0;
			//nodeModule1->setIsInReach(l,false);
		}
	}
	
	//initialization of the travellers
	for (int i=0; i<numberOfTravellers; i++) {
		
		hosts[i].isATraveller=true;
	    //definition of the initial speeds of the travellers	
	    hosts[i].speed=travellerSpeed;
		
	}
	
	for (int i=numberOfTravellers; i<numHosts; i++) {
	
		hosts[i].isATraveller=false;
		hosts[i].speed=generateRandomDouble(minHostSpeed,maxHostSpeed);
	
	}

	double numberOfReconfigurations=0.0;
	double nextReconfigurationTime=0.0;
	
	
	int initialNumberOfGroups=numberOfGroups;
	//Update of the positions

	//adjacency=initialise_int_array(numHosts);
	//interaction=initialise_double_array(numHosts);
	//groups=initialise_int_array(numHosts);

	for (simTime=0.0; simTime<totalSimulationTime;simTime=simTime+stepInterval) {

	//reconfiguration mechanism
		if (simTime==nextReconfigurationTime) {
		
			for (int i=0;i<numberOfRows;i++)
				for (int j=0;j<numberOfColumns;j++) {
	
					cells[i][j].numberOfHosts=0;
				}	
			
			//start correction
			adjacency=initialise_int_array(numHosts);
			print_int_array(adjacency,numHosts);
         	interaction=initialise_double_array(numHosts);
			
			print_double_array(interaction,numHosts);
			
         	groups=initialise_int_array(numHosts);
			//end correction
			
			numberOfReconfigurations=numberOfReconfigurations+1.0;
			nextReconfigurationTime=reconfigurationInterval*numberOfReconfigurations;
			

 			bool splitted=true;
			double previousModth=0.0;
			double modth=0.1;
		
			if (girvanNewmanOn==true) {	
				
				initialise_weight_array_ingroups(interaction,numHosts,initialNumberOfGroups,rewiringProb,threshold,inputSeed);
				generate_adjacency(interaction,adjacency,threshold,numHosts);					
				
				//clustering using the Girvan-Newman algorithm
				do {
					for (int i=0;i<numHosts; i++)
						numberOfMembers[i]=0;
						splitted=false;
						double betw[numHosts];
						for (int i=0;i<numHosts; i++)
							betw[i]=0;
						calculate_betweenness(betw,adjacency,numHosts);
				
						for (int i=0;i<numHosts; i++)
							numberOfMembers[i]=0;
						numberOfGroups=getGroups(adjacency, groups ,numberOfMembers,numHosts);
				
						if (verbose==true)
							printGroups(numberOfGroups, groups, numberOfMembers, numHosts);	
						previousModth=modth;
			
						modth=splitNetwork_Threshold(adjacency, betw, numHosts, modth);
						if (verbose==true)
							cout<<"\nModth is equal to "<<modth;
				
				}
				while ((previousModth<modth)&&(modth>-1));
				print_double_array(interaction, numHosts);
			}//end if (girvanNewmanOn==true)
			else {
				//communities based on the initial number of caves in the Caveman model		
				//i.e., w=0
				initialise_weight_array_ingroups(interaction,numHosts,initialNumberOfGroups,0,threshold,inputSeed);
				generate_adjacency(interaction,adjacency,threshold,numHosts);

				for (int i=0;i<numHosts; i++)
					numberOfMembers[i]=0;
				numberOfGroups=getGroups(adjacency, groups ,numberOfMembers,numHosts);
				if (verbose==true)
					printGroups(numberOfGroups, groups, numberOfMembers, numHosts);
				
				initialise_weight_array_ingroups(interaction,numHosts,initialNumberOfGroups,rewiringProb,threshold,inputSeed);
				generate_adjacency(interaction,adjacency,threshold,numHosts);					
				
				print_double_array(interaction, numHosts);
				
			}//end else
		
			
			//for (int i=0;i<numberOfGroups;i++) {
			//	cout<<"\n";
			//	cout<<"The members of group "<<i+1<<" are: ";
			//	for (int j=0;j<numberOfMembers[i];j++)
			//		cout<<groups[i][j]<<" ";
			//	cout<<"\n";
				
			int pointer=0; 
			for (int i=0;i<numberOfGroups;i++) {
				
				int cellIdX=generateRandomInteger(0,numberOfRows);
				int cellIdY=generateRandomInteger(0,numberOfColumns);
				
				
				for (int j=0;j<numberOfMembers[i];j++) {
					int hostId=groups[i][j];
					hosts[hostId-1].cellIdX=cellIdX;
					hosts[hostId-1].cellIdY=cellIdY;
					
					communities[pointer][0]=hostId;
					communities[pointer][1]=i+1;
					pointer++;
					
					//increment the number of the hosts in that cell
					cells[cellIdX-1][cellIdY-1].numberOfHosts+=1;
				}
	 			
				
			}
			
			if (communitiesTracesOn==true) {
				
				int temp1=0;
				int temp2=0;
				for (int i=0;i<numHosts;i++) {
					for (int j=0;j<numHosts-1;j++) {
						
						if (communities[j][0]>communities[j+1][0]) {
							temp1=communities[j+1][0];
							temp2=communities[j+1][1];
							
							communities[j+1][0]=communities[j][0];
							communities[j+1][1]=communities[j][1];
							
							communities[j][0]=temp1;
							communities[j][1]=temp2;
						}
					}//end for (int j=0
				}// end for (int i=0
				
				
				for (int i=0;i<numHosts;i++) {
					fprintf(outputCommunitiesFile, "%d%s", communities[i][0]," ");
					fprintf(outputCommunitiesFile, "%d%s", communities[i][1],"\n");
			
				}

				fclose(outputCommunitiesFile);			
			
			}//end communitiesTracesOn
			
			
			if (firstTime==true) {
				//definition of the initial position of the hosts
				for (int k=0;k<numHosts;k++) {			
					hosts[k].currentX=cells[hosts[k].cellIdX-1][hosts[k].cellIdY-1].minX+generateRandomDouble(0.0,sideLengthX);	
					hosts[k].currentY=cells[hosts[k].cellIdX-1][hosts[k].cellIdY-1].minY+generateRandomDouble(0.0,sideLengthY);
				}
					
				timeFirstReconfiguration=simTime;
				firstTime=false;
	
				// ************************
				// Generation of the traces
				// ************************
	
				printf("%s","#Mobility Patterns Generator for ns-2 simulator\n");	
				printf("%s","#Based on the Community Based Mobility Model\n");
				printf("%s","#Copyright (C) Mirco Musolesi University College London\n");
				printf("%s","#              m.musolesi@cs.ucl.ac.uk\n");
				printf("%s","#Version 0.3 December 2006\n");                

				printf("%s","set god_ [God instance]\n");
	
				//Generating initial positions of the hosts
				for (int i=0;i<numHosts;i++) {
			 		
					printf("%s%d%s%f%s","$node_(",i,") set X_ ",hosts[i].currentX,"\n");	
	 		 		printf("%s%d%s%f%s","$node_(",i,") set Y_ ",hosts[i].currentY,"\n");	
			 		printf("%s%d%s%f%s","$node_(",i,") set Z_ ",0.0,"\n");
						  	  
				}
  	
				//Generating initial positions of the hosts - XML code 	
				if (generateXML==true) {
  		
					generateXMLHeader(outputFileXML, sideLength_WholeAreaX, sideLength_WholeAreaY, radius, numHosts);
					fprintf(outputFileXML,"\t<node_settings>\n");
  		
					for (int i=0;i<numHosts;i++) {
						generateXMLSetNodePosition(outputFileXML, i, hosts[i].currentX, hosts[i].currentY,sideLength_WholeAreaX,sideLength_WholeAreaY);
					}
  	
					fprintf(outputFileXML,"\t</node_settings>\n");
					fprintf(outputFileXML,"\t<mobility>\n");
  	
				}//end generation XML code
			}//end if (firstTime==true)
			
			//definition of the initial goals
			for (int k=0;k<numHosts;k++) {
				hosts[k].goalCurrentX=cells[hosts[k].cellIdX-1][hosts[k].cellIdY-1].minX+generateRandomDouble(0.0,sideLengthX);
				hosts[k].goalCurrentY=cells[hosts[k].cellIdX-1][hosts[k].cellIdY-1].minY+generateRandomDouble(0.0,sideLengthY);				
			}	
		
		
			//generation of the traces - setting of the goals
			for (int i=0;i<numHosts;i++) {
				
				hosts[i].absSpeed=hosts[i].speed;

				printf("%s %f %s%d%s %f %f %f%s","$ns_ at",simTime,"\"$node_(",i,") setdest",hosts[i].goalCurrentX,hosts[i].goalCurrentY,(hosts[i].absSpeed)/stepInterval,"\"\n");  					
   		 	
   			 if (generateXML==true) 
  		  			generateXMLSetNewGoal(outputFileXML,i,simTime,hosts[i].goalCurrentX,hosts[i].goalCurrentY,hosts[i].absSpeed/stepInterval,sideLength_WholeAreaX,sideLength_WholeAreaY);	
			} //end for (int i=0;i<numHosts;i++)
		}//end reconfiguration mechanism
	
		for (int i=0;i<numHosts;i++) {

				if((hosts[i].currentX>hosts[i].goalCurrentX+hosts[i].speed)||
   				   (hosts[i].currentX<hosts[i].goalCurrentX-hosts[i].speed)||
				   (hosts[i].currentY>hosts[i].goalCurrentY+hosts[i].speed)||
				   (hosts[i].currentY<hosts[i].goalCurrentY-hosts[i].speed)) {
						//move towards the goal					
						if (hosts[i].currentX<(hosts[i].goalCurrentX-hosts[i].speed))
							hosts[i].currentX=hosts[i].currentX+hosts[i].speed;
	  					if (hosts[i].currentX>(hosts[i].goalCurrentX+hosts[i].speed))
							hosts[i].currentX=(hosts[i].currentX)-hosts[i].speed;
						if (hosts[i].currentY<(hosts[i].goalCurrentY-hosts[i].speed))
							hosts[i].currentY=(hosts[i].currentY)+hosts[i].speed;
						if (hosts[i].currentY>(hosts[i].goalCurrentY+hosts[i].speed))
							hosts[i].currentY=(hosts[i].currentY)-hosts[i].speed; 
					}
				else {
					
					int selectedGoalCellX=0;
					int selectedGoalCellY=0;	
					int previousGoalCellX=hosts[i].cellIdX;
					int previousGoalCellY=hosts[i].cellIdY;
					
					if (deterministicOn==true) {
					
						//Algorithm of the selection of the new cell
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++)
								CA[c][r]=0;
							
						for (int n=0;n<numHosts;n++)
							if ((n<i)||(n>i))
								CA[hosts[n].cellIdX-1][hosts[n].cellIdY-1]+=interaction[i][n];	
					
						for (int c=0;c<numberOfRows;c++) 
							for (int r=0;r<numberOfColumns;r++) {
								
								if (cells[c][r].numberOfHosts!=0) 
									CA[c][r]=CA[c][r]/(double)cells[c][r].numberOfHosts;
									//CA[c][r]=CA[c][r];								
								else
									CA[c][r]=0;
								//printf ("Deterministic %f ",CA[c][r]);

								}	

						int selectedGoalCellX2=0;
						int selectedGoalCellY2=0;
						
						double CAMax1=0;
						double CAMax2=0;
					
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++)
								if (CA[c][r]>CAMax1) {
								
									//set the second best
									selectedGoalCellX2=selectedGoalCellX;
									selectedGoalCellY2=selectedGoalCellY;
									CAMax2=CAMax1;
								
									selectedGoalCellX=c+1;
									selectedGoalCellY=r+1;
								
									CAMax1=CA[c][r];

								}
									else if (CA[c][r]>CAMax2) {
										selectedGoalCellX2=c+1;
										selectedGoalCellY2=r+1;
										CAMax2=CA[c][r];
									}	
										
									if ((previousGoalCellX==selectedGoalCellX)&&(previousGoalCellY==selectedGoalCellY)) {
										if (hosts[i].isATraveller==true) {
											if (selectedGoalCellX!=0) 
											{
												selectedGoalCellX=selectedGoalCellX2;
												selectedGoalCellY=selectedGoalCellY2;
											}
										}
									}		
					}//end deterministic
					else {//probabilistic
						
						
						
						//Algorithm of the selection of the new cell
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++)
								CA[c][r]=0.0;
						
						for (int n=0;n<numHosts;n++)
							if ((n<i)||(n>i))
								CA[hosts[n].cellIdX-1][hosts[n].cellIdY-1]+=interaction[i][n];	
						
						for (int c=0;c<numberOfRows;c++) 
							for (int r=0;r<numberOfColumns;r++) {
								
								//if (cells[c][r].numberOfHosts!=0) 
									//CA[c][r]=CA[c][r]/(double)cells[c][r].numberOfHosts;
								//	CA[c][r]=CA[c][r];								
								//else
								//	CA[c][r]=0;
								//printf ("NonDeterministic %f ",CA[c][r]);
								
							}	
						
						
						//Denonmiantor for the normalization of the values
						float denNorm=0.00;
						
						//used for simulations
						//drift=0.01;
						
						drift=0.0;
						
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++) {
								denNorm=denNorm+CA[c][r]+drift;
								//printf ("%f ",CA[c][r]);
							}
						
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++)
								if (CA[c][r]==0)
									a[c][r]=drift/denNorm;
								else
									a[c][r]=(CA[c][r]+drift)/(+denNorm);
						
						float current=0.0;
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++) {
								p[c][r].min=current;
								
								p[c][r].max=p[c][r].min+a[c][r];
								current=p[c][r].max;
								
							}
						
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++) {		
								
								p[c][r].min=p[c][r].min;
								p[c][r].max=p[c][r].max;
								//printf("%s %d %d %s %f %f %s","Square",c,r," : ",p[c][r].min,p[c][r].max,"\n");
							}
							
						//float infiniteDice= (float) random() / (float) 0x7fffffff;;
								float infiniteDice = (float) generateRandomDouble(0.0,1.0);
						printf ("%f %s",infiniteDice,"\n");
						for (int c=0;c<numberOfRows;c++)
							for (int r=0;r<numberOfColumns;r++)
								if ((infiniteDice>p[c][r].min)&&(infiniteDice<p[c][r].max)) {
									
									selectedGoalCellX=c+1;
									selectedGoalCellY=r+1;
								
								}
									
					}
				
						
						
					//Re-definition of the number of hosts in each cell
					cells[previousGoalCellX-1][previousGoalCellY-1].numberOfHosts-=1;
					cells[selectedGoalCellX-1][selectedGoalCellY-1].numberOfHosts+=1;		
					
					double newGoalRelativeX=generateRandomDouble(0,sideLengthX);
					double newGoalRelativeY=generateRandomDouble(0,sideLengthY);						
					
					//refresh of the information
					hosts[i].cellIdX=selectedGoalCellX;
					hosts[i].cellIdY=selectedGoalCellY;
						
				    hosts[i].goalCurrentX=cells[selectedGoalCellX-1][selectedGoalCellY-1].minX+newGoalRelativeX;
				    hosts[i].goalCurrentY=cells[selectedGoalCellX-1][selectedGoalCellY-1].minY+newGoalRelativeY;								
					hosts[i].absSpeed=hosts[i].speed;
					
   		        	printf("%s %f %s%d%s %f %f %f%s","$ns_ at",simTime,"\"$node_(",i,") setdest",hosts[i].goalCurrentX,hosts[i].goalCurrentY,(hosts[i].absSpeed)/stepInterval,"\"\n");  					
   		 	
   			 if (generateXML==true) 
  		  		 	generateXMLSetNewGoal(outputFileXML,i,simTime,hosts[i].goalCurrentX,hosts[i].goalCurrentY,hosts[i].absSpeed/stepInterval,sideLength_WholeAreaX,sideLength_WholeAreaY);	
 
				}	
				
												
	}// end for (int i=0;i<numHosts;i++) {

	for (int i=0;i<numHosts;i++) {
		//update connectivity
		for (int j=0;j<numHosts;j++) {
				//I am connected with me so the follwing condition must hold...:)
				if (i!=j) {
					
					//calculation of the current distance
					double currentDistance=sqrt((hosts[i].currentX-hosts[j].currentX)*(hosts[i].currentX-hosts[j].currentX)+(hosts[i].currentY-hosts[j].currentY)*(hosts[i].currentY-hosts[j].currentY));
				
					//if currentDistance<=radius then the hosts are connected
					if (currentDistance<radius) {
						//if the hosts has been previously disconnected, then they must be connected
						if (isConnected[i][j]==false) {
							
							isConnected[i][j]=true;
							
							lastValues[i][j]=simTime;
						
						} 
					}
					//else they are disconnected
					else {						
						if (isConnected[i][j]==true) {
								if (simTime!=0) {
									//if the hosts has been previously connected, then they must be disconnected
									if (colocationTracesOn)
										fprintf(outputColocFiles[i],"%d %d %d\n",j+1,(long)lastValues[i][j],(long)simTime);
									isConnected[i][j]=false;
								}
							}
					}//end else 

				}//end if (i!=j)

			}

		}
				
	}
	//generation of the XML file - generation of the last lines
	if (generateXML==true)	{	
		fprintf(outputFileXML,"\t</mobility>\n");	
		fprintf(outputFileXML,"\t<statistics>\n");
		fprintf(outputFileXML,"\t\t<stoptime>%f</stoptime>\n",simTime);
		fprintf(outputFileXML,"\t</statistics>\n");
		fprintf(outputFileXML,"</simulation>\n");
	}
	
	// close output file and outputFileXML
	if (generateXML==true) 
		fclose (outputFileXML);

	
}

