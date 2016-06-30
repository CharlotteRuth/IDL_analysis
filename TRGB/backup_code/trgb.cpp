/*Charlotte Christensen
11/07/07
This program is to be called by TRGB.pro to determin the best fit to the luminosity function

g++ -I /astro/apps/boost/include/boost-1_33_1 -o trgb trgb.cpp

trgb Fake3 0
*/
#include <cstdlib> 
#include <boost/random.hpp>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

int binsearch(double target, double array[], int first, int last)
  //simple binary recursive search -- looks for closest value
{
  int index, mid;
  if(last - first > 1){
    mid = int((last + first)/2.0);
    if (array[mid] > target) index = binsearch(target, array, first, mid);
    else index = binsearch(target, array, mid, last);
  }
  else if(last - first == 1){
    if(fabs(array[first] - target) < fabs(array[last] - target)) index = first;
    else index = last;
  }
  else if(last == first){
    index = last;
  }
  else{
    index = -1;
  }
  return index;
}

double erf_m(double m, double m2, double sigma){
  //This is a simple error function
  return  0.398942/abs(sigma) * exp(-1.0*pow((m - m2),2.0)/2.0/sigma/sigma);
}

double broken_powerlaw(double m, double a, double b, double c, double d, double e,double mtrgb)
     //This function gives the functional form of the tip of the RGB
{
  double g, slope, y1, y2;
  m = m*-1.0;
  mtrgb = mtrgb*-1.0;
  if(m > (mtrgb + e/2)) g = pow(10.0,a*(m - mtrgb) + d);
  else if(m < (mtrgb - e/2)) g = pow(10.0,b*(m - mtrgb ) + c); 
  else 
    {
      y1 = a*(e/2.0) + d;
      y2 = b*(-1*e/2.0) + c;
      slope = (y1 - y2)/e;
      g = pow(10.0,(slope*(m - (mtrgb + e/2.0)) + y1));
    } 
  int temp;
  if (g < 0) {
    cout<<"G is below zero: "<<m<<" "<<mtrgb<<" "<<b<<" "<<c;
    cin>>temp;
  }
  return g;
}

double intg_broken_powerlaw(double min, double max, double a, double b, double c, double d, double e,double mtrgb)
  //This function returns the integral of the broken powerlaw
{
  min =   min*-1.0; 
  max =   max*-1.0;
  mtrgb = mtrgb*-1.0;
  double integral;/*, slope, y1, y2, nlog10;     
  y1 = a*(e/2.0) + d;
  y2 = b*(-1*e/2.0) + c;
  slope = (y1 - y2)/e;
  min =   min*-1.0; 
  max =   max*-1.0;
  mtrgb = mtrgb*-1.0;
  nlog10 = 2.302585093; //natural log of 10
 
  integral = (exp(nlog10*(a*(min - mtrgb) + d)) - exp(nlog10*(a*(e/2.0) + d)))/(nlog10*a)
           + (exp(nlog10*y1) - exp(nlog10*(-1.0*slope*e + y1)))/(nlog10*slope)
           + (exp(nlog10*(b*(-1.0*e/2) + c)) - exp(nlog10*(b*(max - mtrgb) + c)))/(nlog10*b);

  cout<<"y1: "<<y1<<", y2: "<<y2<<", slope: "<<slope<<endl;

  cout<<"AGB "<<max<<": "<<mtrgb-e/2.0<<" -> "<<(exp(nlog10*(b*(-1.0*e/2) + c)) - exp(nlog10*(b*(max - mtrgb) + c)))/(nlog10*b)<<endl;
  cout<<"Tip "<<mtrgb-e/2.0<<": "<<mtrgb+e/2<<" -> "<< (exp(nlog10*y1) - exp(nlog10*(-1.0*slope*e + y1)))/(nlog10*slope)<<endl;
  cout<<"RGB "<<mtrgb + e/2.0<<": "<<min<<" -> "<<(exp(nlog10*(a*(min - mtrgb) + d)) - exp(nlog10*(a*(e/2.0) + d)))/(nlog10*a)<<endl;*/

  //Let's try a different meathod and just divide by the greatest value in the distribuiton
  cout<<"Function at min: "<<pow(10.0,(a*(min - mtrgb) + d))<<"\nFunction at mtrgb + e/2: "<<pow(10.0,(a*(e/2.0) + d))<<endl;
  if (a >= 0) integral = pow(10.0,(a*(min - mtrgb) + d));
  else integral = pow(10.0,(a*(e/2.0) + d));

  return integral;
}

void findgen(double array[], int length, double min, double step){
     //fills an array with increasing elements
  for(int i = 0; i < length; i++)  array[i] = i*step + min;
}

double total_array(double array[],int start,int length){
     //returns the total of all elements in the array
  double total = 0;
  for(int i = start; i < length+start; i++) total = total + array[i];
  return total;
}

int findmax(double array[], int length){
  //returns the index of the maximum element (first in the array if there are several)
  int index = 0;
  double maxvalue = array[0];
  for(int i = 0; i < length; i++){
    if (array[i] != 0) cout<<" "<<array[i];
    if((array[i] > maxvalue && isnan(array[i]) == 0 && (array[i] != 0)) || maxvalue == 0){
      if (array[i] !=0 ) cout<<"*";
      maxvalue = array[i];
      index = i;
    }
  }
  cout<<endl;
  return index;
}

int findmin(double array[], int length){
  //returns the index of the maximum element (first in the array if there are several)
  int index = 0;
  double minvalue = array[0];
  for(int i = 0; i < length; i++){
    if((array[i] < minvalue && isnan(array[i]) == 0 && (array[i] != 0)) || minvalue == 0){
      minvalue = array[i];
      index = i;
    }
  }
  cout<<endl;
  return index;
}

double det_trgb(int numi, int numm, int lowerpos, int upperpos, int lowerpos_b, int upperpos_b, double imag[], double ierr[], double m[], double merr[], double lumfunc[], double parameters[])
     //This function does the work of finding the best fit -- To hell with Mendez.  I'm doing it my way
{
  int temp, maxindex, maxindex_b, maxindex_c;
  double a = parameters[0], d = parameters[3], e = parameters[4], bmin = -0.20, bmax = 1.5, db = 0.10, cmin = 1.8, cmax = 3.2, dc = 0.10,dm = m[2] - m[1], littlelum, error;
  //for M81, cmax = 3.2, cmin > 1.8;bmin = -0.20, bmax = 1.5,
  //for NGC, cmas = 5.2, cmin > 1,8, bmin = -0.5, bmax = 1.5
  int numb = int((bmax - bmin)/db)+1, numc = int((cmax - cmin)/dc)+1;
  double b[numb], b_results[numb], c[numc], c_results[numc], phi[numi], phi2[numm], mtrgb_results[numm];

  //*******************************Initializing Vectors and Zeroing out Parameters********************************
  findgen(b, numb, bmin, db);
  findgen(c, numc, cmin, dc);
  for(int mtrgbct = 0; mtrgbct <= numm; mtrgbct++) mtrgb_results[mtrgbct] = 0;
  for(int bct = 0; bct <= numb; bct++) b_results[bct] = 0;
  for(int cct = 0; cct <= numc; cct++) c_results[cct] = 0;
  if(lowerpos_b == -1) lowerpos_b = 0;
  if(upperpos_b == -1) upperpos_b = int(0.9*numm);

  //*******************************Making a luminosity function************************************************
  cout<<"Finding Maximum Likelyhood for Mag of TRGB\nLower limit: "<<m[lowerpos]<<", Upper limit: "<<m[upperpos]<<endl<<endl;
  for(int mct = 0; mct < numm; mct++) {    //Makes a luminosity function with an error spread
    for(int sct = 0; sct < numi ; sct++) phi[sct] = erf_m(imag[sct], m[mct], ierr[sct]);
    lumfunc[mct] = total_array(phi,0,numi);
    if (lumfunc[mct] == 0)lumfunc[mct] = pow(10.0,-3.0);
  }

  //********************************Marginalizing over Parameters to find the tip*********************************
  for(int mtrgbct = lowerpos; mtrgbct <= upperpos; mtrgbct++){ //marginalizes over b and c to find a good value for mtrgb
    for(int cct = 0; cct < numc; cct++){
      for(int bct = 0; bct < numb; bct++){
	for(int mct = lowerpos_b; mct < upperpos_b; mct++){ //Only fit luminosity function over this range of magnitudes
	  if(lumfunc[mct] > 1) error = sqrt(lumfunc[mct]);
	  else error = 1.0;
	  phi2[mct] =  log(erf_m(lumfunc[mct],broken_powerlaw(m[mct],a,b[bct],c[cct],d,e,m[mtrgbct]),error)); //Poisson Error (SQRT(N))
	  if (isnan(phi2[mct])  != 0)  cout<<mct<<" Nan: "<<" "<<m[mct]<<" "<<broken_powerlaw(m[mct],a, b[bct],c[cct],d,e,m[mtrgbct])<<" "<<lumfunc[mct]<<endl;	
	  if (isinf(phi2[mct]) != 0 || isnan(phi2[mct]) != 0) phi2[mct] = phi2[mct-1]; //If phi2 = log of  zero, set to -500
	}
	b_results[bct] = total_array(phi2,lowerpos_b,upperpos_b - lowerpos_b)*dm;
	if (isnan(b_results[bct]) == 1) cout<<"224, Nan: "<<cct<<" "<<total_array(phi2,lowerpos_b,upperpos_b - lowerpos_b+1)<<endl;
      }
      c_results[cct] = total_array(b_results,0,numb)*db;
    }
    mtrgb_results[mtrgbct] = total_array(c_results,0,numc)*dc;
  }
  maxindex = findmax(mtrgb_results,numm);
  cout<<endl<<"Integration Complete, Magnitude of TRGB: "<<m[maxindex]<<" ("<<maxindex-lowerpos<<")"<<endl<<endl;

  //*****************************Marginalizing to find parameters
  if (parameters[1] == -1.0) { //If the best parameters have not already been found
    cout<<"Finding Parameters"<<endl;
    for(int cct = 0; cct < numc; cct++){ //Marginalize over c to find best b
      for(int bct = 0; bct < numb; bct++){
	for(int mct = lowerpos_b; mct < upperpos_b; mct++){
	  if(lumfunc[mct] != 0) error = sqrt(lumfunc[mct]);
	  else error = 1.0;
	  phi2[mct] =  log(erf_m(lumfunc[mct],broken_powerlaw(m[mct],a, b[bct],c[cct],d,e,m[maxindex]),error));
	  if (isnan(phi2[mct]) == 1) cout<<mct<<" Nan: "<<" "<<m[mct]<<" "<<broken_powerlaw(m[mct],a, b[bct],c[cct],d,e,m[maxindex])<<" "<<lumfunc[mct]<<endl;
	  if (isinf(phi2[mct]) != 0) phi2[mct] = -500.0;			
	}
	b_results[bct] = total_array(phi2,lowerpos_b,upperpos_b - lowerpos_b)*dm;
	if (isnan(b_results[bct]) == 1) cout<<"Nan: "<<bct<<" "<<total_array(phi2,lowerpos_b,upperpos_b - lowerpos_b)*dm<<endl;
      }
      c_results[cct] = total_array(b_results,0,numb)*db;
    }
    maxindex_c = findmax(c_results,numc);
    parameters[2] = c[maxindex_c];
    cout<<"Best C: "<<parameters[2]<<endl;

    for(int bct = 0; bct < numb; bct++){
      for(int mct = lowerpos_b; mct < upperpos_b; mct++){
	  if(lumfunc[mct] != 0) error = sqrt(lumfunc[mct]);
	  else error = 1.0;
	phi2[mct] =  log(erf_m(lumfunc[mct],broken_powerlaw(m[mct],a,  b[bct],c[maxindex_c],d,e,m[maxindex]),error));
	if (isnan(phi2[mct]) == 1) cout<<mct<<" Nan: "<<" "<<m[mct]<<" "<<broken_powerlaw(m[mct],a,b[bct],c[maxindex_c] ,d,e,m[maxindex])<<" "<<lumfunc[mct]<<endl;
	if (isinf(phi2[mct]) != 0) phi2[mct] = -500.0;		
      }
      b_results[bct] = total_array(phi2,lowerpos_b,upperpos_b - lowerpos_b)*dm;
      if (isnan(b_results[bct]) == 1) cout<<"Nan: "<<bct<<" "<<total_array(phi2,lowerpos_b,upperpos_b - lowerpos_b)*dm<<endl;
    }
    maxindex_b = findmax(b_results,numb);
    parameters[1] = b[maxindex_b];    
    cout<<"Best B: "<<parameters[1]<<endl<<endl;
  }
  return m[maxindex];
}


//trgb.ymir
int main(int argc, char *argv[])
{
  //Main program.  
  //Reads in the data, finds the overall best fit and then finds the error by multiple calls to that function using a subset of the points
  int numi = 0, numm = 0, nruns = 50, flag = 0, datalength = 60000, lowerpos, upperpos, lowerpos_b, upperpos_b, numi_small, pos;
  double mean = 0, std = 1.0, prob, ranmag, minmag_small, maxmag_small, rangemag_small, univalue, univalue2, univalue1, normalvalue, normalizer;
  double tips[nruns + 1], error_tips[nruns], imag[datalength], ierr[datalength], binval[1000], lumfunc[1000], sigma[1000], parameters[5];
  char round[2], datafile_name[35], magfile_name[35], tipout_name[35], lumfout_name[35];
  boost::lagged_fibonacci19937 engine;
  boost::uniform_01 <boost::lagged_fibonacci19937> uni_dist(engine);
  boost::normal_distribution<double> norm_dist(mean,std);  
  ifstream datafile;//file containing the stars (mag and err)
  ifstream magfile;//file containing the range of mags and corresponding sigmas
  ofstream tipout;
  ofstream lumfout;

  //****************************Creating File Names********************************************
  strcat(datafile_name,"../Datafiles/trgbstars");
  strcat(datafile_name,argv[1]);
  strcat(datafile_name,"_");
  strcat(datafile_name,argv[2]);
  strcat(datafile_name,".dat");

  strcat(magfile_name,"../Datafiles/magANDerr");
  strcat(magfile_name,argv[1]);
  strcat(magfile_name,"_");
  strcat(magfile_name,argv[2]);
  strcat(magfile_name,".dat");

  strcat(tipout_name,"../Datafiles/found_tips");
  strcat(tipout_name,argv[1]);
  strcat(tipout_name,"_");
  strcat(tipout_name,argv[2]);
  strcat(tipout_name,".dat");

  strcat(lumfout_name,"../Datafiles/lum_func");
  strcat(lumfout_name,argv[1]);
  strcat(lumfout_name,"_");
  strcat(lumfout_name,argv[2]);
  strcat(lumfout_name,".dat");

  //*****************************Reading in Information***********************************
  cout<<"Opening "<<datafile_name<<endl;
  datafile.open(datafile_name);
  if (!datafile.is_open()){
    cout<<"Data file not read";
    exit(1);
  }  
  while (! datafile.eof()&& numi < datalength ){
    datafile >> imag[numi];
    datafile >> ierr[numi];
    numi++;
  }
  numi = numi - 1;
  datafile.close();

  cout<<"Opening "<<magfile_name<<endl;;  
  magfile.open(magfile_name);
  if (!magfile.is_open()){
      cout<<"Magnitude and sigma file not read";
      exit(1);
  }
  while (! magfile.eof() ){
    if (flag == 0){
      flag = 1;
      magfile >> parameters[0] >> parameters[3]>>parameters[4];
    }
    else if (flag == 1){
      flag = 2;
      magfile >> lowerpos >> upperpos;
    }
    else if (flag == 2){
      flag = 3;
      magfile >> lowerpos_b >> upperpos_b>>numi_small;
    }      
    else {
      magfile >> binval[numm]>> sigma[numm];
      numm++;
    }
  }
  numm = numm - 1;
  magfile.close();  
  binval[0] = binval[1];
  sigma[0] = sigma[1];

  //******************************Finding Tip and Luminosity Function********************************************************
  cout<<"Read in files, numm: "<<numm<<", numi: "<<numi<<", lowerpos: "<<lowerpos<<", upperpos: "<<upperpos<<endl;
  parameters[1] = -1.0;
  tips[0] = det_trgb(numi, numm, lowerpos, upperpos, lowerpos_b, upperpos_b,imag, ierr, binval, sigma, lumfunc, parameters);

  lumfout.open(lumfout_name); //Write luminosity function
  for(int i = 0; i < numm; i++)
    {
      lumfout << binval[i] << " "<<lumfunc[i] << endl;
    } 
  lumfout<<parameters[0]<<" "<<parameters[1]<<endl<<parameters[2]<<" "<<parameters[3]<<endl<<parameters[4]<<"    0"<<endl;
  lumfout.close();


  //*******************************Finds the tip using all the stars******************************
  //  cout<<"Numi: "<<numi;
  double imag_mc[numi], ierr_mc[numi];
  //  minmag = binval[0]; maxmag = binval[numm-1]; rangemag = maxmag - minmag; ranmag; //Initialize ranges of files and so on
  minmag_small = binval[lowerpos_b]; maxmag_small = binval[upperpos_b]; rangemag_small = maxmag_small - minmag_small;
//*************MC check for error***********************************************
  /*for(int i = 0; i < nruns; i++)
    {
      for(int i2 = 0; i2 < numi; i2++)
     	{
  	  normalvalue = norm_dist.operator()(engine);
	  univalue = uni_dist.operator()();
  	  pos = int(univalue*numi);
  	  ierr_mc[i2] = ierr[pos];
  	  imag_mc[i2] = imag[pos] + normalvalue*ierr_mc[i2]; 
  	}
      cout<<"Trial "<<i<<" of "<<nruns<<endl;
      tips[i+1] = det_trgb(numi, numm, lowerpos, upperpos,lowerpos_b, upperpos_b, imag_mc, ierr_mc, binval, sigma, lumfunc, parameters);
      }
*/
  /*  ofstream lumtest;
  char lumtestfile[35];
  strcat(lumtestfile,"../Datafiles/lumtest_");
  strcat(lumtestfile,argv[1]);
  strcat(lumtestfile,"_");
  strcat(lumtestfile,argv[2]);
  strcat(lumtestfile,".dat");
  lumtest.open(lumtestfile);

  for(int i = 0; i < numi; i++)
    {
      lumtest << imag_mc[i] <<endl;
    }
    lumtest.close();*/

  //************New check of errors*****************************************
  normalizer = intg_broken_powerlaw(minmag_small, maxmag_small, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], tips[0]); 
  cout<<"Min magnitude: "<<minmag_small<<", Max mag: "<<maxmag_small<<", Range: "<<rangemag_small<<" numi_small: "<<numi_small<<endl;
  cout<<"Parameters A: "<<parameters[0]<<", B: "<<parameters[1]<<", C: "<<parameters[2]<<", D: "<<parameters[3]<<", E: "<<parameters[4]<<", MTRGB: "<<tips[0]<<"; Normalization: "<<normalizer<<endl;
  for(int i = 0; i < nruns; i++) 
    { 
      int i2 = 0;
      while (i2 < numi_small)
  	{
	  normalvalue = norm_dist.operator()(engine);
	  univalue1 = uni_dist.operator()();
	  univalue2 = uni_dist.operator()();
  	  ranmag = univalue1*rangemag_small + minmag_small;
	  prob = broken_powerlaw(ranmag, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], tips[0])/normalizer;
	  if(univalue2 < prob)
	    {
	      pos = binsearch(ranmag,binval,0,numm);
	      ierr_mc[i2] = sigma[pos];
	      imag_mc[i2] = ranmag + normalvalue*ierr_mc[i2]; 
	      i2++;
	    }
  	}
      tips[i+1] = det_trgb(numi_small, numm, lowerpos, upperpos,lowerpos_b, upperpos_b, imag_mc, ierr_mc, binval, sigma, lumfunc, parameters);
      }
  /*  ofstream tipout_test;
  char tipout_name_test[35];
  strcat(tipout_name_test,"../Datafiles/found_tipsTest_");
  strcat(tipout_name_test,argv[1]);
  strcat(tipout_name_test,"_");
  strcat(tipout_name_test,argv[2]);
  strcat(tipout_name_test,".dat"); 
  tipout_test.open(tipout_name_test);
  for(int i = 0; i < nruns + 1; i++)
    {
      tipout_test << error_tips[i] <<endl;
    }
    tipout_test.close();*/


  //******Testing the luminosity function generated
  ofstream lumtest;
  char lumtestfile[35];
  strcat(lumtestfile,"../Datafiles/lumtest_");
  strcat(lumtestfile,argv[1]);
  strcat(lumtestfile,"_");
  strcat(lumtestfile,argv[2]);
  strcat(lumtestfile,".dat");
  lumtest.open(lumtestfile);
  for(int i = 0; i < numm; i++)
    {
      lumtest <<binval[i] <<" "<< lumfunc[i]<<endl;
    }
    lumtest.close();


//***************Print results to a file***********************
  tipout.open(tipout_name);
  for(int i = 0; i < nruns + 1; i++)
    {
      tipout << tips[i] <<endl;
    }
  tipout.close();
  return 0;
}

