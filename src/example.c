/*
A simple example of using the gfx library.
CSE 20211
9/7/2011
by Prof. Thain
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gfx.h"

double Gamma = 0.80;
double IntensityMax = 255;

/** Taken from Earl F. Glynn's web page:
* <a href="http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm">Spectra Lab Report</a>
* */
void waveLengthToRGB(double Wavelength){
	double factor;
	double Red,Green,Blue;

	if((Wavelength >= 380) && (Wavelength<440)){
		Red = -(Wavelength - 440) / (440 - 380);
		Green = 0.0;
		Blue = 1.0;
	}else if((Wavelength >= 440) && (Wavelength<490)){
		Red = 0.0;
		Green = (Wavelength - 440) / (490 - 440);
		Blue = 1.0;
	}else if((Wavelength >= 490) && (Wavelength<510)){
		Red = 0.0;
		Green = 1.0;
		Blue = -(Wavelength - 510) / (510 - 490);
	}else if((Wavelength >= 510) && (Wavelength<580)){
		Red = (Wavelength - 510) / (580 - 510);
		Green = 1.0;
		Blue = 0.0;
	}else if((Wavelength >= 580) && (Wavelength<645)){
		Red = 1.0;
		Green = -(Wavelength - 645) / (645 - 580);
		Blue = 0.0;
	}else if((Wavelength >= 645) && (Wavelength<781)){
		Red = 1.0;
		Green = 0.0;
		Blue = 0.0;
	}else{
		Red = 0.0;
		Green = 0.0;
		Blue = 0.0;
	};

	// Let the intensity fall off near the vision limits

	if((Wavelength >= 380) && (Wavelength<420)){
		factor = 0.3 + 0.7*(Wavelength - 380) / (420 - 380);
	}else if((Wavelength >= 420) && (Wavelength<701)){
		factor = 1.0;
	}else if((Wavelength >= 701) && (Wavelength<781)){
		factor = 0.3 + 0.7*(780 - Wavelength) / (780 - 700);
	}else{
		factor = 0.0;
	};
	gfx_color(
		Red==0.0 ? 0 : (int) (IntensityMax * pow(Red * factor, Gamma)),
		Green==0.0 ? 0 : (int) (IntensityMax * pow(Green * factor, Gamma)),
		Blue==0.0 ? 0 : (int) (IntensityMax * pow(Blue * factor, Gamma)));

}
//380-781

void setcolor(double x)
{
	waveLengthToRGB(380+x*(781-380));
}

int main()
{
	int ysize = 300;
	int xsize = 300;
	int i,j;

	char c;

	// Open a new window for drawing.
	gfx_open(xsize,ysize,"Example Graphics Program");

	// Set the current drawing color to green.
	gfx_color(0,200,100);

	// Draw a triangle on the screen.
	gfx_line(100,100,200,100);
	gfx_line(200,100,150,150);
	gfx_line(150,150,100,100);
	for(i=0;i<xsize;i++)
	{
		setcolor(((double)i)/xsize);
		for(j=0;j<10;j++)
			gfx_point(i,j);
	}

	while(1) {
		// Wait for the user to press a character.
		c = gfx_wait();

		// Quit if it is the letter q.
		if(c=='q') break;
	}

	return 0;
}
