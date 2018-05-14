//=====================================================================================================
//Image and video dehazing based on the paper: An Investigation in Dehazing Compressed Images and Video
//Authors: Kristofor Gibson, D˜ung V˜o and Truong Nguyen
//Code by Yuze Chen
//chenyuze1988@gmail.com
//=====================================================================================================

#include <cv.h>
#include <highgui.h>

#include <opencv2/opencv.hpp>

#include <stdio.h>
#include <iostream>

using namespace cv;
using namespace std;


//median filtered dark channel
Mat getMedianDarkChannel(Mat src, int patch)
{
        Mat rgbmin = Mat::zeros(src.rows, src.cols, CV_8UC1);
        Mat MDCP;
        Vec3b intensity;

	int wsize, hsize;
	wsize = 5;
	hsize = 5;
	for(int m=hsize; m<src.rows-hsize; m++)
	{
		for(int n=wsize; n<src.cols-wsize; n++)
		{
			int minr = 0;
			int ming = 0;
			int minb = 0;
			intensity = src.at<Vec3b>(m,n);
			rgbmin.at<uchar>(m,n) = min(min(intensity.val[0], intensity.val[1]), intensity.val[2]);
		}
	}
	medianBlur(rgbmin, MDCP, patch);
	int i = 15;
	// bilateralFilter ( rgbmin, MDCP, i, i*2, i/2 );
	return MDCP;

}


//estimate airlight by the brightest pixel in dark channel (proposed by He et al.)
int estimateA(Mat DC)
{
	double minDC, maxDC;
	minMaxLoc(DC, &minDC, &maxDC);
	cout<<"estimated airlight is:"<<maxDC<<endl;
	//return 200;
	return maxDC;
}


//estimate transmission map
Mat estimateTransmission(Mat DCP, int ac)
{
	double w = 0.85;
	Mat transmission = Mat::zeros(DCP.rows, DCP.cols, CV_8UC1);
	Scalar intensity;

	for (int m=0; m<DCP.rows; m++)
	{
		for (int n=0; n<DCP.cols; n++)
		{
			intensity = DCP.at<uchar>(m,n);
			transmission.at<uchar>(m,n) = (1 - w * intensity.val[0] / ac) * 255;
		}
	}
	return transmission;
}


//dehazing foggy image
Mat getDehazed(Mat source, Mat t, int al)
{
	double tmin = 0.1;
	double tmax;
	
	Scalar inttran;
	Vec3b intsrc;
	Mat dehazed = Mat::zeros(source.rows, source.cols, CV_8UC3);

	for(int i=0; i<source.rows; i++)
	{
		for(int j=0; j<source.cols; j++)
		{
			inttran = t.at<uchar>(i,j);
			intsrc = source.at<Vec3b>(i,j);
			tmax = (inttran.val[0]/255) < tmin ? tmin : (inttran.val[0]/255);
			for(int k=0; k<3; k++)
			{
				dehazed.at<Vec3b>(i,j)[k] = abs((intsrc.val[k] - al) / tmax + al) > 255 ? 255 : abs((intsrc.val[k] - al) / tmax + al);
			}
		}
	}
	return dehazed;
}


#define USE_VIDEO
// #define USE_DEHAZE

int main(int argc, char** argv)
{
	//for video defogging
	
	VideoCapture vid;
	
	#ifdef USE_VIDEO
	if( argc == 1 )
		vid.open( 0 );
	else if( argc > 1 )
		vid.open( argv[1] );
		
	//VideoCapture vid(argv[1]);
	if(!vid.isOpened())
		return -1;
	double rate = vid.get(CV_CAP_PROP_FPS);
	int delay = 1000/rate;
        bool stop(false);
		#else
		
		double rate = 30.0;
			int delay = 1000/rate;
        bool stop(false);
		#endif
		
	Mat frame;
	Mat darkChannel;
	Mat T;
	Mat fogfree;
	double alpha = 0.05;    //alpha smoothing
	int Airlightp;          //airlight value of previous frame
	int Airlight;           //airlight value of current frame
	int FrameCount = 0;     //frame number
        int ad;                 //temp airlight value
        namedWindow("before and after", CV_WINDOW_AUTOSIZE);
	
	for(;;)
	{
	#ifdef USE_VIDEO
		vid >> frame;
	#else
		frame = imread( argv[1] );
		#endif
		
		//
		// invert frame
		//
		#ifndef USE_DEHAZE
		bitwise_not ( frame, frame );
		#endif
		
		FrameCount++;
		
		#ifdef USE_VIDEO
		if(vid.get(CV_CAP_PROP_POS_AVI_RATIO)==1)
			break;
		#endif
		
		//create mat for showing the frame before and after processing
	        Mat beforeafter = Mat::zeros(frame.rows, 2 * frame.cols, CV_8UC3);
	        Rect roil (0, 0, frame.cols, frame.rows);
	        Rect roir (frame.cols, 0, frame.cols, frame.rows);

		//first frame, without airlight smoothing
		if (FrameCount == 1)
		{
			darkChannel = getMedianDarkChannel(frame, 5);
		        Airlight = estimateA(darkChannel);
		        T = estimateTransmission(darkChannel, Airlight);
			ad = Airlight;
                        fogfree = getDehazed(frame, T, Airlight);
		}

		//other frames, with airlight smoothing
		else
		{
			double t = (double)getTickCount();

			Airlightp = ad;
			darkChannel = getMedianDarkChannel(frame, 9);
		        Airlight = estimateA(darkChannel);
		        // Airlight = 250	;
		        T = estimateTransmission(darkChannel, Airlight);

		        imshow( "transmission", T );
                        cout<<"previous:"<<Airlightp<<"--current:"<<Airlight<<endl;
		        ad = int(alpha * double(Airlight) + (1 - alpha) * double(Airlightp));//airlight smoothing
		        		        
		        int dmax = (255 - ad) * 6;
		        cv::normalize(T, T, dmax, 255 - dmax, NORM_MINMAX, CV_8UC1);
		        
		        cout<<"smoothed airlight is:"<<ad<<endl;
		        fogfree = getDehazed(frame, T, ad);

			t = (double)cvGetTickCount() - t;
			printf( "=============Execution time per frame = %gms\n", t/((double)cvGetTickFrequency()*1000.) );
		}
		#ifndef USE_DEHAZE
			        bitwise_not ( frame, frame );
			        			
	        bitwise_not ( fogfree, fogfree );
	        #endif
	        
	        frame.copyTo(beforeafter(roil));
	      	        
	        fogfree.copyTo(beforeafter(roir));
		imshow("before and after", beforeafter);

		if(waitKey(delay) >= 0)
			stop = true;
	}


	//for image defogging

	//Mat fog = imread("tiananmen1.bmp");
	//Mat darkChannel;
	//Mat T;
	//Mat fogfree;
	//Mat beforeafter = Mat::zeros(fog.rows, 2 * fog.cols, CV_8UC3);
	//Rect roil (0, 0, fog.cols, fog.rows);
	//Rect roir (fog.cols, 0, fog.cols, fog.rows);
	//int Airlight;
	//namedWindow("before and after", CV_WINDOW_AUTOSIZE);

	//darkChannel = getMedianDarkChannel(fog, 5);
	//Airlight = estimateA(darkChannel);
	//T = estimateTransmission(darkChannel, Airlight);
	//fogfree = getDehazed(fog, T, Airlight);

	//fog.copyTo(beforeafter(roil));
	//fogfree.copyTo(beforeafter(roir));
	//imwrite("./dehazed.jpg", fogfree);
	//imshow("before and after", beforeafter);
	//waitKey();

	return 0;
}
