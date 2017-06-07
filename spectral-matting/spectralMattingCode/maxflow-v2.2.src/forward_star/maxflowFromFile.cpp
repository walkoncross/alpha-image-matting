// main function of maxflow. To be used with conjunction with the
// MATLAB function continuesVideo.m. For father information read the
// README file in the videoTextures directory, and the comments in
// findCut.m.
// This program uses the graph cut files made by the authors of [2]
// as described in the readme.pdf file.
//
// Writen By Yaniv Alon, yanival@cs.huji.ac.il.

#include <stdio.h>
#include <fstream>
#include <vector>
#include "graph.h"
using namespace std;

static const double c_infinity = 1e10;


int main()
{

    // define streams for files reading/writing
    ifstream dTime("graphTime.txt");
    ifstream dRight("graphRight.txt");
    ifstream dDown("graphDown.txt");
    ofstream outFile("minCutResults.txt");
    
    // read image size
    double imWidthT, imHeightT, imWidthR, imHeightR, imWidthD, imHeightD, imColor, framesNum;
    dTime>>imHeightT>>imWidthT>>imColor>>framesNum;
    dRight>>imHeightR>>imWidthR>>framesNum;
    dDown>>imHeightD>>imWidthD>>framesNum;
    int imSize=(int)(imWidthT*imHeightT);
    int imWidth=(int)(imWidthT);

    // for every pixel in the stream build node. Attach 6 edges to the node: Left, right, up, down, next frame, previous frame.
    // special cases: All pixels in first frame are connected to source. All pixels in last frame are connected to sink.
    // Pixels in image borders does not have edges going outside the image.
	vector<Graph::node_id> nodes((int)(imHeightT*imWidthT*framesNum));
	Graph *g = new Graph();
    int x,y,f,i;
    int firstPixelI=0;
    double dist;
    for (f=1; f<=framesNum; f++) {
      i=firstPixelI;
      // build nodes, and add edges to next frame if not in first frame
      if (f==1) // first frame 
        // add node for each pixel in the frame, and connect it to the source node with weight c_infinity
        for (y=1; y<=imHeightT; y++) 
          for (x=1; x<=imWidthT; x++, i++) { 
            nodes[i] = g -> add_node(); 
            g -> set_tweights(nodes[i], c_infinity, 0);
          }
      else if (f==framesNum) // last frame. connect to sink
        for (y=1; y<=imHeightT; y++) 
          for (x=1; x<=imWidthT; x++, i++) {
            nodes[i] = g -> add_node(); 
            g -> set_tweights(nodes[i], 0, c_infinity);
            dTime>>dist;
            g -> add_edge(nodes[i], nodes[i-imSize], dist, dist);
          }
      else   // not first or last frame
        // add node for each pixel in the frame, connect it to the source node with weight c_infinity, and add time edges
        for (y=1; y<=imHeightT; y++) 
          for (x=1; x<=imWidthT; x++, i++) {
            nodes[i] = g -> add_node(); 
            g -> set_tweights(nodes[i], 0, 0);
            dTime>>dist;
            g -> add_edge(nodes[i], nodes[i-imSize], dist, dist);
          }
      
      // add edges to the right
      i=firstPixelI;
      for (y=1; y<=imHeightR; y++) {
        for (x=1; x<=imWidthR; x++, i++) {  // note that miWidthR is usual 1 pixel smaller then imWidthT
          dRight>>dist;
          g -> add_edge(nodes[i], nodes[i+1], dist, dist);
        }
        i+=(int)(imWidthT-imWidthR); // compensate for the fact that imWidthR is smaller
      }
    
      // add edges to the bottom
      i=firstPixelI;
      for (y=1; y<=imHeightD; y++)
        for (x=1; x<=imWidthD; x++, i++) {  
          dDown>>dist;
          g -> add_edge(nodes[i], nodes[i+imWidth], dist, dist);
        }
      firstPixelI+=imSize;
    }

    // calc min cut
	Graph::flowtype flow = g -> maxflow();

    // write results into file
    i=0;
    for (f=1; f<=framesNum; f++) {
      for (y=1; y<=imHeightT; y++) {
        for (x=1; x<=imWidthT; x++, i++) 
          outFile<<(g->what_segment(nodes[i]) == Graph::SOURCE)<<" ";
        outFile<<endl;
      }
      outFile<<endl;
    }
        
	delete g;
	return 0;
}
