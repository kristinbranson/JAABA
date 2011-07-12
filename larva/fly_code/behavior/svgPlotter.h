/* 
Copyright (C) 2010-11 Steve Branson

This code may not be redistributed without the consent of the authors.
*/

#include <stdio.h>
#include "../blob.h"

typedef struct {
  char *lineClass;
  char *pointClass;
  double *x, *y; 
  int numPts;
  char *legendName;
} SVGPlot;

class SVGPlotter {
  char *cssFile;
  char *title;
  char *xLabel, *yLabel;
  double xMin, yMin, xMax, yMax;
  double sx, sy;
  int width, height;
  SVGPlot *plots;
  int numPlots; 
  bool stretchable;
  int pad, padAxis, padTitle;
  int axisLabelGap;
  int ticLength;

public:
  SVGPlotter();
  ~SVGPlotter();

  void AddPlot(double *x, double *y, int numPts, double xOffset=0, double xScale=1, const char *legend_name=NULL, const char *lineClass="default", const char *pointClass=NULL);
  bool Save(const char *fname);

  void SetCSS(const char *css) { cssFile = css ? StringCopy(css) : NULL; }
  void SetSize(int w, int h) { width = w; height = h; }
  void SetXLabel(const char *l) { xLabel = StringCopy(l); }
  void SetYLabel(const char *l) { yLabel = StringCopy(l); }
  void SetTitle(const char *l) { title = StringCopy(l); }

  void SetXMin(double x) { xMin = x; }
  void SetXMax(double x) { xMax = x; }
  void SetYMin(double y) { yMin = y; }
  void SetYMax(double y) { yMax = y; }

private:
  bool SaveStyles(FILE *fout);
  bool SaveHeader(FILE *fout);
  bool SaveFooter(FILE *fout);
  bool DrawPlots(FILE *fout);
  bool DrawAxis(FILE *fout);
  bool DrawLegend(FILE *fout);
  
  bool SaveHeaderMatlab(FILE *fout);
  bool DrawPlotsMatlab(FILE *fout);
  bool DrawLegendMatlab(FILE *fout);
};
