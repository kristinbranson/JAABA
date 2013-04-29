/* 
Copyright (C) 2010-11 Steve Branson

This code may not be redistributed without the consent of the authors.
*/

#include "util.h"

/**
 * @file svgPlotter.h
 * @brief Simple, light-weight curve plotting tool that produces matlab plots or SVG plots that can be embedded in a web page.
 */

/**
 * @struct SVGPlot
 * @brief Stores data for a single curve on a plot
 */
typedef struct {
  char *lineClass;  /**< Class name defining the color and style used to draw this curve */
  char *pointClass;  /**< Class name defining the color and style used to draw points on this curve */
  double *x;  /**< An array of size numPts defining the x-coordinates of each point on this curve */
  double *y;  /**< An array of size numPts defining the y-coordinates of each point on this curve */
  int numPts; /**< the number of points in this curve */
  char *legendName;  /**< The name of this curve, as will be displayed in the legend for this plot */
} SVGPlot;

/**
 * @class SVGPlotter
 * @brief Simple, light-weight curve plotting tool that produces SVG plots that can be embedded in a web page.  The plots
 * can be generated in either matlab format or svg format
 */
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
  int maxPoints;

public:
  /**
   * @brief Create a new plotter
   */
  SVGPlotter();

  ~SVGPlotter();

  /**
   * @brief Add a curve to this plot
   * @param x An array of size numPts defining the x-coordinates of each point on this curve 
   * @param y An array of size numPts defining the y-coordinates of each point on this curve 
   * @param numPts The number of points in this curve
   * @param xOffset If x is null, x[i] is defined implicitly as xOffset+i*xScale
   * @param xScale If x is null, x[i] is defined implicitly as xOffset+i*xScale
   * @param legend_name The name of this curve, as will be displayed in the legend for this plot
   * @param lineClass Class name defining the color and style used to draw this curve
   * @param pointClass Class name defining the color and style used to draw points on this curve
   */
  void AddPlot(double *x, double *y, int numPts, double xOffset=0, double xScale=1, const char *legend_name=NULL, const char *lineClass="default", const char *pointClass=NULL);


  /**
   * @brief Save a plot to a .svg or .m file, as determined by the file extension
   * @param fname the filename in which to save the plot to
   */
  bool Save(const char *fname);

  /**
   * @brief Set the filename of the css file defining the appearance classes for curves and points for svg plotting
   * @param css the filename of the css file
   */
  void SetCSS(const char *css) { cssFile = css ? StringCopy(css) : NULL; }

  /**
   * @brief Set the size in pixels of the svg plot
   * @param w the width of the plot
   * @param h the height of the plot
   */
  void SetSize(int w, int h) { width = w; height = h; }

  /**
   * @brief Set the label of the x-axis
   * @param l the label of the x-axis
   */
  void SetXLabel(const char *l) { xLabel = StringCopy(l); }

  /**
   * @brief Set the label of the y-axis
   * @param l the label of the y-axis
   */
  void SetYLabel(const char *l) { yLabel = StringCopy(l); }

  /**
   * @brief Set the title of the plot
   * @param l the title of the plot
   */
  void SetTitle(const char *l) { title = StringCopy(l); }

  /**
   * @brief Manually set the starting value of the x-axis
   * @param x the starting value of the x-axis
   */
  void SetXMin(double x) { xMin = x; }

  /**
   * @brief Manually set the ending value of the x-axis
   * @param x the ending value of the x-axis
   */
  void SetXMax(double x) { xMax = x; }

  /**
   * @brief Manually set the starting value of the y-axis
   * @param y the starting value of the y-axis
   */
  void SetYMin(double y) { yMin = y; }

  /**
   * @brief Manually set the ending value of the y-axis
   * @param y the ending value of the y-axis
   */
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
