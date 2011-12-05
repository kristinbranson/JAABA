/* 
Copyright (C) 2010-11 Steve Branson

This code may not be redistributed without the consent of the authors.
*/

#include "svgPlotter.h"

SVGPlotter::SVGPlotter() {
  width = 600;
  height = 400;
  cssFile = title = xLabel = yLabel = NULL;
  xMin = yMin = HUGE_VAL;
  xMax = yMax = -HUGE_VAL;
  sx = sy = 1;
  plots = NULL;
  numPlots = 0;
  stretchable = true;

  pad = 10;
  padAxis = 80;
  padTitle = 50;
  axisLabelGap = 50;
  ticLength = 5;
  maxPoints = 5000;
}

SVGPlotter::~SVGPlotter() {
  if(cssFile) free(cssFile);
  if(title) free(title);
  if(xLabel) free(xLabel);
  if(yLabel) free(yLabel);

  if(plots) {
    for(int i = 0; i < numPlots; i++) {
      free(plots[i].x);  free(plots[i].y); 
      if(plots[i].lineClass) free(plots[i].lineClass);
      if(plots[i].pointClass) free(plots[i].pointClass);
      if(plots[i].legendName) free(plots[i].legendName);
    }
    free(plots);
  }
}


void SVGPlotter::AddPlot(double *x, double *y, int numPts, double xOffset, double xScale, const char *legendName, const char *lineClass, const char *pointClass) {
  int stride = 1;
  while(numPts/stride > maxPoints) 
    stride++;
  plots = (SVGPlot*)realloc(plots, sizeof(SVGPlot)*(numPlots+1));
  plots[numPlots].lineClass = lineClass ? StringCopy(lineClass) : NULL;
  plots[numPlots].pointClass = pointClass ? StringCopy(pointClass) : NULL;
  plots[numPlots].legendName = legendName ? StringCopy(legendName) : NULL;
  plots[numPlots].numPts = numPts/stride;
  plots[numPlots].x = (double*)malloc(sizeof(double)*plots[numPlots].numPts);
  plots[numPlots].y = (double*)malloc(sizeof(double)*plots[numPlots].numPts);
  for(int i = 0, j = 0; i < numPts; i += stride, j++) {
    plots[numPlots].x[j] = x ? x[i] : xOffset+i*xScale;   
    plots[numPlots].y[j] = y[i]; 
    if(plots[numPlots].x[j] < xMin) xMin = plots[numPlots].x[j];
    if(plots[numPlots].x[j] > xMax) xMax = plots[numPlots].x[j];
    if(plots[numPlots].y[j] < yMin) yMin = plots[numPlots].y[j];
    if(plots[numPlots].y[j] > yMax) yMax = plots[numPlots].y[j];
  }
  numPlots++;
}

bool SVGPlotter::DrawPlots(FILE *fout) {
  for(int i = 0; i < numPlots; i++) {
    if(plots[i].lineClass) {
      if(plots[i].legendName) fprintf(fout, "  <!-- %s -->", plots[i].legendName);
      fprintf(fout, "  <polyline fill=\"none\" class=\"%s\" points=\"", plots[i].lineClass);
      int num = 0;
      for(int j = 0; j < plots[i].numPts; j++) {
	if(plots[i].x[j] >= xMin && plots[i].x[j] <= xMax && plots[i].y[j] >= yMin && plots[i].y[j] <= yMax) {
	  fprintf(fout, "%s%f,%f", num ? " " : "", (float)((plots[i].x[j]-xMin)*sx), (float)(height-(plots[i].y[j]-yMin)*sy));
	  num++;
	}
      }
      fprintf(fout, "\" />\n");
    }
    if(plots[i].pointClass) {
      if(plots[i].legendName) fprintf(fout, "  <!-- %s -->", plots[i].legendName);
      int num = 0;
      for(int j = 0; j < plots[i].numPts; j++) {
	if(plots[i].x[j] >= xMin && plots[i].x[j] <= xMax && plots[i].y[j] >= yMin && plots[i].y[j] <= yMax) {
	  fprintf(fout, "  <circle class=\"%s\" cx=%f cy=%f /> ", plots[i].pointClass, 
		  (float)((plots[i].x[j]-xMin)*sx), (float)(height-(plots[i].y[j]-yMin)*sy));
	  num++;
	}
      }
      fprintf(fout, "\" />\n");
    }
  }
  return true;
}

bool SVGPlotter::DrawAxis(FILE *fout) {
  char text[1000], format[1000];
  fprintf(fout, "\n  <rect class=\"axis\" x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" />\n", width, height);

  int xd = (int)(LOG_B((xMax-xMin), 10));
  float xTicSpacing = pow(10.0f, xd);
  if((xMax-xMin) / xTicSpacing < 5) { xd--; xTicSpacing /= 2; }
  float x = ceil(xMin/xTicSpacing)*xTicSpacing;
  while(x <= xMax) {
    if(xd >= 0) sprintf(text, "%d", (int)x);
    else { sprintf(format, "%%.%df", -xd); sprintf(text, format, x); }
    fprintf(fout, "  <line class=\"tic\" x1=\"%f\" y1=\"%d\" x2=\"%f\" y2=\"%d\" />", (x-xMin)*sx, height, (x-xMin)*sx, height+ticLength);
    fprintf(fout, "  <text class=\"xTic\" x=\"%f\" y=\"%d\" text-anchor=\"middle\" >%s</text>", (x-xMin)*sx, height+ticLength+1, text);
    x += xTicSpacing;
  }
  fprintf(fout, "\n");

  int yd = (int)(LOG_B((yMax-yMin), 10));
  float yTicSpacing = pow(10.0f, yd);
  if((yMax-yMin) / yTicSpacing < 5) { yd--; yTicSpacing /= 2; }
  float y = ceil(yMin/yTicSpacing)*yTicSpacing;
  while(y <= yMax) {
    if(yd >= 0) sprintf(text, "%d", (int)y);
    else { sprintf(format, "%%.%df", -yd); sprintf(text, format, y); }
    fprintf(fout, "  <line class=\"tic\" x1=\"%d\" y1=\"%f\" x2=\"%d\" y2=\"%f\" />", 0, height-(y-yMin)*sy, -ticLength, height-(y-yMin)*sy);
    fprintf(fout, "  <text class=\"yTic\" x=\"%d\" y=\"%f\" text-anchor=\"end\" >%s</text>", -ticLength-1, height-(y-yMin)*sy, text);
    y += yTicSpacing;
  }
  fprintf(fout, "\n");

  if(title) fprintf(fout, "  <text class=\"title\" text-anchor=\"middle\" x=\"%d\" y=\"%d\" >%s</text>\n", width/2, -4, title);
  if(xLabel) fprintf(fout, "  <text class=\"xLabel\" text-anchor=\"middle\" x=\"%d\" y=\"%d\" >%s</text>\n", width/2, height+ticLength+1+axisLabelGap, xLabel);
  if(yLabel) fprintf(fout, "  <text class=\"yLabel\" transform=\"rotate(-90)\" x=\"%d\" y=\"%d\" >%s</text>\n", -height/2, -ticLength-1-axisLabelGap, yLabel);

  return true;
}
bool SVGPlotter::DrawLegend(FILE *fout) {
  fprintf(fout, "\n<g id=\"legend\">\n  <rect id=\"legend_rect\" class=\"legend\" />\n");
  for(int i = 0; i < numPlots; i++)
    fprintf(fout, "  <polyline id=\"legend_line%d\" class=\"%s\" points=\"%d,%d %d,%d\"/> <text id=\"legend_text%d\" class=\"legend\" text-anchor=\"start\" x=\"%d\" y=\"%d\" >%s</text>\n", 
            i, plots[i].lineClass, 5, i*20+5, 25, i*20+5, i, 30, i*20+5, plots[i].legendName);
  fprintf(fout, "</g>\n");
  return true;
}

bool SVGPlotter::Save(const char *fname) {
  sx = width/(xMax-xMin);
  sy = height/(yMax-yMin);
  bool retval = false;
  const char *ext = GetFileExtension(fname);
  FILE *fout = fopen(fname, "w");
  if(!fout) { fprintf(stderr, "Failed to open '%s' for writing plot\n", fname); return false; }

  if(!strcmp(ext, "svg") || !strcmp(ext, "SVG")) {
    retval = SaveHeader(fout) &&
      SaveStyles(fout) &&
      DrawPlots(fout) &&
      DrawAxis(fout) &&
      DrawLegend(fout) &&
      SaveFooter(fout);
  } else if(!strcmp(ext, "m") || !strcmp(ext, "M")) {
    retval = SaveHeaderMatlab(fout) &&
      DrawPlotsMatlab(fout) &&
      DrawLegendMatlab(fout);
  } else {
    fprintf(stderr, "Error: plotter only supports file format .m or .svg\n"); 
    retval = false;
  }
  fclose(fout);

  return retval;
}


bool SVGPlotter::SaveHeader(FILE *fout) {
  return fprintf(fout, "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"") && 
    (stretchable || fprintf(fout, "width=\"%dpx\" height=\"%dpx\"", width, height)) &&
    fprintf(fout, " version=\"1.1\" viewBox=\"%f %f %f %f\" preserveAspectRatio=\"none\"", 
	    (float)(-padAxis), (float)(-padTitle), (float)(width+pad+padAxis), (float)(height+padTitle+padAxis)) &&
	    //xMin*sx-padAxis, yMin*sy-padTitle, xMax*sx+pad - (xMin*sx-padAxis), yMax*sy+padAxis - (yMin*sy-padTitle)) &&
    fprintf(fout, " onload=\"init()\" onmousemove=\"mouseMove(evt)\" >\n");
}
bool SVGPlotter::SaveFooter(FILE *fout) {
  return fprintf(fout, "</svg>\n");
}
bool SVGPlotter::SaveStyles(FILE *fout) {
  char line[1000];
  FILE *fin = fopen(cssFile ? cssFile : "plot.css", "r");
  assert(fin);
  while(fgets(line, 999, fin))
    fprintf(fout, "%s", line);
  fclose(fin);
  return true;
}


bool SVGPlotter::SaveHeaderMatlab(FILE *fout) {
  fprintf(fout, "figure(1); clf; hold on;\n");
  if(title) fprintf(fout, "title('%s');\n", title);
  if(xLabel) fprintf(fout, "xlabel('%s');\n", xLabel);
  if(yLabel) fprintf(fout, "ylabel('%s');\n", yLabel);
  fprintf(fout, "plotStyle = { '-', '--', '-.', '-v', '-p', '-s', '-x', '*', '+' };\n");
  fprintf(fout, "plotColor = { 'r', 'b', 'g', 'm', 'c', 'y', 'k' };\n");
  fprintf(fout, "labels = { };\n");
  return true;
}
bool SVGPlotter::DrawPlotsMatlab(FILE *fout) {
  for(int i = 0; i < numPlots; i++) {
    if(plots[i].legendName) fprintf(fout, "\n\nlabels{%d}='%s';\n", i+1, plots[i].legendName);
    fprintf(fout, "X%d=[", i);
    for(int j = 0; j < plots[i].numPts; j++) {
      fprintf(fout, "%s%f,%f", j ? ";" : "", (float)(plots[i].x[j]), (float)(plots[i].y[j]));
    }
    fprintf(fout, "];\n");
    fprintf(fout, "plot(X%d(:,1), X%d(:,2), plotStyle{%d},'LineWidth',3,'Color',plotColor{%d});\n", i, i, i+1, i+1);
  }
  return true;
}
bool SVGPlotter::DrawLegendMatlab(FILE *fout) {
  fprintf(fout, "legend(labels,'Location','Southeast');\n");
  return true;
}
