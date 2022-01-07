#ifndef CanvasTools_hh
#define CanvasTools_hh

#include <TCanvas.h>
#include <TPad.h>

void PartitionCanvas(TCanvas* canvas, int xDivisions, int yDivisions, float leftMargin, float rightMargin, float topMargin, float bottomMargin, float verticalSpacing = 0.01, float horizontalSpacing = 0.01) {

  if (!canvas)
    return;

  float verticalStep  = (1.0 - bottomMargin - topMargin - (yDivisions - 1) * verticalSpacing) / yDivisions;
  float horizontalStep  = (1.0 - leftMargin - rightMargin - (xDivisions - 1) * horizontalSpacing) / xDivisions;

  float yLow, yUp, padBottomMargin, padTopMargin;
  float xLow, xUp, padLeftMargin, padRightMargin;

  for (int i = 0; i < xDivisions; ++i) {
    if (i == 0) {
      xLow = 0.0;
      xUp = leftMargin + horizontalStep;
      padLeftMargin = leftMargin / (xUp - xLow);
      padRightMargin = 0.0;
    } else if (i == xDivisions - 1) {
      xLow = xUp + horizontalSpacing;
      xUp = xLow + horizontalStep + rightMargin;
      padLeftMargin = 0.0;
      padRightMargin = rightMargin / (xUp - xLow);
    } else {
      xLow = xUp + horizontalSpacing;
      xUp = xLow + horizontalStep;
      padLeftMargin = 0.0;
      padRightMargin = 0.0;
    }

    for (int j = 0; j < yDivisions; ++j) {
      if (j == 0) {
        yLow = 0.0;
        yUp = bottomMargin + verticalStep;
        padBottomMargin = bottomMargin / (yUp - yLow);
        padTopMargin = 0.0;
      } else if (j == yDivisions - 1) {
        yLow = yUp + verticalSpacing;
        yUp = yLow + verticalStep + topMargin;
        padBottomMargin = 0.0;
        padTopMargin = topMargin / (yUp - yLow);
      } else {
        yLow = yUp + verticalSpacing;
        yUp = yLow + verticalStep;
        padBottomMargin = 0.0;
        padTopMargin = 0.0;
      }

      canvas->cd(0);

      const char* name = Form("pad_%i_%i", i, j);
      TPad* pad = dynamic_cast<TPad*>(gROOT->FindObject(name));
      if (pad)
        delete pad;

      pad = new TPad(name, "", xLow, yLow, xUp, yUp);
      pad->SetLeftMargin(padLeftMargin);
      pad->SetRightMargin(padRightMargin);
      pad->SetTopMargin(padTopMargin);
      pad->SetBottomMargin(padBottomMargin);

      pad->SetFrameBorderMode(0);
      pad->SetBorderMode(0);
      pad->SetBorderSize(0);

      pad->Draw();
    }
  }
}

#endif
