#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


////////////////////////////////////////////////////
// functions used to subset matrices and vectors
// All called from within C++ functions. 

Rcpp::NumericMatrix submat_currentScan(Rcpp::NumericMatrix X, Rcpp::NumericVector T, double testValue) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::mat y = Xmat.rows(arma::find(tIdx<testValue+0.0001 && tIdx>testValue-0.0001));
  return Rcpp::wrap(y);
}

Rcpp::NumericMatrix submat_intensities(Rcpp::NumericMatrix X, Rcpp::NumericVector T, double testValue, int testValue2) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::mat y = Xmat.rows(arma::find(tIdx > (testValue - testValue2/2) && tIdx < testValue+testValue2/2));
  return Rcpp::wrap(y);
}

Rcpp::NumericMatrix submat_leftIntens(Rcpp::NumericMatrix X, Rcpp::NumericVector T, double lastMz, double halfWindow, double mz) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::mat y = Xmat.rows(arma::find(tIdx > lastMz - halfWindow && (tIdx < mz - halfWindow || (tIdx > mz - halfWindow - 0.00001 && tIdx < mz - halfWindow + 0.00001))));
  return Rcpp::wrap(y);
}

Rcpp::NumericMatrix submat_rightIntens(Rcpp::NumericMatrix X, Rcpp::NumericVector T, double testValue, double testValue2, double testValue3) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::mat y = Xmat.rows(arma::find((tIdx >testValue+(testValue2) || (tIdx<testValue+(testValue2)+0.00001 && tIdx>testValue+(testValue2)-0.00001)) && tIdx < testValue3+(testValue2)));
  return Rcpp::wrap(y);
}

Rcpp::NumericVector subset_vector(Rcpp::NumericMatrix X, int column){
  NumericVector Z = X( _, column);
  return Z;
}

template <typename T>
inline bool rows_equal(const T& lhs, const T& rhs, double tol = 0.00000001) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

Rcpp::NumericVector unique_rowsVec(NumericVector X) {
  arma::vec Xvec(X.begin(), X.size(), false);
  unsigned int count = 1, i = 1, j = 1, nr = Xvec.n_elem;
  arma::vec result(nr);
  result.row(0) = Xvec.row(0);
  
  for ( ; i < nr; i++) {
    bool matched = false;
    if (rows_equal(Xvec.row(i), result.row(0))) continue;
    
    for (j = i + 1; j < nr; j++) {
      if (rows_equal(Xvec.row(i), Xvec.row(j))) {
        matched = true;
        break;
      }
    }
    if (!matched) result.row(count++) = Xvec.row(i);
  }
  arma::mat returnVec = result.rows(0, count - 1);
  return Rcpp::wrap(returnVec);
}

Rcpp::NumericMatrix remove_row( Rcpp::NumericMatrix X, arma::uword e) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  Xmat.shed_row(e);
  return Rcpp::wrap(Xmat);
}


Rcpp::NumericMatrix submat_singleScan(Rcpp::NumericMatrix X, Rcpp::NumericVector T, double testValue) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::mat y = Xmat.rows(arma::find(tIdx < testValue+0.00001 && tIdx > testValue-0.00001));
  return Rcpp::wrap(y);
}

Rcpp::NumericMatrix submat_reduced(Rcpp::NumericMatrix X, Rcpp::NumericVector T, double testValue, double testValue2) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::mat y = Xmat.rows(arma::find(tIdx > testValue2 && tIdx<testValue+10.0f));
  return Rcpp::wrap(y);
}

Rcpp::NumericMatrix subset_reducedMzCollectedMZ(Rcpp::NumericMatrix X, Rcpp::NumericVector T, Rcpp::NumericVector P, double mzRange0, double mzRange1, double mzTolerance, double rtRange1, double rtTolerance) {
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::colvec pIdx(P.begin(), P.size(),false);
  arma::mat y = Xmat.rows(arma::find((tIdx > (mzRange0 - (mzTolerance * 5))) && (tIdx < (mzRange1 + (mzTolerance * 5))) &&(pIdx < (rtRange1 + rtTolerance))));
  return Rcpp::wrap(y);
}

Rcpp::NumericMatrix submat_final(Rcpp::NumericMatrix X, int start, Rcpp::LogicalVector condition) { 
  int n=X.nrow(), k=X.ncol();
  int maxsize=0;

  for (int i = start, j = 0; i < n ; i++) {
    if(condition[i]){
      j = j+1;
      maxsize = maxsize+1;
    }
  }
  Rcpp::NumericMatrix out(maxsize,k);

  for (int i = start, j=0; i<n; i++) {
    if(condition[i]){
      for(int p = 0; p < k; p++  ){
        out(j,p) = X(i,p);
      }
      j = j+1;
    }
  }
  return(out);
}

//@title addScansRT 
//@description adds the retention times to the mass table for each entry in the run
//@param ionsTable the mz and intensities table
//@param scans - the list with each rt and scan for each mz.intensity entry to be added
//@return ionsTable

// [[Rcpp::export]]
Rcpp::NumericMatrix addScansCPP(Rcpp::NumericMatrix iontable, Rcpp::NumericMatrix scans){
  int iontableMaxRow = iontable.nrow();
  int scansIndex = 0;
  int scansValue = scans(scansIndex+1,0);
  
  for(int i = 0; i < iontableMaxRow; i++){
    if(iontable(i,2)<scansValue || scansIndex+1 ==scans.nrow()){
      iontable(i,2) = scans(scansIndex,1);
    }else{
      scansIndex++;
      scansValue = scans(scansIndex+1,0);
      iontable(i,2) = scans(scansIndex,1);
    }
  }
  colnames(iontable) = Rcpp::CharacterVector::create("intensity", "mz", "rt");
  return iontable;
}

//@title removeBackgroundCPP
//@description gets removes all ions below a sliding window cutoff
//@param ionsTable the input ions table with intensity mz and rt
//@param window.Size the sliding window size in mz
//@param aboveAverageMultiplier how many times above the background an ion must be
//@param minimumIntensity lowest intensity ions to consider
//@param mzRange range of ion masses to consider
//@return significantIonsTable input table but without ions that didn't pass filter

// [[Rcpp::export]]
Rcpp::NumericMatrix removeBackgroundCPP(Rcpp::NumericMatrix ionsTable, double windowSize, double aboveAverageMultiplier, double minIntensity, Rcpp::NumericVector mzRange){
  int nrows = ionsTable.nrow();
  int ncols = ionsTable.ncol();
  Rcpp::NumericMatrix significantIonsTable(nrows,ncols);
  int ionsCount = 0;
  double lastMz = 0.0f;
  int j = 0;
  double intenAvgLog = 0.0f;
  double oldRT = 0.0f;
  double halfWindow = windowSize/2;
  Rcpp::NumericVector armaTestVec = subset_vector(ionsTable, 2);
  Rcpp::NumericMatrix currentScan = submat_currentScan(ionsTable, armaTestVec, oldRT);
  
  for(int i = 0; i<nrows; i++){
    double currentRT = ionsTable(i,2);
    double mz = ionsTable(i,1);
    double intens = ionsTable(i,0);
    if(mz < mzRange(0) || mz > mzRange(1)){
      
    }else{
      
      if(oldRT != currentRT){
        currentScan = submat_currentScan(ionsTable, armaTestVec, currentRT);
        Rcpp::NumericVector mzIntensites = subset_vector(currentScan,1);
        Rcpp::NumericMatrix intesMat = submat_intensities(currentScan,mzIntensites, mz, windowSize);
        Rcpp::NumericVector intensities = subset_vector(intesMat, 0);
        intenAvgLog = mean(log(intensities));
        ionsCount = intensities.size();
        oldRT = currentRT;
      }else{
        Rcpp::NumericVector mzIntensities = subset_vector(currentScan,1);
        Rcpp::NumericMatrix leftIonsIntensityMatrix = submat_leftIntens(currentScan, mzIntensities, lastMz, halfWindow, mz);
        Rcpp::NumericVector leftIonsIntensities = subset_vector(leftIonsIntensityMatrix,0);
        Rcpp::NumericMatrix rightIonsIntensityMatrix = submat_rightIntens(currentScan, mzIntensities, lastMz, halfWindow, mz);
        Rcpp::NumericVector rightIonsIntensities = subset_vector(rightIonsIntensityMatrix,0);
        int ionsCountCurrent = ionsCount + rightIonsIntensities.size() - leftIonsIntensities.size();
        intenAvgLog = (intenAvgLog*ionsCount + sum(log(rightIonsIntensities)) - sum(log(leftIonsIntensities)))/ionsCountCurrent;
        ionsCount = ionsCountCurrent;
      }
      lastMz = mz;
      if((intens > minIntensity) && (log(intens) > (intenAvgLog + log(aboveAverageMultiplier)))){
        significantIonsTable(j,_) = ionsTable(i,_);
        j = j+1;
      }
    }
  } 
  Rcpp::LogicalVector finalCondition = significantIonsTable(_,0) > 0;
  Rcpp::NumericMatrix significantIonsTable2 = submat_final(significantIonsTable, 0, finalCondition);
  colnames(significantIonsTable2) = Rcpp::CharacterVector::create("intensity", "mz", "rt");
  return significantIonsTable2;
}

//@title groupIonsSameScanCPP
//@description groups MZ ions that are 1 mass apart in the same scan
//@param ionsTable
//@param mzTolerance is the mass value that the different peaks can differ by
//@return mzIntensities.mzCollected

// [[Rcpp::export]]
Rcpp::NumericVector groupIonsSameScanCPP(Rcpp::NumericMatrix ionsTable, double mzTolerance){
  int currentRow = 0;
  Rcpp::NumericMatrix mzIntensitiesMzCollected(ionsTable.nrow(),5);
  Rcpp::NumericVector retentionTimes = subset_vector(ionsTable, 2);
  Rcpp::NumericVector RTs = unique_rowsVec(retentionTimes);

  for(int i = 0; i<RTs.size(); i++){
    double rt = RTs(i);
    Rcpp::NumericMatrix singleScan = submat_singleScan(ionsTable, retentionTimes, rt);

    while(singleScan.nrow()>0){
      Rcpp::NumericMatrix rowQuery(1,3);
      rowQuery(0,_) = singleScan(0,_);
      arma::uword remove = 0;
      singleScan = remove_row(singleScan, remove);
      double intensityQuery = rowQuery(0,0);
      double mzQuery = rowQuery(0,1);
      double mzHighest = mzQuery;
      double intensityMax = intensityQuery;
      double mzMaxIntensity = mzQuery;
      Rcpp::NumericVector singleScanMZ = subset_vector(singleScan,1);
      double mzQT = mzQuery-mzTolerance;
      Rcpp::NumericMatrix reducedTable = submat_reduced(singleScan, singleScanMZ, mzQuery, mzQT);

      if(!NumericVector::is_na(reducedTable(0,0))){
        int toDelete [reducedTable.rows()];
        int arrayIndex = 0;

        for(int j = 0; j< reducedTable.rows(); j++){
          Rcpp::NumericMatrix rowNext(1,3);
          rowNext(0,_) = reducedTable(j,_);
          double intensityNext = rowNext(0,0);
          double mzNext = rowNext(0,1);
          double nextHighest = fabs(mzNext - mzHighest);

          if(((nextHighest>(1-mzTolerance*2)-0.00001) && (nextHighest<(1+mzTolerance*2)+0.00001))
            || ((nextHighest>(0.5-mzTolerance*2)-0.00001) && (nextHighest<(0.5+mzTolerance*2)+0.00001))
            || ((nextHighest>(0.33-mzTolerance*2)-0.00001) && (nextHighest<(0.33+mzTolerance*2)+0.00001))
            || ((nextHighest>(0.25-mzTolerance*2)-0.00001) && (nextHighest<(0.25+mzTolerance*2)+0.00001))){
              mzHighest = mzNext;
              if(intensityNext > intensityMax){
                intensityMax = intensityNext;
                mzMaxIntensity = mzNext;
              }
              toDelete[arrayIndex] = j;
              arrayIndex++;              
          }
        }
        for(int p = arrayIndex-1; p > -1; p--){
          singleScan = remove_row(singleScan, toDelete[p]);
        }
      }
      mzIntensitiesMzCollected(currentRow,0) = intensityMax;
      mzIntensitiesMzCollected(currentRow,1) = mzMaxIntensity;
      mzIntensitiesMzCollected(currentRow,2) = rt;
      mzIntensitiesMzCollected(currentRow,3) = mzQuery;
      mzIntensitiesMzCollected(currentRow,4) = mzHighest;
      currentRow = currentRow +1;
    }
  }
  Rcpp::LogicalVector finalCondition = mzIntensitiesMzCollected(_,0) > 0;
  Rcpp::NumericMatrix mzIntensitiesMzCollected2 = submat_final(mzIntensitiesMzCollected, 0, finalCondition);
  colnames(mzIntensitiesMzCollected2) = Rcpp::CharacterVector::create("highestIntensity", "mz", "rt", "minMZ", "maxMZ");
  return mzIntensitiesMzCollected2;
}

//@title collectFinalPeaks
//@description group overlapping peaks together
//@param ionsDataFrame that has passed through the mass and time collections
//@param mzTolerance is the mass value that the different peaks can differ by to be grouped
//@param rtTolerance is the retention time value that the different peaks can differ by to be grouped in seconds
//@return reduced

// [[Rcpp::export]]
Rcpp::NumericMatrix collectFinalPeaksCPP(Rcpp::NumericMatrix ionsDataFrame, double mzTolerance, double rtTolerance){
  NumericMatrix reduced = ionsDataFrame;

  int previous_nrow = 0;
  while(reduced.nrow() != previous_nrow){
    previous_nrow = reduced.nrow();
    for(int i = 0; i<reduced.nrow(); i++){
      Rcpp::NumericMatrix row(1,reduced.ncol());
      row(0,_) = reduced(i,_);
  
      if(i+1 == reduced.nrow()){
        break;
      }
      for(int j=i+1; j<reduced.nrow();j++){
        Rcpp::NumericMatrix other(1,reduced.ncol());
        other(0,_) = reduced(j,_);
  
        if(NumericVector::is_na(other(0,0))){
          break;
        }
        if(row(0,1) > other(0,1) - mzTolerance 
             && row(0,1) < other(0,1) + mzTolerance 
             && ((row(0,5) > other(0,5) - rtTolerance
                    && row(0,5) < other(0,5) + rtTolerance)
                ||(row(0,5) > other(0,6) - rtTolerance
                    && row(0,5) < other(0,6) + rtTolerance)
                ||(row(0,6) > other(0,6) - rtTolerance
                    && row(0,6) < other(0,6) + rtTolerance)
                ||(row(0,6) > other(0,5) - rtTolerance
                    && row(0,6) < other(0,5) + rtTolerance))){
          if(row(0,5) > other(0,5)){
            reduced(i,5) = other(0,5);
          }
          if(row(0,6) < other(0,6)){
            reduced(i,6) = other(0,6);
          }
          reduced = remove_row(reduced, j);
        }
      }
    }
  }
  colnames(reduced) = Rcpp::CharacterVector::create("highestIntensity","mz","rt","minMZ","maxMZ","minRT","maxRT","RTfraction","NumScanRange");
  return(reduced);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// sameMZCPP and groupPeaksAcrossScansCPP are currently not being used. The R version was one of the
// faster functions and the data from the C++ version never matched to the R version. More effort
// would be needed to fix up these functions to get them to produce identical output to the R version
//////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
bool sameMZCPP(double mzRange10, double mzRange11, double mzRange20, double mzRange21, double mzTolerance){
  bool overlap = FALSE;

  if((mzRange10-mzTolerance) < mzRange20 && (mzRange11+mzTolerance) > mzRange20
       ||(mzRange10-mzTolerance) < mzRange21 && (mzRange11+mzTolerance) > mzRange21
       ||(mzRange10+mzTolerance) > mzRange20 && (mzRange11-mzTolerance) < mzRange21){
    double intpart = 0.0 ;
    double decimal1 = modf(mzRange10, &intpart);
    double decimal2 = modf(mzRange20, &intpart);
   
    if(fabs(decimal1 - decimal2) < mzTolerance*2 || (fabs(decimal1 - decimal2) < mzTolerance*2+0.0001 && fabs(decimal1 - decimal2) > mzTolerance*2+0.0001)){
      overlap = TRUE;
    }
  }

  return(overlap);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// NOT FUNCTIONAL AND NOT USED
// TODO: GET FUNCTIONAL, WILL RESULT IN LARGE SPEED UP
//////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::NumericMatrix groupPeaksAcrossScansCPP(Rcpp::NumericMatrix ionsDataFrame, double mzTolerance, double rtTolerance, double meanTimeBetweenScans){
  int currentRow = 0;
  Rcpp::NumericMatrix mzIntensitiesRTCollected(ionsDataFrame.nrow(),9);

  while(ionsDataFrame.nrow()>0){
    int toDelete [ionsDataFrame.rows()];
    int arrayIndex = 0;
    Rcpp::NumericMatrix rowQuery(1, ionsDataFrame.ncol());
    rowQuery(0,_) = ionsDataFrame(0,_);
    ionsDataFrame = remove_row(ionsDataFrame, 0);
    double rtRange[2] = {0,0};
    rtRange[0] = rowQuery(0,2);
    rtRange[1] = rowQuery(0,2);
    double mzRange[2] = {0,0};
    mzRange[0] = rowQuery(0,3);
    mzRange[1] = rowQuery(0,4);
    double intensMzRt[3] = {0,0,0};
    intensMzRt[0] = rowQuery(0,0);
    intensMzRt[1] = rowQuery(0,1);
    intensMzRt[2] = rowQuery(0,2);
    int rtCount = 1;
    Rcpp::NumericVector minMzVector = subset_vector(ionsDataFrame, 3);    
    Rcpp::NumericVector RtVector  = subset_vector(ionsDataFrame, 2);
    Rcpp::NumericMatrix reducedMzCollected = subset_reducedMzCollectedMZ(ionsDataFrame, minMzVector, RtVector, mzRange[0], mzRange[1], mzTolerance,rtRange[1], rtTolerance);
    
    if(reducedMzCollected(0,0) > 0){
      for(int j=0; j<reducedMzCollected.nrow(); j++){
        Rcpp::NumericMatrix rowNext(1,reducedMzCollected.ncol());
        rowNext(0,_) = reducedMzCollected(j,_);
        double intensityNext = rowNext(0,0);
        double mzIntensNext = rowNext(0,1);
        double rtNext = rowNext(0,2);
        double mzRangeNext[2] = {0,0};
        mzRangeNext[0] = rowNext(0,3);
        mzRangeNext[1] = rowNext(0,4);
        
        if(sameMZCPP(mzRange[0],mzRange[1], mzRangeNext[0], mzRangeNext[1],mzTolerance)){
          rtRange[1] = rtNext;
          
          if(mzRangeNext[0]<mzRange[0]){
            mzRange[0]= mzRangeNext[0];
          }
          
          if(mzRangeNext[1]<mzRange[1]){
            mzRange[1]= mzRangeNext[1];
          }
          
          if(intensityNext > intensMzRt[0]){
            intensMzRt[0] = intensityNext;
            intensMzRt[1] = mzIntensNext;
            intensMzRt[2] = rtNext;
          }
          toDelete[arrayIndex] = j;
          arrayIndex++;
          rtCount = rtCount + 1;
        }
        arma::mat Xmat(ionsDataFrame.begin(), ionsDataFrame.nrow(), ionsDataFrame.ncol(), false);
        arma::mat Ymat(reducedMzCollected.begin(), reducedMzCollected.nrow(), reducedMzCollected.ncol(), false);
        
        for(int p = arrayIndex-1; p > -1; p--){
          for(int k = 0; k <ionsDataFrame.nrow(); k++){
            if(rows_equal(Xmat.row(k),Ymat.row(p))){
              ionsDataFrame = remove_row(ionsDataFrame, k);
              break;
            }
          }
        }
      }
    }
    double numScans = ((rtRange[1] - rtRange[0])/meanTimeBetweenScans)+1;
    double rtFractiontemp = rtCount/numScans;
    double rtFraction = floor(rtFractiontemp*100)/100;
    mzIntensitiesRTCollected(currentRow, 0) = intensMzRt[0];
    mzIntensitiesRTCollected(currentRow, 1) = intensMzRt[1];
    mzIntensitiesRTCollected(currentRow, 2) = intensMzRt[2];
    mzIntensitiesRTCollected(currentRow, 3) = mzRange[0];
    mzIntensitiesRTCollected(currentRow, 4) = mzRange[1];
    mzIntensitiesRTCollected(currentRow, 5) = rtRange[0];
    mzIntensitiesRTCollected(currentRow, 6) = rtRange[1];
    mzIntensitiesRTCollected(currentRow, 7) = rtFraction;
    mzIntensitiesRTCollected(currentRow, 8) = rtCount;
    currentRow = currentRow + 1;
  }
  Rcpp::LogicalVector finalCondition = mzIntensitiesRTCollected(_,0) > 0;
  Rcpp::NumericMatrix mzIntensitiesRTCollected2 = submat_final(mzIntensitiesRTCollected, 0, finalCondition);
  colnames(mzIntensitiesRTCollected2) = Rcpp::CharacterVector::create("highestIntensity", "mz", "rt", "minMZ", "maxMZ", "minRT", "maxRT", "RTfraction", "NumScanRange");
  return mzIntensitiesRTCollected2;
}