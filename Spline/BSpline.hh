/**
 *  B spline fitting class using GSL
 *  
 *  @author James Monk <jmonk@cern.ch>
 * 
 */

#ifndef SPLINE_BSPLINE_HH
#define SPLINE_BSPLINE_HH

#include "Spline/DataPoint.hh"
#include "Spline/Matrix.hh"

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

#include <ostream>
#include <vector>

namespace Spline{
    
  using std::vector;
  
  class BSpline{
    
  public:
    
    /**
     *  Default constructor.
     *
     *  \param degree the degree of the splines.  Default 4
     */
    BSpline(size_t degree = 4);
    
    /**
     * Set the knot vector
     * This allows for non-uniform knot placement
     */
    void setKnotVector(const vector<double> &knots);
    
    /**
     *  Set the number of coefficients using a uniform knot vector in the range min to max
     *  The min and max values should tpically be min and max x range of the histogram to be fitted, 
     *  or within that range if only a subset of the histogram is to be used.
     *
     *  \param nCoeffs the number of coefficients
     *  \param min the minimum value of the affine parameter
     *  \param max the maximum value of the affine parameter
     *
     */
    void setNCoefficientsUniform(size_t nCoeffs, double min, double max);
    
    /**
     *  \return The number of coefficients
     */
    size_t nCoefficients()const;
    
    /**
     *  Set the knot vector and coefficients.  Used if the spline is to be defined rather than fitted to data
     */
    void define(const vector<double> &knots, const vector<double> &coefficients);
    
    /**
     *  Set the knot vector, coefficients and the covariance matrix of the coefficients.
     *  Used if the spline is to be defined rather than fitted to data
     */
    void define(const vector<double> &knots, const vector<double> &coefficients, const Matrix<double> &covariance);
    
    /**
     *  Evaluate one of the basis functions for a given affine parameter.
     *  For a normal 1D histogram the affine parameter is simply the distance along the x axis
     *  This method is less efficient than evaluating all of the functions simultaneously.
     *  \param index the label of the basis function
     *  \param affine the value of the affine parameter  
     *  
     *  \return the value of the basis function
     *
     */
    double basisFunction(size_t index, double affine)const;
    
    /**
     *  Evaluate all of the basis functions at a give affine parameter.
     *  For a normal 1D histogram the affine parameter is simply the distance along the x axis.
     *  This method is more efficient that evaluating each basis function individually
     *
     *  \param affine the affine parameter at which to evaluate the bases
     *  \return a vector of the function values
     */
    vector<double> basisFunctions(double affine)const;
    
    /**
     * Fit this B-spline curve to a set of data points
     */
    void fitLeastSquares(const DataSet &data);
    
    /**
     * Fit this B-spline curve to a set of data points and include a signal template in the fit
     */    
    void fitLeastSquares(const DataSet &data, const DataSet &signalTemplate);
    
    /**
     *  \return the chi squared of this fitted curve compared to the input data
     */
    double chi2()const;
    
    /**
     * \return the chi squared per degree of freedom of this fitted curve compared to the input data
     */
    double chi2_dof()const;
    
    /**
     *  Evaluate the curve at the give affine parameter
     *  For a normal 1D histogram the affine parameter is simply the distance along the x axis.
     * 
     *  \param affine the affine parameter
     *  \return DataPoint giving the x and y value of the curve and the error on the y value
     */
    DataPoint evaluate(double affine)const;
    
    
    /**
     *  \return The knot vector of this B-spline
     */
    const vector<double> &knots()const;
    
    /**
     * \return the coefficients describing this curve
     */
    const vector<double> &coefficients()const;
    
    /**
     * Give the estimate of the signal size if a signal template was included in the fit
     * Will throw an error if no signal was included!
     * \return the signal template contribution to the fit
     */
    double signalCoefficientSize()const;
    
    /**
     *  Give the covariance matrix of the fit coefficients
     *
     *  \return a N*N vector<vector<double> >, where N is the number of coefficients
     */
    const Matrix<double> &covarianceMatrix()const;
    
    /**
     *  The coefficients and knots describing this curve can be interpretted geometrically 
     *  as a set of control points describing a polygon.  
     *  The curve can be derived by interpolating between these control points using the de-boor algorithm
     * 
     *  \return The set of control points describing this curve
     */
    const PointSet &controlPoints()const;
    
  private:
    
    size_t m_nCoeffs;
    size_t m_nBreakPoints;
    size_t m_degree;
    
    bool m_canEvaluate;
    
    vector<double> m_knots;
    
    gsl_vector *m_gslCoefficients;
    gsl_matrix *m_covarianceMatrix;
    
    gsl_vector *_basisFuncCoeffs()const;
    mutable gsl_vector *m_basisFuncCoeffs;
    
    bool m_haveSignalTemplate;
    
    // The same coefficients in a more usable std::vector
    vector<double> m_coefficients;
    // The same covariance matrix in a more friendly std::vector<std::vector>
    Matrix<double> m_covariance;
    
    mutable PointSet m_controlPoints;
    
    double m_chi2;
    double m_chi2_dof;
    
    gsl_bspline_workspace *m_bSplineWorkspace;
    
  };
}

#endif

