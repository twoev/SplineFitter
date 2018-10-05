#include "Spline/BSpline.hh"
#include "Spline/Exceptions.hh"

#define BOOST_UBLAS_NDEBUG 1 

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace Spline{
  
  using namespace boost::numeric;
  
  BSpline::BSpline(size_t degree): m_degree(degree), m_canEvaluate(false),
  m_gslCoefficients(0), m_covarianceMatrix(0), m_basisFuncCoeffs(0), m_haveSignalTemplate(false),
  m_covariance(0,0),
  m_bSplineWorkspace(0){
    
  }
  
  DataPoint BSpline::evaluate(double affine)const{
    
    if(!m_canEvaluate) throw EvaluationException();
        
    double y, er;
    
    gsl_bspline_eval(affine, _basisFuncCoeffs(), m_bSplineWorkspace);
    
    if(m_haveSignalTemplate){
      gsl_vector_view coefficients = gsl_vector_subvector(m_gslCoefficients, 0, nCoefficients());
      gsl_matrix_view covariance = gsl_matrix_submatrix(m_covarianceMatrix, 0, 0, nCoefficients(), nCoefficients());
      gsl_multifit_linear_est(_basisFuncCoeffs(), &(coefficients.vector), &(covariance.matrix), &y, &er);
    }else{
      gsl_multifit_linear_est(_basisFuncCoeffs(), m_gslCoefficients, m_covarianceMatrix, &y, &er);
    }
    
    DataPoint result(affine, y, er);
    
    return result;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  vector<double> BSpline::basisFunctions(double affine)const{
   
    if(!m_bSplineWorkspace){
      throw BasisException();
    }
    
    gsl_bspline_eval(affine, _basisFuncCoeffs(), m_bSplineWorkspace);
    vector<double> result;
    
    for(size_t ii=0; ii != nCoefficients(); ++ii){
      result.push_back(gsl_vector_get(_basisFuncCoeffs(), ii));
    }
    
    return result;
  }

  //////////////////////////////////////////////////////////////////////////////
  double BSpline::basisFunction(size_t index, double affine)const{
    
    if(!m_bSplineWorkspace){
      throw BasisException();
    }
    
    gsl_bspline_eval(affine, _basisFuncCoeffs(), m_bSplineWorkspace);
    
    return gsl_vector_get(_basisFuncCoeffs(), index);
  }
  
  //////////////////////////////////////////////////////////////////////////////
  double BSpline::chi2()const{
    return m_chi2;
  }
  
  //////////////////////////////////////////////////////////////////////////////  
  double BSpline::chi2_dof()const{
    return m_chi2_dof;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  const vector<double> &BSpline::coefficients()const{
    return m_coefficients;
  }
  
  //////////////////////////////////////////////////////////////////////////////  
  double BSpline::signalCoefficientSize()const{
    if(!m_haveSignalTemplate)throw std::runtime_error("Signal coefficient size not available - no signal template was used in the fit!");
    return m_coefficients.back();
  }
  
  //////////////////////////////////////////////////////////////////////////////  
  const Matrix<double> &BSpline::covarianceMatrix()const{
    return m_covariance;
  }
  
  //////////////////////////////////////////////////////////////////////////////  
  size_t BSpline::nCoefficients()const{
    return m_nCoeffs;
  }

  //////////////////////////////////////////////////////////////////////////////  
  const PointSet &BSpline::controlPoints()const{
   
    size_t nDiff = nCoefficients() - knots().size();
    
    size_t t = 0;
    
    ublas::vector<double> xt(nCoefficients());
    
    for(vector<double>::const_iterator k = knots().begin(); k != knots().end(); ++k){
      
      xt(t) = *k;
      ++t;
      if(t < 2 * nDiff){
        vector<double>::const_iterator kp = k;
        ++kp;
        xt(t) = 0.5*(*k + *kp);
        ++t;
      }
    }
    
    ublas::matrix<double> basisMatrix(nCoefficients(), nCoefficients());
    ublas::vector<double> yt(nCoefficients());
    
    for(size_t t=0; t != nCoefficients(); ++t){
      const vector<double> &bFuncs = basisFunctions(xt(t));
      
      double y, er;
      
      if(m_haveSignalTemplate){
        gsl_vector_view coefficients = gsl_vector_subvector(m_gslCoefficients, 0, nCoefficients());
        gsl_matrix_view covariance = gsl_matrix_submatrix(m_covarianceMatrix, 0, 0, nCoefficients(), nCoefficients());
        gsl_multifit_linear_est(_basisFuncCoeffs(), &(coefficients.vector), &(covariance.matrix), &y, &er);
      }else{
        gsl_multifit_linear_est(_basisFuncCoeffs(), m_gslCoefficients, m_covarianceMatrix, &y, &er);
      }      
      yt(t) = y;
      
      for(size_t i=0; i != nCoefficients(); ++i){
        basisMatrix(t,i) = bFuncs[i];
      }
      
    }
    
    ublas::matrix<double> toTest(basisMatrix);
        
    // here comes some boost-voodoo for inverting the basisMatrix
    
    ublas::permutation_matrix<size_t> pmatrix(basisMatrix.size1());
    int res = ublas::lu_factorize(basisMatrix, pmatrix);
    
    if(res != 0) throw std::runtime_error("Unable to invert basis matrix for determining control points!");
    
    ublas::matrix<double> invBasis(ublas::identity_matrix<double>(basisMatrix.size1()));
    ublas::lu_substitute(basisMatrix, pmatrix, invBasis);
    
    ublas::matrix<double> test = prod(toTest, invBasis) - ublas::identity_matrix<double>(basisMatrix.size1());
    
    /**
     * Here we test the inversion of the matrix.
     * Since larger matrices can suffer from numerical instabilities, we lower the test threshold 
     * (hence making failure less likely) for larger matrices.  Nevertheless, it may still fail
     * for larger matrices, or for those with difficult knot vectors or sampling points.
     *
     */
    
    double thresh = 1.e-7;
    
    if(nCoefficients() > 28){
      thresh = 1.e-3;
    }else if(nCoefficients() > 20){
      thresh = 1.e-5;
    }else if(nCoefficients() > 15){
      thresh = 1.e-6;
    }
    
    for(size_t i=0; i != test.size1(); ++i){
      for(size_t j=0; j !=test.size2(); ++j){
        if(fabs(test(i,j)) > thresh){
          
          std::cout<<"Matrix inversion failed."<<std::endl;
          std::cout<<"This is probably due to numerical instabilities in the inversion if you have a large (~30) number of coefficients or knots"<<std::endl;
          std::cout<<"It may also happen if you have an especially pathalogical knot vector!"<<std::endl;
          std::cout<<"Useful debug info follows..."<<std::endl;
          
          std::cout<<"Basis matrix is..."<<std::endl;
          
          for(size_t m=0; m != toTest.size1(); ++m){
            for(size_t n=0; n !=toTest.size2(); ++n){
              std::cout.width(12);
              std::cout<<std::fixed<<toTest(m,n);
              std::cout.fill(' ');
            }
            std::cout<<std::endl;
          }
          
          std::cout<<std::endl<<"B.B^-1 is..."<<std::endl;
          
          for(size_t m=0; m != test.size1(); ++m){
            for(size_t n=0; n !=test.size2(); ++n){
              std::cout.width(12);
              std::cout<<std::fixed<<test(m,n);
              std::cout.fill(' ');
            }
            std::cout<<std::endl;
          }
          
          throw std::runtime_error("Inverse basis matrix has numerical errors!");
        }
      }
    }
    
    ublas::vector<double> xi = ublas::prod(invBasis, xt);
    ublas::vector<double> yi = ublas::prod(invBasis, yt);
    
    m_controlPoints.clear();
    
    for(size_t i=0; i!=xi.size(); ++i){
      m_controlPoints.push_back(Point(xi(i), yi(i)));
    }
    
    return m_controlPoints;
  }
  
  //////////////////////////////////////////////////////////////////////////////  
  const vector<double> &BSpline::knots()const{
    return m_knots;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  gsl_vector *BSpline::_basisFuncCoeffs()const{
        
    if(! m_basisFuncCoeffs){
      m_basisFuncCoeffs = gsl_vector_alloc(nCoefficients());
    }
    
    return m_basisFuncCoeffs;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  void BSpline::setNCoefficientsUniform(size_t nCoeffs, double min, double max){
   
    m_nCoeffs = nCoeffs;
    size_t nbreak = nCoeffs + 2 - m_degree;
    
    if(m_bSplineWorkspace){
      gsl_bspline_free(m_bSplineWorkspace);
    }
    
    m_bSplineWorkspace = gsl_bspline_alloc(m_degree, nbreak);
    
    gsl_bspline_knots_uniform(min, max, m_bSplineWorkspace);
    
    gsl_vector *knots = m_bSplineWorkspace->knots;
    
    m_knots.clear();
    
    for(size_t i= m_degree - 1; i < knots->size + 1 - m_degree; ++i){
      m_knots.push_back(gsl_vector_get(knots, i));
    }
    
    m_canEvaluate = false;
    
    return;
  }
 
  //////////////////////////////////////////////////////////////////////////////
  void BSpline::setKnotVector(const vector<double> &knots){
  
    m_nCoeffs = knots.size() + m_degree - 2;
    
    m_knots = knots;
    
    gsl_vector *breakPts = gsl_vector_alloc(knots.size());
    
    for(size_t i=0; i != knots.size(); ++i){
      gsl_vector_set(breakPts, i, knots[i]);
    }
    
    if(m_bSplineWorkspace){
      gsl_bspline_free(m_bSplineWorkspace);
    }
    
    m_bSplineWorkspace = gsl_bspline_alloc(m_degree, knots.size());
        
    gsl_bspline_knots(breakPts, m_bSplineWorkspace);
    
    gsl_vector_free(breakPts);
    
    if(m_basisFuncCoeffs){
      gsl_vector_free(m_basisFuncCoeffs);
      m_basisFuncCoeffs = 0;
    }
    
    m_canEvaluate = false;
    
    /*
    for(size_t i=0; i != m_bSplineWorkspace->knots->size; ++i){
      std::cout<<gsl_vector_get(m_bSplineWorkspace->knots, i)<<"  ";
    }
    */    
    return;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  void BSpline::define(const vector<double> &knots, const vector<double> &coefficients){
    
    setKnotVector(knots);
    
    if(coefficients.size() != nCoefficients()) throw BadDefinition();
    
    m_coefficients = coefficients;
    
    if(m_gslCoefficients){
      gsl_vector_free(m_gslCoefficients);
    }
    
    if(m_covarianceMatrix){
      gsl_matrix_free(m_covarianceMatrix);
    }
    
    m_gslCoefficients = gsl_vector_alloc(nCoefficients());
    m_covarianceMatrix = gsl_matrix_alloc(nCoefficients(), nCoefficients());

    for(size_t ii=0; ii != nCoefficients(); ++ii){
      gsl_vector_set(m_gslCoefficients, ii, m_coefficients[ii]);
      for(size_t jj = 0; jj != nCoefficients(); ++jj){
        gsl_matrix_set(m_covarianceMatrix, ii, jj, 0);
      }
    }
    
    m_canEvaluate = true;
    m_haveSignalTemplate = false;
    
    return;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  void BSpline::define(const vector<double> &knots,
                       const vector<double> &coefficients, 
                       const Matrix<double> &covariance){
    
    define(knots, coefficients);
    
    m_covariance = covariance;
    
    for(size_t ii=0; ii != nCoefficients(); ++ii){
      for(size_t jj = 0; jj != nCoefficients(); ++jj){
        gsl_matrix_set(m_covarianceMatrix, ii, jj, covariance.element(ii, jj));
      }
    }
    
    return;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  void BSpline::fitLeastSquares(const DataSet &data){
   
    if(!m_bSplineWorkspace){
      throw BasisException();
    }
    
    gsl_vector *basisCoeffs = gsl_vector_alloc(nCoefficients());
    gsl_matrix *bases = gsl_matrix_alloc(data.size(), nCoefficients());
    
    gsl_vector *ys = gsl_vector_alloc(data.size());
    gsl_vector *weights = gsl_vector_alloc(data.size());
    
    for(size_t i = 0; i != data.size(); ++i){
      gsl_vector_set(ys, i, data[i].y());
      double w = 1. / (data[i].error());
      gsl_vector_set(weights, i, w*w);
      gsl_bspline_eval(data[i].x(), basisCoeffs, m_bSplineWorkspace);
      
      for(size_t j = 0; j != nCoefficients(); ++j){
        gsl_matrix_set(bases, i, j, gsl_vector_get(basisCoeffs, j));
      }
    }
    
    gsl_multifit_linear_workspace *fitWorkSpace = gsl_multifit_linear_alloc(data.size(), nCoefficients());
    
    if(m_gslCoefficients){
      gsl_vector_free(m_gslCoefficients);
    }
    
    if(m_covarianceMatrix){
      gsl_matrix_free(m_covarianceMatrix);
    }
    
    m_gslCoefficients = gsl_vector_alloc(nCoefficients());
    m_covarianceMatrix = gsl_matrix_alloc(nCoefficients(), nCoefficients());
    
    gsl_multifit_wlinear(bases, weights, ys, m_gslCoefficients, m_covarianceMatrix, &m_chi2, fitWorkSpace);
    
    m_chi2_dof = m_chi2 / ((double)data.size());
    
    gsl_multifit_linear_free(fitWorkSpace);
    
    gsl_vector_free(basisCoeffs);
    gsl_matrix_free(bases);
    gsl_vector_free(ys);
    gsl_vector_free(weights);
    
    m_coefficients.clear();
    m_covariance = Matrix<double>(nCoefficients(), nCoefficients());
    
    for(size_t i = 0; i != nCoefficients(); ++i){
      m_coefficients.push_back(gsl_vector_get(m_gslCoefficients, i));
      for(size_t j = 0; j != nCoefficients(); ++j){
        m_covariance.setElement(i, j, gsl_matrix_get(m_covarianceMatrix, i, j));
      }
    }
    
    
    
    m_canEvaluate = true;
    m_haveSignalTemplate = false;
    
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  
  void BSpline::fitLeastSquares(const DataSet &data, const DataSet &signalTemplate){
    
    if(!m_bSplineWorkspace){
      throw BasisException();
    }
    
    gsl_vector *basisCoeffs = gsl_vector_alloc(nCoefficients());
    gsl_matrix *bases = gsl_matrix_alloc(data.size(), nCoefficients() + 1);
    
    gsl_vector *ys = gsl_vector_alloc(data.size());
    gsl_vector *weights = gsl_vector_alloc(data.size());
    
    for(size_t i = 0; i != data.size(); ++i){
      gsl_vector_set(ys, i, data[i].y());
      double w = 1. / (data[i].error());
      gsl_vector_set(weights, i, w*w);
      gsl_bspline_eval(data[i].x(), basisCoeffs, m_bSplineWorkspace);
      
      for(size_t j = 0; j != nCoefficients(); ++j){
        gsl_matrix_set(bases, i, j, gsl_vector_get(basisCoeffs, j));
      }
      
      gsl_matrix_set(bases, i, nCoefficients(), signalTemplate[i].y());
      
    }
    
    gsl_multifit_linear_workspace *fitWorkSpace = gsl_multifit_linear_alloc(data.size(), nCoefficients()+1);

    if(m_gslCoefficients){
      gsl_vector_free(m_gslCoefficients);
    }
    
    if(m_covarianceMatrix){
      gsl_matrix_free(m_covarianceMatrix);
    }
    
    m_gslCoefficients = gsl_vector_alloc(nCoefficients()+1);
    m_covarianceMatrix = gsl_matrix_alloc(nCoefficients()+1, nCoefficients()+1);
    
    gsl_multifit_wlinear(bases, weights, ys, m_gslCoefficients, m_covarianceMatrix, &m_chi2, fitWorkSpace);
    
    m_chi2_dof = m_chi2 / ((double)data.size());
    
    gsl_multifit_linear_free(fitWorkSpace);
    
    gsl_vector_free(basisCoeffs);
    gsl_matrix_free(bases);
    gsl_vector_free(ys);
    gsl_vector_free(weights);
    
    m_coefficients.clear();
    m_covariance = Matrix<double>(nCoefficients()+1, nCoefficients()+1);
    
    for(size_t i = 0; i != nCoefficients()+1; ++i){
      m_coefficients.push_back(gsl_vector_get(m_gslCoefficients, i));
      for(size_t j = 0; j != nCoefficients()+1; ++j){
        m_covariance.setElement(i, j, gsl_matrix_get(m_covarianceMatrix, i, j));
      }
    }
    
    m_canEvaluate = true;
    m_haveSignalTemplate = true;
    
    return;
  }
  
  
}

