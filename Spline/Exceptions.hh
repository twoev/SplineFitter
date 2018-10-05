#ifndef SPLINE_EXCEPTIONS_HH
#define SPLINE_EXCEPTIONS_HH

#include <stdexcept>

namespace Spline{

  //////////////////////////////////////////////////////////////////////////////
  class EvaluationException: public std::runtime_error{
    
  public:
    EvaluationException(): std::runtime_error
    ("Unable to evaluate b-spline: control points or coefficients not set!"){};
    
  };
  
  //////////////////////////////////////////////////////////////////////////////  
  class BasisException: public std::runtime_error{
    
  public:
    BasisException(): std::runtime_error
    ("Unable to evaluate basis functions: knot vector not set!"){};
    
  };

  //////////////////////////////////////////////////////////////////////////////  

  class BadDefinition: public std::runtime_error{
  public:
    BadDefinition(): std::runtime_error
    ("Spline definition is bad"){};
  };
  
  //////////////////////////////////////////////////////////////////////////////  
  class MatrixException: public std::runtime_error{
  public:
    MatrixException(): std::runtime_error
    ("row or column index is beyond matrix bounds!"){};
  };
  
}
#endif