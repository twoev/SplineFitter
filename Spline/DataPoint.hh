#ifndef SPLINE_DATAPOINT_HH
#define SPLINE_DATAPOINT_HH

#include <vector>

namespace Spline{
 
  using std::vector;
  
  //////////////////////////////////////////////////////////////////////////////
  class Point{
  public:
    Point(double x, double y): m_x(x), m_y(y){}
    
    double x()const{return m_x;}
    double y()const{return m_y;}
    
  private:
    
    double m_x;
    double m_y;
  };
  
  //////////////////////////////////////////////////////////////////////////////
  class DataPoint: public Point{
  public:
    DataPoint(double x, double y, double er): Point(x, y), m_err(er)
    {}
    
    double error()const{return m_err;}
    
  private:
    
    double m_err;
    
  };
  //////////////////////////////////////////////////////////////////////////////
  
  typedef vector<Point> PointSet;
  typedef vector<DataPoint> DataSet;
}

#endif
