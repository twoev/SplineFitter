#ifndef SPLINE_MATRIX_HH
#define SPLINE_MATRIX_HH

#include "Spline/Exceptions.hh"
#include <ostream>
#include <iostream>

namespace Spline{
  
  /**
   *  Basic matrix class to hold the covariance matrix
   */
  template <class T>
  class Matrix{
    
  public:
    
    /**
     *  Construct a matrix of size rows * columns
     *
     *  \param rows the number of rows
     *  \param columns the number of columns
     */
    Matrix(size_t rows, size_t columns): m_elements(rows, vector<T>(columns)),
    m_nRows(rows), m_nColumns(columns){}
    
    /**
     *  \param row the row of the element to return
     *  \param column the column of the element to return
     *  \return the element at row, column
     */
    T element(size_t row, size_t column)const{
      if(row >= m_nRows || column >= m_nColumns) throw MatrixException();
      return m_elements[row][column];
    };
    
    /**
     *  Set the element at location row, column
     *  \param row the row of the element
     *  \param column the column of the element
     */
    void setElement(size_t row, size_t column, T value){
      if(row >= m_nRows || column >= m_nColumns) throw MatrixException();
      m_elements[row][column] = value;
      return;
    };
    
    /**
     * \return the number of rows
     */
    size_t nRows()const{
      return m_nRows;
    }
    
    /**
     * \return the number of columns
     */
    size_t nColumns()const{
      return m_nColumns;
    }
    
  private:
    
    vector<vector<T> > m_elements;
    size_t m_nRows;
    size_t m_nColumns;
    
  };
  
  /**
   *  Write a vector to a stream
   */
  template<class T> inline std::ostream &writeVector(std::ostream &out, const std::vector<T> &row){
        
    for(typename std::vector<T>::const_iterator it = row.begin();
        it != row.end(); ++it){
//      out.width(25);
  //    out<<std::fixed<<*it;
    //  out.fill(' ');
      out<<*it<<" ";
    }
    out<<std::endl;
    return out;
  };
  
  
  /**
   *  Write a matrix to a stream
   */
  template<class T> inline std::ostream &operator << (std::ostream &out, const Matrix<T> &matrix){
    
    for(size_t r=0 ; r != matrix.nRows(); ++r){
      for(size_t c=0; c != matrix.nColumns(); ++c){
        out.width(12);
        out<<std::fixed<<matrix.element(r, c);
        out.fill(' ');
      }
      out<<std::endl;
    }
    
    return out;
  };
  
}
#endif
