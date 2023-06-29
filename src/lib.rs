//! matrix_lib is liblary giving easy way to do calcualtions on matrices
//! Currently implemented operations are:
//! sum, scalar multiplication, matrix multiplication, transposition, gaussian elimination, determinant, vectorization, devectorization and inversion
//!
//! example usage:
//! ```
//! use matrix_lib::*;
//!
//! let m1 = matrix![[1,2,3],[4,5,6],[7,8,9]];
//! assert_eq!(m1.scalar_mul(2), matrix![[2,4,6],[8,10,12],[14,16,18]]);
//!
//!
//! let m2 = matrix![[1,2,3],[4,5,6],[7,8,9]];
//! let m3 = matrix![[1,2,3],[4,5,6],[7, 8, 9]];
//! assert_eq!(m2.mul(m3).unwrap(), matrix![[30, 36, 42], [66, 81, 96], [102, 126, 150]]);
//! ```

pub mod functions;
pub mod parallel;
pub mod tests;
pub mod staticmatrix;

use std::{fmt, ops::*};

use functions::Size;
use num::ToPrimitive;

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<T> {
    matrix: Vec<Vec<T>>,
    __det_mul: f64,
}

impl<T> Matrix<T>
where
    T: std::ops::Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Neg<Output = T>
        + AddAssign
        + SubAssign
        + MulAssign
        + DivAssign
        + Clone
        + num::One
        + num::Zero
        + std::cmp::PartialEq
        + Copy
        + ToPrimitive,
{
    /// Creates a new [`Matrix<T>`].
    ///
    /// # Panics
    ///
    /// Panics if rows of matrix are different size.
    /// example:
    /// `
    /// matrix![[1, 2, 3], [4,5]]
    /// `
    ///
    pub fn new(matrix: Vec<Vec<T>>) -> Matrix<T> {
        let len_0 = &matrix[0].len();
        for i in &matrix {
            if i.len() != *len_0 {
                panic!("Error occured creating matrix: rows different size");
            }
        }
        Matrix {
            matrix,
            __det_mul: 1.0,
        }
    }

    /// Sums two matrices of same size.
    ///
    /// # Errors
    ///
    /// This function will return an error if matrices are different sizes.
    ///
    /// examples:
    /// ```
    /// use matrix_lib::*;
    /// let m1 = matrix![[1,2,3], [4,5,6]];
    /// let m2 = matrix![[1,2], [3,4], [5,6]];
    /// assert!(matches!(m1.sum(m2), Err(String)));
    /// ```
    /// will not sum successfully
    ///
    /// # Examples of correct usage
    /// ```
    /// use matrix_lib::*;
    /// let m1 = matrix![[1,2,3], [4,5,6]];
    /// let m2 = matrix![[1,2,3], [4,5,6]];
    /// assert_eq!(matrix![[2,4,6], [8,10,12]], m1.sum(m2).unwrap());
    /// ```
    ///
    pub fn sum(&self, other: Matrix<T>) -> Result<Self, String> {
        if self.matrix.len() != other.matrix.len() || self.matrix[0].len() != other.matrix[0].len()
        {
            return Err("Size of matrices not maching".to_string());
        }
        let mut result_matrix: Vec<Vec<T>> = Vec::new();
        for i in 0..self.matrix.len() {
            let mut row: Vec<T> = Vec::new();
            for j in 0..self.matrix[0].len() {
                row.push(self.matrix[i][j] + other.matrix[i][j]);
            }
            result_matrix.push(row);
        }

        Ok(Matrix::new(result_matrix))
    }

    /// Multiplies two matrices.
    ///
    /// # Errors
    /// ERRORS WIP IN THIS FUNCTION
    /// This function will return an error if ...
    pub fn mul(&self, other: Matrix<T>) -> Result<Self, String> {
        if self.matrix[0].len() != other.matrix.len() {
            return Err("Size of matrices not maching".to_string());
        }

        let mut result_matrix: Vec<Vec<T>> = Vec::new();

        for i in 0..self.matrix.len() {
            let mut result_row: Vec<T> = Vec::new();
            for j in 0..other.matrix[0].len() {
                let mut sum = T::zero();
                for k in 0..other.matrix.len() {
                    sum += self.matrix[i][k] * other.matrix[k][j];
                }
                result_row.push(sum);
            }
            result_matrix.push(result_row);
        }

        Ok(Matrix::new(result_matrix))
    }

    /// Multiplies matrice by number.
    pub fn scalar_mul(&self, scalar: T) -> Self {
        let mut result_matrix: Vec<Vec<T>> = Vec::new();
        for i in self.matrix.clone() {
            let mut row: Vec<T> = Vec::new();
            for j in i {
                row.push(j * scalar);
            }
            result_matrix.push(row);
        }
        Matrix::new(result_matrix)
    }

    /// Returns the transpose of this [`Matrix<T>`].
    pub fn transpose(&self) -> Matrix<T> {
        let mut result_matrix: Vec<Vec<T>> = Vec::new();
        for i in 0..self.matrix[0].len() {
            let mut result_row: Vec<T> = Vec::new();
            for j in 0..self.matrix.len() {
                result_row.push(self.matrix[j][i]);
            }
            result_matrix.push(result_row);
        }
        Matrix::new(result_matrix)
    }

    /// Returns the gaussian elimination of this [`Matrix<T>`].
    ///
    /// # Panics
    ///
    /// Panics if type T is not convertable to f64 as f64 is needed to calculate this correctly.
    pub fn gaussian_elimination(&mut self) -> Matrix<f64> {
        let mut cloned = Vec::new();
        for i in 0..self.matrix.len() {
            let mut row = Vec::new();
            for j in 0..self.matrix[0].len() {
                row.push(self.matrix[i][j].to_f64().unwrap());
            }
            cloned.push(row);
        }

        for i in 0..cloned.len() {
            for j in (i + 1)..cloned.len() {
                if cloned[i][i] == 0.0 {
                    for r in (i + 1)..cloned.len() {
                        if cloned[r][i] != 0.0 {
                            let temp = cloned[j].clone();
                            cloned[j] = cloned[i].clone();
                            cloned[i] = temp;
                            if r - i % 2 != 0 {
                                self.__det_mul = self.__det_mul * -1.0;
                            }
                            break;
                        }
                    }
                }
                if cloned[i][i] != 0.0 {
                    let scalar = cloned[j][i] / cloned[i][i];
                    for k in 0..cloned[0].len() {
                        cloned[j][k] = cloned[j][k] - scalar * cloned[i][k];
                    }
                }
            }
        }
        Matrix::new(cloned)
    }

    /// Returns the determinant of this [`Matrix<T>`].
    ///
    /// # Errors
    ///
    /// This function will return an error if matrix is not square.
    pub fn determinant(&mut self) -> Result<f64, String> {
        if self.matrix.len() != self.matrix[0].len() {
            return Err("matrix not square".to_string());
        } else if self.matrix.len() == 2 {
            let mat = &self.matrix;
            let det_t = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
            let det = det_t.to_f64().ok_or("error occured on type conversion")?;
            return Ok(det * self.__det_mul);
        } else if self.matrix.len() == 1 {
            let det = self.matrix[0][0]
                .to_f64()
                .ok_or("error occured on type conversion")?;
            return Ok(det * self.__det_mul);
        } else if self.matrix.len() >= 3 {
            let temporary_matrix = Matrix::gaussian_elimination(self);
            let mut det = 1.0;
            for i in 0..temporary_matrix.matrix.len() {
                det *= temporary_matrix.matrix[i][i];
            }
            return Ok(det * self.__det_mul);
        }
        Err("something went wrong while calculating determinant".to_string())
    }

    /// Returns the vectorize of this [`Matrix<T>`].
    pub fn vectorize(&self) -> Vec<T> {
        let mut result_vec = Vec::new();
        for i in 0..self.matrix.len() {
            for j in 0..self.matrix[0].len() {
                result_vec.push(self.matrix[j][i]);
            }
        }
        result_vec
    }

    /// Generates [`Matrix<T>`] from provided vector.
    /// # Example
    /// ```
    /// use matrix_lib::*;
    /// let v = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
    /// assert_eq!(matrix![1,4,7;2,5,8;3,6,9], Matrix::from_vectorized(v, 3));
    /// ```
    pub fn from_vectorized(v: Vec<T>, columns: usize) -> Matrix<T> {
        let mut result_matrix: Vec<Vec<T>> = Vec::new();
        for i in 0..columns {
            let mut row: Vec<T> = Vec::new();
            for j in 0..v.len() / columns {
                row.push(v[j * columns + i].clone());
            }
            result_matrix.push(row);
        }
        Matrix::new(result_matrix)
    }

    /// Returns the inverse of this [`Matrix<T>`].
    ///
    /// # Panics
    ///
    /// Panics if type T is not convertable to f64 as f64 is needed to calculate this correctly.
    ///
    /// # Errors
    ///
    /// This function will return an error if matrix is not suqare or determinant from given matrix is equal to 0.
    pub fn inverse(&self) -> Result<Matrix<f64>, String> {
        if self.matrix.len() != self.matrix[0].len() || self.clone().determinant()? == 0.0 {
            return Err(String::from("matrix not square or det = 0"));
        }

        let len = self.matrix.len();

        let mut result_matrix = Vec::new();
        let mut cloned = Vec::new();

        for i in 0..len {
            let mut row = Vec::new();
            let mut result_row = Vec::new();
            for j in 0..len {
                row.push(self.matrix[i][j].to_f64().unwrap());
                result_row.push(0.0);
            }
            cloned.push(row);
            result_matrix.push(result_row);
        }

        for i in 0..len {
            result_matrix[i][i] = 1.0;
        }

        for i in 0..len {
            for j in (i + 1)..len {
                if cloned[i][i] == 0.0 {
                    for r in (i + 1)..len {
                        if cloned[r][i] != 0.0 {
                            let temp = cloned[j].clone();
                            cloned[j] = cloned[i].clone();
                            cloned[i] = temp;
                            break;
                        }
                    }
                }
                if cloned[i][i] != 0.0 {
                    let scalar = cloned[j][i] / cloned[i][i];
                    for k in 0..len {
                        cloned[j][k] = cloned[j][k] - scalar * cloned[i][k];
                        result_matrix[j][k] = result_matrix[j][k] - scalar * result_matrix[i][k];
                    }
                }
            }
            let scalar = 1.0 / cloned[i][i];
            for j in 0..len {
                result_matrix[i][j] = result_matrix[i][j] * scalar;
                cloned[i][j] = cloned[i][j] * scalar;
            }
        }

        for i in 0..len {
            for j in (i + 1)..len {
                if cloned[i][i] == 0.0 {
                    for r in (i + 1)..len {
                        if cloned[r][i] != 0.0 {
                            let temp = cloned[len - 1 - j].clone();
                            cloned[len - 1 - j] = cloned[len - 1 - i].clone();
                            cloned[len - 1 - i] = temp;
                            break;
                        }
                    }
                }
                if cloned[i][i] != 0.0 {
                    let scalar =
                        cloned[len - 1 - j][len - 1 - i] / cloned[len - 1 - i][len - 1 - i];
                    for k in 0..len {
                        cloned[len - 1 - j][k] =
                            cloned[len - 1 - j][k] - scalar * cloned[len - 1 - i][k];
                        result_matrix[len - 1 - j][k] =
                            result_matrix[len - 1 - j][k] - scalar * result_matrix[len - 1 - i][k];
                    }
                }
            }
        }
        Ok(Matrix::new(result_matrix))
    }

    /// Returns the get size of this [`Matrix<T>`].
    ///
    /// # Example
    /// ```
    /// use matrix_lib::*;
    /// use matrix_lib::functions::Size;
    /// let m = matrix![1,2,3;4,5,6;7,8,9;10,11,12];
    /// assert_eq!(m.get_size(), Size{x: 3, y: 4});
    /// ```
    pub fn get_size(&self) -> Size {
        return Size {
            x: self.matrix[0].len(),
            y: self.matrix.len(),
        };
    }
}

impl<T> fmt::Display for Matrix<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.matrix.len() {
            for j in 0..self.matrix[i].len() {
                write!(f, "{} ", self.matrix[i][j])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

/// Hels with matrix creation
///
/// examples:
/// ```
/// use matrix_lib::*;
/// let m1 = matrix![[1,2,3], [4,5,6], [7,8,9]]; // each array provided in this macro is separate row
///
/// let m2 = matrix![1,2,3;4,5,6;7,8,9]; // rows are separated by semicolon `;`
///
/// let m3 = matrix![1 2 3; 4 5 6; 7 8 9]; // rows are separated by semicolon but with elements separated by space
///
/// let m4 = matrix!([1,2,3,4,5,6,7,8,9] => 3); // splits vectorized matrix into columns of matrix with 3 elements by each
/// assert_eq!(m1, m2);
/// assert_eq!(m2, m3);
/// assert_eq!(m4, matrix![[1, 4, 7], [2, 5, 8], [3, 6, 9]]);
/// ```
/// 
#[macro_export]
macro_rules! matrix {
    [$($e:expr),*] => {
        $crate::Matrix::new(vec![$($e.to_vec()),*])
    };

    [ $($( $f:expr ),*);*] => {
        $crate::Matrix::new(vec![$(vec![$($f),*],)*])
    };

    [$($( $f:expr) *);*] => {
        $crate::Matrix::new(vec![$(vec![$($f),*],)*])
    };

    ($e:expr => $f:expr) => {
        $crate::Matrix::from_vectorized($e.to_vec(), $f)
    }
}
