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
#[cfg(target_feature = "parallel")]
pub mod parallel;
pub mod staticmatrix;
pub mod tests;

use std::{fmt, ops::*};

use functions::Size;
use num::ToPrimitive;

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<T> {
    cols: usize,
    rows: usize,
    matrix: Vec<T>,
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
        let rows = matrix.len();
        let cols = *(&matrix[0].len());
        for i in &matrix {
            if i.len() != cols {
                panic!("Error occured creating matrix: rows different size");
            }
        }

        Matrix {
            matrix: matrix.into_iter().flatten().collect(),
            __det_mul: 1.0,
            cols,
            rows,
        }
    }

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
    pub fn new_from_vec(matrix: Vec<T>, cols: usize, rows: usize) -> Matrix<T> {
        Matrix {
            matrix,
            __det_mul: 1.0,
            cols,
            rows,
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
    pub fn sum(&self, other: &Matrix<T>) -> Result<Self, String> {
        if self.cols != other.cols || self.rows != other.rows {
            return Err("Size of matrices not maching".to_string());
        }

        let mut data = Vec::new();

        for i in 0..self.cols * self.rows {
            data.push(self.matrix[i] + other.matrix[i])
        }

        Ok(Matrix::new_from_vec(data, self.cols, self.rows))
    }

    /// Multiplies two matrices.
    ///
    /// # Errors
    /// ERRORS WIP IN THIS FUNCTION
    /// This function will return an error if ...
    pub fn mul(&self, other: Matrix<T>) -> Result<Self, String> {
        if self.cols != other.rows {
            return Err("Size of matrices not maching".to_string());
        }

        let mut data = vec![T::zero(); self.rows * other.cols];

        for i in 0..self.rows {
            for j in 0..other.cols {
                let mut sum = T::zero();
                for k in 0..self.cols {
                    sum += self.matrix[i * self.cols + k] * other.matrix[k * other.cols + j];
                }

                data[i * other.cols + j] = sum;
            }
        }

        Ok(Matrix::new_from_vec(data, self.rows, other.cols))
    }

    /// Multiplies matrice by number.
    pub fn scalar_mul(&self, scalar: T) -> Self {
        Matrix::new_from_vec(
            self.matrix.iter().map(|x| *x * scalar).collect(),
            self.cols,
            self.rows,
        )
    }

    /// Returns the transpose of this [`Matrix<T>`].
    pub fn transpose(&self) -> Matrix<T> {
        let mut data = vec![T::zero(); self.cols * self.rows];

        for i in 0..self.rows {
            for j in 0..self.cols {
                data[j * self.rows + i] = self.matrix[i * self.cols + j];
            }
        }

        Matrix::new_from_vec(data, self.rows, self.cols)
    }

    /// Returns the gaussian elimination of this [`Matrix<T>`].
    ///
    /// # Panics
    ///
    /// Panics if type T is not convertable to f64 as f64 is needed to calculate this correctly.
    pub fn gaussian_elimination(&mut self) -> Matrix<f64> {
        let mut cloned: Vec<f64> = self.matrix.iter().map(|x| x.to_f64().unwrap()).collect();

        for i in 0..self.rows {
            if cloned[i * self.cols + i] == 0.0 {
                for r in (i + 1)..self.rows {
                    if cloned[r * self.cols + i] != 0.0 {
                        for col in 0..self.cols {
                            let temp = cloned[i * self.cols + col];
                            cloned[i * self.cols + col] = cloned[r * self.cols + col];
                            cloned[r * self.cols + col] = temp;
                        }
                        if (r - i) % 2 != 0 {
                            self.__det_mul *= -1.0;
                        }
                        break;
                    }
                }
            }

            if cloned[i * self.cols + i] != 0.0 {
                for j in (i + 1)..self.rows {
                    let scalar = cloned[j * self.cols + i] / cloned[i * self.cols + i];
                    for k in 0..self.cols {
                        cloned[j * self.cols + k] -= scalar * cloned[i * self.cols + k];
                    }
                }
            }
        }

        Matrix::new_from_vec(cloned, self.cols, self.rows)
    }

    /// Returns the determinant of this [`Matrix<T>`].
    ///
    /// # Errors
    ///
    /// This function will return an error if matrix is not square.
    pub fn determinant(&mut self) -> Result<f64, String> {
        if self.cols != self.rows {
            return Err("matrix not square".to_string());
        } else if self.cols == 2 {
            let mat = &self.matrix;
            let det_t = mat[0] * mat[3] - mat[2] * mat[1];
            let det = det_t.to_f64().ok_or("error occured on type conversion")?;
            return Ok(det * self.__det_mul);
        } else if self.cols == 1 {
            let det = self.matrix[0]
                .to_f64()
                .ok_or("error occured on type conversion")?;
            return Ok(det * self.__det_mul);
        } else if self.matrix.len() >= 3 {
            let temporary_matrix = Matrix::gaussian_elimination(self);
            let mut det = 1.0;
            for i in 0..temporary_matrix.cols {
                det *= temporary_matrix.matrix[i * self.cols + i];
            }
            return Ok(det * self.__det_mul);
        }
        Err("something went wrong while calculating determinant".to_string())
    }

    /// Returns the vectorize of this [`Matrix<T>`].
    pub fn vectorize(&self) -> Vec<T> {
        self.matrix.clone()
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
        if self.cols != self.rows || self.clone().determinant()? == 0.0 {
            return Err(String::from("matrix not square or det = 0"));
        }

        let len = self.rows;
        let n = len;
        
        let mut cloned: Vec<f64> = self.matrix.iter().map(|x| x.to_f64().unwrap()).collect();
        let mut result_matrix = vec![0.0; self.rows * self.rows];
        for i in 0..self.rows {
            result_matrix[i * self.rows + i] = 1.0;
        }

        // for i in 0..len {
        //     for j in (i + 1)..len {
        //         if cloned[i * self.rows + i] == 0.0f64 {
        //             for r in (i + 1)..len {
        //                 if cloned[r * self.rows + i] != 0.0 {
        //                     let temp = cloned[j].clone();
        //                     cloned[j] = cloned[i].clone();
        //                     cloned[i] = temp;
        //                     break;
        //                 }
        //             }
        //         }
        //         if cloned[i * self.rows + i] != 0.0 {
        //             let scalar = cloned[j*self.rows+i] / cloned[i*self.rows+i];
        //             for k in 0..len {
        //                 cloned[j * self.rows + k] = cloned[j * self.rows + k] - scalar * cloned[i * self.rows +k];
        //                 result_matrix[j*self.rows+k] = result_matrix[j*self.rows+k] - scalar * result_matrix[i*self.rows+k];
        //             }
        //         }
        //     }
        //     let scalar = 1.0 / cloned[i*self.rows+i];
        //     for j in 0..len {
        //         result_matrix[i*self.rows+j] = result_matrix[i*self.rows+j] * scalar;
        //         cloned[i*self.rows+j] = cloned[i*self.rows+j] * scalar;
        //     }
        // }

        // for i in 0..len {
        //     for j in (i + 1)..len {
        //         if cloned[i*self.rows+i] == 0.0 {
        //             for r in (i + 1)..len {
        //                 if cloned[r*self.rows+i] != 0.0 {
        //                     let temp = cloned[len - 1 - j].clone();
        //                     cloned[len - 1 - j] = cloned[len - 1 - i].clone();
        //                     cloned[len - 1 - i] = temp;
        //                     break;
        //                 }
        //             }
        //         }
        //         if cloned[i*self.rows+i] != 0.0 {
        //             let scalar =
        //                 cloned[(len - 1 - j)*self.rows+(len - 1 - i)] / cloned[(len - 1 - i)*self.rows+(len - 1 - i)];
        //             for k in 0..len {
        //                 cloned[(len - 1 - j)*self.rows+k] =
        //                     cloned[(len - 1 - j)*self.rows+k] - scalar * cloned[(len - 1 - i)*self.rows+k];
        //                 result_matrix[(len - 1 - j)*self.rows+k] =
        //                     result_matrix[(len - 1 - j)*self.rows+k] - scalar * result_matrix[(len - 1 - i)*self.rows+k];
        //             }
        //         }
        //     }
        // }

           // Perform Gaussian elimination on the augmented matrix
           for i in 0..n {
            // Pivoting: Ensure that the diagonal element is non-zero
            let pivot_index = i * n + i;
            if cloned[pivot_index].abs() < 1e-10 {
                let mut swapped = false;
                for row in (i + 1)..n {
                    if cloned[row * n + i].abs() > 1e-10 {
                        // Swap rows in both original and identity matrix
                        for col in 0..n {
                            cloned.swap(i * n + col, row * n + col);
                            result_matrix.swap(i * n + col, row * n + col);
                        }
                        swapped = true;
                        break;
                    }
                }
                if !swapped {
                    panic!("Matrix is singular and cannot be inverted.");
                }
            }

            // Normalize the pivot row
            let pivot_value = cloned[pivot_index];
            for col in 0..n {
                cloned[i * n + col] /= pivot_value;
                result_matrix[i * n + col] /= pivot_value;
            }

            // Eliminate the current column in other rows
            for row in 0..n {
                if row != i {
                    let factor = cloned[row * n + i];
                    for col in 0..n {
                        cloned[row * n + col] -= factor * cloned[i * n + col];
                        result_matrix[row * n + col] -= factor * result_matrix[i * n + col];
                    }
                }
            }
        }
        Ok(Matrix::new_from_vec(result_matrix, self.rows, self.rows))
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
            x: self.rows,
            y: self.cols,
        };
    }
}

impl<T> fmt::Display for Matrix<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.rows {
            for j in 0..self.cols {
                write!(f, "{} ", self.matrix[i * self.rows + j])?;
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
