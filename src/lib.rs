mod parallel;
mod tests;

use std::{fmt, ops::*};

use num::ToPrimitive;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Matrix<T> {
    matrix: Vec<Vec<T>>,
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
    pub fn new(matrix: Vec<Vec<T>>) -> Matrix<T> {
        let len_0 = &matrix[0].len();
        for i in &matrix {
            if i.len() != *len_0 {
                panic!("Error occured creating matrix: rows different size");
            }
        }
        Matrix { matrix }
    }

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

    pub fn mul(&self, other: Matrix<T>) -> Result<Self, String> {
        if self.matrix[0].len() == other.matrix.len() {
            // return Err("Size of matrices not maching".to_string());
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

    pub fn gaussian_elimination(&self) -> Matrix<f64> {
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

    pub fn determinant(&self) -> Result<f64, String> {
        if self.matrix.len() != self.matrix[0].len() {
            return Err("matrix not square".to_string());
        } else if self.matrix.len() == 2 {
            let mat = &self.matrix;
            let det_t = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
            let det = det_t.to_f64().ok_or("error occured on type conversion")?;
            return Ok(det);
        } else if self.matrix.len() == 1 {
            let det = self.matrix[0][0]
                .to_f64()
                .ok_or("error occured on type conversion")?;
            return Ok(det);
        } else if self.matrix.len() >= 3 {
            let temporary_matrix = Matrix::gaussian_elimination(&self.clone());
            let mut det = 1.0;
            for i in 0..temporary_matrix.matrix.len() {
                det *= temporary_matrix.matrix[i][i];
            }
            return Ok(det);
        }
        Err("something went wrong while calculating determinant".to_string())
    }

    pub fn vectorize(&self) -> Vec<T> {
        let mut result_vec = Vec::new();
        for i in 0..self.matrix.len() {
            for j in 0..self.matrix[0].len() {
                result_vec.push(self.matrix[j][i]);
            }
        }
        result_vec
    }

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

    pub fn inverse(&self) -> Result<Matrix<f64>, String> {
        if self.matrix.len() != self.matrix[0].len() || self.clone().determinant().unwrap() == 0.0 {
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
                if cloned[i][i] != 0.0 {
                    let scalar = cloned[j][i] / cloned[i][i];
                    for k in 0..len {
                        cloned[j][k] = cloned[j][k] - scalar * cloned[i][k];
                        result_matrix[j][k] = result_matrix[j][k] - scalar * result_matrix[i][k];
                    }
                }
            }
        }

        for i in 0..len {
            let scalar = 1.0 / cloned[i][i];
            for j in 0..len {
                result_matrix[i][j] = result_matrix[i][j] * scalar;
                cloned[i][j] = cloned[i][j] * scalar;
            }
        }

        for i in 0..len {
            for j in (i + 1)..len {
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
        $crate::Matrix::from_vectorized($e, $f)
    }
}
