//! Provides basic parallelism for calcualtions on matrices
//! 
//! Available operations are:
//! par_scalar_mul -> multiplication of matrix by number
//! par_mul -> multiplication of matrix by matrix 
//! par_sum -> sum of matrices
 
 

use std::ops::*;

use num::ToPrimitive;

use crate::Matrix;

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
        + ToPrimitive
        + std::marker::Send
        + 'static,
{
    pub async fn par_scalar_mul(self, scalar: T) -> Matrix<T> {
        let futures = {
            let mut handles = vec![];
            for i in 0..self.matrix.len() {
                let t_matrix = self.matrix[i].clone();
                handles.push(tokio::spawn(async move {
                    let mut result_row = Vec::new();
                    for j in 0..t_matrix.len() {
                        result_row.push(t_matrix[j] * scalar.clone());
                    }
                    result_row
                }));
            }
            handles
        };

        let mut result = Vec::new();
        for f in futures {
            result.push(f.await.unwrap());
        }

        Matrix::new(result)
    }

    pub async fn par_mul(self, other: Matrix<T>) -> Result<Matrix<T>, String> {
        let mut result_matrix: Vec<Vec<T>> = Vec::new();
        let futures = {
            let mut handles = Vec::new();
            for i in 0..self.matrix.len() {
                let t_matrix = self.matrix.clone();
                let t_other = other.matrix.clone();
                handles.push(tokio::spawn(async move {
                    let mut result_row: Vec<T> = Vec::new();
                    for j in 0..t_other[0].len() {
                        let mut sum = T::zero();
                        for k in 0..t_other.len() {
                            sum += t_matrix[i][k] * t_other[k][j];
                        }
                        result_row.push(sum);
                    }
                    result_row
                }));
            }
            handles
        };
        for f in futures {
            result_matrix.push(f.await.unwrap());
        }

        Ok(Matrix::new(result_matrix))
    }

    pub async fn par_sum(self, other: Matrix<T>) -> Result<Matrix<T>, String> {
        if self.matrix.len() != other.matrix.len() && self.matrix[0].len() != other.matrix[0].len()
        {
            return Err("matrix dimensions do not match".to_string());
        }

        let futures = {
            let mut handles = Vec::new();
            for i in 0..self.matrix.len() {
                let t_matrix = self.matrix.clone();
                let t_other = other.matrix.clone();
                handles.push(tokio::spawn(async move {
                    let mut result_row: Vec<T> = Vec::new();
                    for j in 0..t_other[0].len() {
                        result_row.push(t_matrix[i][j] + t_other[i][j]);
                    }
                    result_row
                }));
            }
            handles
        };
        let mut result_matrix = Vec::new();
        for f in futures {
            result_matrix.push(f.await.unwrap());
        }

        Ok(Matrix::new(result_matrix))
    }
}
