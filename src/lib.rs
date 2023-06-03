use std::{fmt, ops::*};

use num::ToPrimitive;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Matrix<T>
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
        + fmt::Display
        + Clone
        + num::One
        + num::Zero
        + std::cmp::PartialEq
        + Copy
        + ToPrimitive,
{
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
        + fmt::Display
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
                panic!("Error occured creating matrix: matrices different size");
            }
        }
        Matrix { matrix }
    }

    pub fn sum(self, other: Matrix<T>) -> Result<Self, String> {
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

        Ok(Matrix {
            matrix: result_matrix,
        })
    }

    pub fn mul(self, other: Matrix<T>) -> Result<Self, String> {
        if self.matrix.len() == other.matrix[0].len() {
            //return Err("Size of matrices not maching".to_string());
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

        Ok(Matrix {
            matrix: result_matrix,
        })
    }

    pub fn scalar_mul(self, scalar: T) -> Self {
        let mut result_matrix: Vec<Vec<T>> = Vec::new();
        for i in self.matrix {
            let mut row: Vec<T> = Vec::new();
            for j in i {
                row.push(j * scalar);
            }
            result_matrix.push(row);
        }
        Matrix {
            matrix: result_matrix,
        }
    }

    pub fn transpose(self) -> Matrix<T> {
        let mut result_matrix: Vec<Vec<T>> = Vec::new();
        for i in 0..self.matrix[0].len() {
            let mut result_row: Vec<T> = Vec::new();
            for j in 0..self.matrix.len() {
                result_row.push(self.matrix[j][i]);
            }
            result_matrix.push(result_row);
        }
        Matrix {
            matrix: result_matrix,
        }
    }

    pub fn gaussian_elimination(self) -> Matrix<f64> {
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
        Matrix { matrix: cloned }
    }

    // laplace expansion is not needed any more
    // pub fn laplace_expansion_determinant(m: Matrix) -> T {
    //     if m.matrix.len() == 1 {
    //         return m.matrix[0][0];
    //     }

    //     let mut cloned = m.clone();
    //     // let cloned = cloned.gaussian_elimination();
    //     println!("{}", cloned);
    //     let mut result = 0;
    //     for i in 0..cloned.matrix.len() {
    //         let mut temp_matrix: Vec<Vec<T>> = Vec::new();
    //         if cloned[i][0] != 0 {
    //             println!("{:?}", cloned[i]);
    //             for j in 0..cloned.matrix[0].len() {
    //                 let mut row: Vec<T> = Vec::new();
    //                 if j != i {
    //                     for k in 1..cloned.matrix.len() {
    //                         row.push(cloned[j][k]);
    //                     }
    //                     temp_matrix.push(row);
    //                 }
    //             }
    //             println!("{:?}", temp_matrix);
    //             println!("{}", Matrix::determinant(Matrix::new(temp_matrix.clone())).unwrap());
    //             result += cloned[i][0] * pow(cloned[i][0], i+1) * Matrix::determinant(Matrix::new(temp_matrix)).unwrap();
    //         }
    //     }
    //     result
    // }

    pub fn determinant(self) -> Result<f64, String> {
        if self.matrix.len() != self.matrix[0].len() {
            return Err("matrix not square".to_string());
        }

        let len = self.matrix.len();

        if self.matrix.len() > 3 {
            let temporary_matrix = Matrix::gaussian_elimination(self.clone());
            let mut det = 1.0;
            for i in 0..temporary_matrix.matrix.len() {
                det *= temporary_matrix.matrix[i][i];
            }
            return Ok(det);
        }

        let mut temporary_matrix: Vec<Vec<T>> = self.matrix.clone();
        for m in &self.matrix {
            temporary_matrix.push(m.clone());
        }

        let mut sum_part = T::zero();
        let mut sub_part = T::zero();
        for j in 0..len {
            let mut sum_temp = T::one();
            let mut sub_temp = T::one();
            for i in 0..len {
                sum_temp *= temporary_matrix[i + j][i];
                sub_temp *= temporary_matrix[len - j + i][len - i - 1];
            }
            sum_part += sum_temp;
            sub_part += sub_temp;
        }
        let det = sum_part - sub_part;

        Ok(det.to_f64().unwrap())
    }
}

impl<T> fmt::Display for Matrix<T>
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
        + fmt::Display
        + Clone
        + num::One
        + num::Zero
        + std::cmp::PartialEq
        + Copy
        + ToPrimitive,
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
        Matrix::new(vec![$($e),*])
    }
}

#[cfg(test)]
mod tests {
    use crate::Matrix;

    #[test]
    fn sum() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        let m2 = matrix![vec![3, 2, 1], vec![6, 5, 4]];
        assert_eq!(
            matrix![vec![4, 4, 4], vec![10, 10, 10]],
            m1.sum(m2).unwrap()
        )
    }

    #[test]
    fn sum_fail() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        let m2 = matrix![vec![3, 2], vec![6, 5]];
        assert!(m1.sum(m2).is_err());
    }

    #[test]
    fn mul() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        let m2 = matrix![vec![3, 2], vec![6, 5], vec![6, 5]];
        assert_eq!(matrix![vec![33, 27], vec![78, 63]], m1.mul(m2).unwrap())
    }
    #[test]
    fn scalar_mul() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        assert_eq!(
            matrix![vec![10, 20, 30], vec![40, 50, 60]],
            m1.scalar_mul(10)
        )
    }

    #[test]
    fn transpose() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        assert_eq!(matrix![vec![1, 4], vec![2, 5], vec![3, 6]], m1.transpose())
    }

    #[test]
    fn determinant() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]];
        assert_eq!(0 as f64, m1.determinant().unwrap());
        let m2 = matrix![vec![5, 4, 1], vec![9, 2, 2], vec![4, 4, 0]];
        assert_eq!(20 as f64, m2.determinant().unwrap());
    }

    #[test]
    fn gaussian_elimination() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]];
        assert_eq!(
            matrix![
                vec![1.0, 2.0, 3.0],
                vec![0.0, -3.0, -6.0],
                vec![0.0, 0.0, 0.0]
            ],
            m1.gaussian_elimination()
        );
    }

    #[test]
    fn gaussian_elimination_with_division_by_zero() {
        let m1 = matrix![vec![1, -2, 3], vec![-5, 10, 14], vec![3, -6, -9]];
        assert_eq!(
            matrix![
                vec![1.0, -2.0, 3.0],
                vec![0.0, 0.0, 29.0],
                vec![0.0, 0.0, -18.0]
            ],
            m1.gaussian_elimination()
        )
    }

    #[test]
    fn determinant_cmp_gaussian_elimination_determinant() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]];
        let m2 = m1.clone();
        assert_eq!(
            m1.determinant().unwrap() as f64,
            m2.gaussian_elimination().determinant().unwrap()
        );
    }
}
