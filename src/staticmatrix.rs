use std::{
    fmt,
    ops::{AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use num::ToPrimitive;

use crate::functions::Size;

#[derive(Debug, Copy, Clone)]
pub struct StaticMatrix<T, const N: usize, const M: usize> {
    matrix: [[T; M]; N],
    __det_mul: f64,
}

impl<T, const N: usize, const M: usize> StaticMatrix<T, N, M>
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
        + std::fmt::Debug,
{
    pub fn init() -> Self {
        StaticMatrix {
            matrix: [[T::zero(); M]; N],
            __det_mul: 1.0,
        }
    }
    pub fn new(matrix: [[T; M]; N]) -> Self {
        StaticMatrix {
            matrix,
            __det_mul: 1.0,
        }
    }
    pub fn new_from_vec(vec: Vec<Vec<T>>) -> Self {
        let mut matrix: [[T; M]; N] = [[T::zero(); M]; N];

        let len_0 = &vec[0].len();

        for i in 0..vec.len() {
            if vec[i].len() != *len_0 {
                panic!("Error occured creating matrix: rows different size");
            }
            matrix[i] = vec[i]
                .clone()
                .try_into()
                .expect("Error occured on converting `Vec<T>` to `[T; M]`");
        }
        StaticMatrix {
            matrix,
            __det_mul: 1.0,
        }
    }

    pub fn sum(&self, other: StaticMatrix<T, N, M>) -> Self {
        let mut result_matrix = [[T::zero(); M]; N];
        for i in 0..N {
            let mut row = [T::zero(); M];
            for j in 0..M {
                row[j] = self.matrix[i][j] + other.matrix[i][j];
            }
            result_matrix[i] = row;
        }

        StaticMatrix::new(result_matrix)
    }

    pub fn mul<const K: usize>(&self, other: StaticMatrix<T, M, K>) -> Self {
        let mut result_matrix = [[T::zero(); M]; N];
        for i in 0..self.matrix.len() {
            let mut result_row = [T::zero(); M];
            for j in 0..other.matrix[0].len() {
                let mut sum = T::zero();
                for k in 0..other.matrix.len() {
                    sum += self.matrix[i][k] * other.matrix[k][j];
                }
                result_row[j] = sum;
            }
            result_matrix[i] = result_row;
        }

        StaticMatrix::new(result_matrix)
    }

    pub fn scalar_mul(&self, scalar: T) -> Self {
        let mut result_matrix = [[T::zero(); M]; N];
        for i in 0..N {
            let mut row = [T::zero(); M];
            for j in 0..M {
                row[j] *= scalar;
            }
            result_matrix[i] = row;
        }
        StaticMatrix::new(result_matrix)
    }

    pub fn transpose(&self) -> StaticMatrix<T, M, N> {
        let mut result_matrix = [[T::zero(); N]; M];
        for i in 0..self.matrix[0].len() {
            let mut row = [T::zero(); N];
            for j in 0..self.matrix.len() {
                row[j] = self.matrix[j][i];
            }
            result_matrix[i] = row;
        }
        StaticMatrix::new(result_matrix)
    }

    pub fn gaussian_elimination(&mut self) -> StaticMatrix<f64, N, M> {
        let mut cloned: [[f64; M]; N] = [[0.0; M]; N];

        for i in 0..self.matrix.len() {
            let mut row = [0.0; M];
            for j in 0..self.matrix[0].len() {
                row[j] = self.matrix[i][j].to_f64().unwrap();
            }
            cloned[i] = row;
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

        StaticMatrix::new(cloned)
    }

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
            let temporary_matrix = StaticMatrix::gaussian_elimination(self);
            let mut det = 1.0;
            for i in 0..temporary_matrix.matrix.len() {
                det *= temporary_matrix.matrix[i][i];
            }
            return Ok(det * self.__det_mul);
        }
        Err("something went wrong while calculating determinant".to_string())
    }

    pub fn inverse(&self) -> Result<StaticMatrix<f64, N, M>, String> {
        if self.matrix.len() != self.matrix[0].len() || self.clone().determinant()? == 0.0 {
            return Err(String::from("matrix not square or det = 0"));
        }

        let len = self.matrix.len();

        let mut result_matrix = [[0.0; M]; N];
        let mut cloned = [[0.0; M]; N];

        for i in 0..self.matrix.len() {
            let mut row = [0.0; M];
            for j in 0..self.matrix[0].len() {
                row[j] = self.matrix[i][j].to_f64().unwrap();
            }
            cloned[i] = row;
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

        Ok(StaticMatrix::new(result_matrix))
    }

    pub fn get_size(&self) -> Size {
        return Size {
            x: self.matrix[0].len(),
            y: self.matrix.len(),
        };
    }
}

impl<T, const N: usize, const M: usize> fmt::Display for StaticMatrix<T, N, M>
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
macro_rules! static_matrix {
    [$($e:expr),*] => {
        $crate::staticmatrix::StaticMatrix::new([$($e),*])
    };

    [ $($( $f:expr ),*);*] => {
        $crate::staticmatrix::StaticMatrix::new([$([$($f),*],)*])
    };

    [$($( $f:expr) *);*] => {
        $crate::staticmatrix::StaticMatrix::new([$([$($f),*],)*])
    };
}
