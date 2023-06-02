use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Matrix {
    matrix: Vec<Vec<i32>>,
}

impl Matrix {
    pub fn new(matrix: Vec<Vec<i32>>) -> Self {
        let len_0 = &matrix[0].len();
        for i in &matrix {
            if i.len() != *len_0 {
                panic!("Error occured creating matrix: matrices different size");
            }
        }
        Matrix { matrix }
    }

    pub fn sum(self, other: &Matrix) -> Result<Self, String> {
        if self.matrix.len() != other.matrix.len() || self.matrix[0].len() != other.matrix[0].len()
        {
            return Err("Size of matrices not maching".to_string());
        }
        let mut result_matrix: Vec<Vec<i32>> = Vec::new();
        for i in 0..self.matrix.len() {
            let mut row: Vec<i32> = Vec::new();
            for j in 0..self.matrix[0].len() {
                row.push(self.matrix[i][j] + other.matrix[i][j]);
            }
            result_matrix.push(row);
        }

        Ok(Matrix {
            matrix: result_matrix,
        })
    }

    pub fn mul(self, other: &Matrix) -> Result<Self, String> {
        if self.matrix.len() == other.matrix[0].len()
        {
            //return Err("Size of matrices not maching".to_string());
        }

        let mut result_matrix: Vec<Vec<i32>> = Vec::new();

        for i in 0..self.matrix.len() {
            let mut result_row: Vec<i32> = Vec::new();
            for j in 0..other.matrix[0].len() {
                let mut sum = 0;
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

    pub fn scalar_mul(self, scalar: i32) -> Self {
        let mut result_matrix: Vec<Vec<i32>> = Vec::new();
        for i in self.matrix {
            let mut row: Vec<i32> = Vec::new();
            for j in i {
                row.push(j * scalar);
            }
            result_matrix.push(row);
        }
        Matrix {
            matrix: result_matrix,
        }
    }

    pub fn transpose(self) -> Matrix {
        let mut result_matrix: Vec<Vec<i32>> = Vec::new();
        for i in 0..self.matrix[0].len() {
            let mut result_row: Vec<i32> = Vec::new();
            for j in 0..self.matrix.len() {
                result_row.push(self.matrix[j][i]);
            }
            result_matrix.push(result_row);
        }
        Matrix {
            matrix: result_matrix,
        }
    }

    pub fn gaussian_elimination(self) -> Self {
        if self.matrix.len() == 1 {
            return self;
        }
        let mut cloned = self.clone();

        for i in 0..cloned.matrix.len() {
            for j in (i + 1)..cloned.matrix.len() {
                if cloned.matrix[i][i] != 0 {
                    let scalar = cloned.matrix[j][i] / cloned.matrix[i][i];
                    for k in 0..cloned.matrix[0].len() {
                        cloned.matrix[j][k] -= scalar * cloned.matrix[i][k];
                    }
                }
            }
        }
        cloned
    }

    pub fn determinant(self) -> Result<i32, String> {
        if self.matrix.len() != self.matrix[0].len() {
            return Err("matrix not square".to_string());
        }

        let mut temporary_matrix: Vec<Vec<i32>> = self.matrix.clone();
        for m in &self.matrix {
            temporary_matrix.push(m.clone());
        }

        let len = self.matrix.len();

        if self.matrix.len() > 3 {
            todo!()
        }

        let mut sum_part = 0;
        let mut sub_part = 0;
        for j in 0..len {
            let mut sum_temp = 1;
            let mut sub_temp = 1;
            for i in 0..len {
                sum_temp *= temporary_matrix[i + j][i];
                sub_temp *= temporary_matrix[len - 1 - i + j][len - 1 - i];
            }
            sum_part += sum_temp;
            sub_part += sub_temp;
        }
        let det = sum_part - sub_part;

        Ok(det)
    }
}

impl fmt::Display for Matrix {
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
            m1.sum(&m2).unwrap()
        )
    }

    #[test]
    fn sum_fail() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        let m2 = matrix![vec![3, 2], vec![6, 5]];
        assert!(m1.sum(&m2).is_err());
    }

    #[test]
    fn mul() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        let m2 = matrix![vec![3, 2], vec![6, 5], vec![6, 5]];
        assert_eq!(matrix![vec![33, 27], vec![78, 63]], m1.mul(&m2).unwrap())
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
        assert_eq!(0, m1.determinant().unwrap())
    }

    #[test]
    fn gaussian_elimination() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]];
        assert_eq!(
            matrix![vec![1, 2, 3], vec![0, -3, -6], vec![0, 0, 0]],
            m1.gaussian_elimination()
        );
    }

    #[test]
    fn gaussian_elimination_with_division_by_zero() {
        let m1 = matrix![vec![1, -2, 3], vec![-5, 10, 14], vec![3, -6, -9]];
        assert_eq!(
            matrix![vec![1, -2, 3], vec![0, 0, 29], vec![0, 0, -18]],
            m1.gaussian_elimination()
        )
    }
}
