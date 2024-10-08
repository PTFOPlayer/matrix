#[cfg(test)]
mod classic_tests {
    use crate::matrix;

    #[test]
    fn sum() {
        let m1 = matrix![[1, 2, 3], [4, 5, 6]];
        let m2 = matrix![3, 2, 1; 6, 5, 4];
        assert_eq!(matrix![[4, 4, 4], [10, 10, 10]], m1.sum(&m2).unwrap())
    }

    #[test]
    fn sum_fail() {
        let m1 = matrix![[1, 2, 3], [4, 5, 6]];
        let m2 = matrix![[3, 2], [6, 5]];
        assert!(m1.sum(&m2).is_err());
    }

    #[test]
    fn mul() {
        let m1 = matrix![[1, 2, 3], [4, 5, 6]];
        let m2 = matrix![[3, 2], [6, 5], [6, 5]];
        assert_eq!(matrix![[33, 27], [78, 63]], m1.mul(m2).unwrap())
    }
    #[test]
    fn scalar_mul() {
        let m1 = matrix![[1, 2, 3], [4, 5, 6]];
        assert_eq!(matrix![[10, 20, 30], [40, 50, 60]], m1.scalar_mul(10))
    }

    #[test]
    fn transpose() {
        let m1 = matrix![1, 2, 3; 4, 5, 6];
        assert_eq!(matrix![[1, 4], [2, 5], [3, 6]], m1.transpose())
    }

    // #[test]
    // fn determinant() {
    //     let mut m1 = matrix![[1, 2, 3], [4, 5, 6], [7, 8, 9]];
    //     assert_eq!(0 as f64, m1.determinant().unwrap());
    //     let mut m2 = matrix![[5, 4, 1], [9, 2, 2], [4, 4, 0]];
    //     assert_eq!(20 as f64, m2.determinant().unwrap());
    // }

    // #[test]
    // fn big_determinant() {
    //     let mut m1 = matrix![[4, 2, 7, 4], [4, 5, 7, 4], [8, 3, 4, 0], [8, 1, 9, 3]];
    //     assert_eq!(120.0, m1.determinant().unwrap());
    //     let mut m2 = matrix![
    //         [6, 7, 6, 2, 7],
    //         [2, 5, 4, 5, 5],
    //         [5, 2, 1, 8, 1],
    //         [9, 9, 4, 3, 8],
    //         [8, 1, 6, 5, 6]
    //     ];
    //     assert_eq!(3685.0, m2.determinant().unwrap().round());
    // }

    #[test]
    fn gaussian_elimination() {
        let mut m1 = matrix![[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        assert_eq!(
            matrix![[1.0, 2.0, 3.0], [0.0, -3.0, -6.0], [0.0, 0.0, 0.0]],
            m1.gaussian_elimination()
        );
    }

    #[test]
    fn gaussian_elimination_with_division_by_zero() {
        let mut m1 = matrix![[1, -2, 3], [-5, 10, 14], [3, -6, -9]];
        assert_eq!(
            matrix![[1.0, -2.0, 3.0], [0.0, 0.0, 29.0], [0.0, 0.0, -18.0]],
            m1.gaussian_elimination()
        )
    }

    #[test]
    fn gaussian_elimination_with_row_swap() {
        let mut m1 = matrix![
            1 1 1 1;
            0 0 1 1;
            0 1 1 1;
            0 0 0 1
        ];
        assert_eq!(
            matrix![
                1.0 1.0 1.0 1.0;
                0.0 1.0 1.0 1.0;
                0.0 0.0 1.0 1.0;
                0.0 0.0 0.0 1.0
            ],
            m1.gaussian_elimination()
        );
    }

    #[test]
    fn determinant_cmp_gaussian_elimination_determinant() {
        let mut m1 = matrix![[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        let mut m2 = m1.clone();
        assert_eq!(
            m1.determinant().unwrap() as f64,
            m2.gaussian_elimination().determinant().unwrap()
        );
    }

    #[test]
    fn determinant_with_row_swap() {
        let mut m1 = matrix![
            1 1 1 1;
            0 0 1 1;
            0 1 1 1;
            0 0 0 1
        ];
        assert_eq!(-1.0, m1.determinant().unwrap());
    }

    #[test]
    fn vectorize() {
        let m1 = matrix![[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        assert_eq!(vec![1, 2, 3, 4, 5, 6, 7, 8, 9], m1.vectorize());
    }

    #[test]
    fn from_vec() {
        let m1 = matrix![1, 2, 3; 4, 5, 6; 7, 8, 9];
        let v = [1, 4, 7, 2, 5, 8, 3, 6, 9];
        assert_eq!(m1, matrix!(v.clone() => 3));
    }

    #[test]
    fn inverse() {
        let m1 = matrix![
        2,1;
        5,3
        ];

        assert_eq!(matrix![3.0, -1.0; -5.0, 2.0], m1.inverse().unwrap());
    }

    #[test]
    fn inverse_2() {
        let m1 = matrix![
            1 1 1 1;
            0 0 1 1;
            0 1 1 1;
            0 0 0 1
        ];

        assert_eq!(
            matrix![
                [1.0, 0.0, -1.0, 0.0],
                [0.0, -1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0, -1.0],
                [0.0, 0.0, 0.0, 1.0]
            ],
            m1.inverse().unwrap()
        );
    }
}
