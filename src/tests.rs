#[cfg(test)]
mod classic_tests {
    use crate::matrix;

    #[test]
    fn sum() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6]];
        let m2 = matrix![3, 2, 1; 6, 5, 4];
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
        let m1 = matrix![1, 2, 3; 4, 5, 6];
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
    fn big_determinant() {
        let m1 = matrix![
            vec![4, 2, 7, 4],
            vec![4, 5, 7, 4],
            vec![8, 3, 4, 0],
            vec![8, 1, 9, 3]
        ];
        assert_eq!(120.0, m1.determinant().unwrap());
        let m2 = matrix![
            vec![6, 7, 6, 2, 7],
            vec![2, 5, 4, 5, 5],
            vec![5, 2, 1, 8, 1],
            vec![9, 9, 4, 3, 8],
            vec![8, 1, 6, 5, 6]
        ];
        assert_eq!(3685.0, m2.determinant().unwrap().round());
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

    #[test]
    fn vectorize() {
        let m1 = matrix![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]];
        assert_eq!(vec![1, 4, 7, 2, 5, 8, 3, 6, 9], m1.vectorize());
    }

    #[test]
    fn from_vec() {
        let m1 = matrix![1, 2, 3; 4, 5, 6; 7, 8, 9];
        let v = vec![1, 4, 7, 2, 5, 8, 3, 6, 9];
        assert_eq!(m1, matrix!(v.clone() => 3));
    }
}

#[cfg(test)]
mod multi_threaded_tests {
    use crate::matrix;

    #[tokio::test]
    async fn scalar_mul() {
        let m1 = matrix![
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0;
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
        10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0;
        19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0
        ];
        let res = m1.clone().par_scalar_mul(4000.123456).await;
        assert_eq!(m1.clone().scalar_mul(4000.123456), res);
    }

    #[tokio::test]
    async fn mul() {
        let m1 = matrix![[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]];
        let m2 = matrix![[3, 2], [6, 5], [6, 5], [3, 2], [6, 5], [6, 5]];
        assert_eq!(matrix![[111, 90], [111, 90]], m1.par_mul(m2).await.unwrap());
    }

    #[tokio::test]
    async fn sum() {
        let m1 = matrix![[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]];
        let m2 = matrix![[6, 5, 4, 3, 2, 1], [6, 5, 4, 3, 2, 1]];
        assert_eq!(m1.clone().sum(m2.clone()).unwrap(), m1.par_sum(m2).await.unwrap());
    }
}
