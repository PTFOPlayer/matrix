#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Size {
    pub x: usize,
    pub y: usize,
}

/// Returns ammount of multiplication operations necessary to multiply matrices.
pub fn non_optimal_multiplication_cost(matrices: Vec<Size>) -> usize {
    let mut counter = 0;

    let mut temp_sizes = Size {
        x: matrices[0].x,
        y: matrices[0].y,
    };

    for i in 1..matrices.len() {
        counter += temp_sizes.x * matrices[i].y;
        temp_sizes.x = temp_sizes.x;
        temp_sizes.y = matrices[i].y;
    }

    counter
}

fn format_optimal(i: i64, j: i64, n: usize, brackets: Vec<Vec<i64>>, num: &mut i64) -> String {
    let mut fmt = String::new();
    if i == j {
        *num += 1;
        fmt += format!("A{}", *num).as_str();
        return fmt;
    }
    let temp = *num;
    if temp != 0 {
        fmt += "(";
    }

    fmt += format_optimal(
        i.clone(),
        brackets.clone()[i as usize][j as usize],
        n,
        brackets.clone(),
        num,
    )
    .as_str();
    fmt += " * ";
    fmt += format_optimal(
        brackets.clone()[i as usize][j as usize] + 1,
        j,
        n,
        brackets.clone(),
        num,
    )
    .as_str();

    if temp != 0 {
        fmt += ")";
    }

    fmt
}


/// Returns string with layout for optimal brackets setup to multiply matrices most efficiently.
pub fn optimal_multiplication_brackets(matrices: Vec<Size>) -> String {
    let len = matrices.len();

    let mut p = Vec::new();
    p.push(matrices[0].x);
    p.push(matrices[0].y);
    for i in 1..len {
        p.push(matrices[i].y);
    }

    let mut m = Vec::new();
    let mut brackets = Vec::new();

    for _ in 0..len {
        let mut t_vec = Vec::new();
        let mut t2_vec = Vec::new();
        for _ in 0..len {
            t_vec.push(0);
            t2_vec.push(0);
        }
        m.push(t_vec.clone());
        brackets.push(t_vec);
    }

    for l in 2..len {
        for i in 1..len - l + 1 {
            let j = l + i - 1;
            m[i][j] = i64::MAX;
            for k in i..j {
                let q = m[i][k] + m[k + 1][j] + p[i - 1] as i64 * p[k] as i64 * p[j] as i64;
                if q < m[i][j] {
                    m[i][j] = q;
                    brackets[i][j] = k as i64;
                }
            }
        }
    }

    let num = &mut 0;
    return format_optimal(0, (len - 1) as i64, len, brackets, num);
}
