use matrix_lib::matrix;

fn main() {
    let mut m1 = matrix![
    2,1;
    5,3
    ];
    println!("{}", m1);
    m1.inverse();
}
