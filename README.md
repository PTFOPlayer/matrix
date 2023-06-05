# matrix_lib
A library implementing the matrices and basic operations on them

# Want to contribute?:
### My github:
<a href="https://github.com/PTFOPlayer">
<img src="https://cdn.jsdelivr.net/npm/simple-icons@3.0.1/icons/github.svg" height="45px" alt="github" />
</a>

### Project:
<a href="https://github.com/PTFOPlayer/euclides">
<img src="https://cdn.jsdelivr.net/npm/simple-icons@3.0.1/icons/github.svg" height="45px" alt="github" />
</a>

### Support:
<p><a href="https://www.buymeacoffee.com/WhiskyAKM"> <img align="left" src="https://cdn.buymeacoffee.com/buttons/v2/default-yellow.png" height="50" width="210" alt="https://www.buymeacoffee.com/WhiskyAKM" /></a></p>
<br/>
<br/>
<br/>
<br/>

## Macors 
```rust
let m1 = matrix![
        vec![1,2,3], 
        vec![1,2,3],
        vec![1,2,3]
    ];

// this macro uses ',' to separate columns and ';' to separate rows
let m2 = matrix![
        1, 2, 3; 
        4, 5, 6; 
        7, 8, 9
    ];

// this macro is equivalent to matrix
// 1 2 3 
// 4 5 6 
// 7 8 9
// it splits the vector by length provided after "=>" operator
let vector = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
let m3 = matrix![vector => 3];
```