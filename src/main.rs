#![allow(unused)]
//pub mod fluidstuff;
//pub mod fluidstuffbetter;
pub mod bestfluid;

//use fluidstuff::fluidobj;
//use fluidstuffbetter::*;
use bestfluid::*;



fn print_arr(arr: Vec<f64>, dim: usize) {
    for i in 0..(arr.len()) {
        //print!("({}) ", i);
        if (i == 0)
        {
            print!("[{:.5}, ", arr[i]);
        } else if ((i + 1) % dim == 0)
        {
            print!("{:.5}] \n", arr[i]);
        } else if ((i) % dim == 0)
        {
            print!("[{:.5}, ", arr[i]);
        } else
        {
            print!("{:.5}, ", arr[i]);
        }
    }
    
}

fn arr_sum(arr: Vec<f64>) -> f64 {
    let mut sum: f64 = 0_f64;
    for val in arr {
        sum += val;
    }
    return sum;

}

fn main() {


    //print_arr(&test.den, ((N+2) as usize));
    //test.full_step(dt);

    // /*

    use std::time::Instant;
    let now = Instant::now();

    let iters = 100;
    {


    let size: i32 = N;
    let dim: usize = (size + 2) as usize;
    let dt: f64 = 0.1;

    let mut test_fluid: b_fluobj = b_fluobj::new(N);

    test_fluid.x[12] = 10_f64;

    print_arr(test_fluid.x.clone(), dim);


    for i in 0..iters {
        test_fluid.update_self(dt);
        let sum = arr_sum(test_fluid.x.clone());
        println!("Sum of densities: {}", sum);
        print_arr(test_fluid.x.clone(), dim);
    }

    }
    let elapsed = now.elapsed();
    println!("For {} interations, time elapsed was: {:.2?}", iters, elapsed);
    // */

}
