pub const N: i32 = 3;
pub const Nex: i32 = (N + 2) * (N + 2);
//const BOUNDARY_INDICES: [usize; (((N+2)*(N+2)) - (N*N)) as usize] = [0; (((N+2)*(N+2)) - (N*N)) as usize];

pub struct fluidobj {
    pub p_den: [f64; Nex as usize],
    pub den: [f64; Nex as usize],
    pub u: [f64; Nex as usize],
    pub v: [f64; Nex as usize],
    pub p_u: [f64; Nex as usize],
    pub p_v: [f64; Nex as usize],
    pub diff: f64,
    pub visc: f64,
}

impl fluidobj {
    pub fn full_step(&mut self, dt: f64) {
        vel_step(
            &mut self.u,
            &mut self.v,
            &mut self.p_u,
            &mut self.p_v,
            self.visc,
            dt,
        );
        dens_step(
            &mut self.den,
            &mut self.p_den,
            &mut self.u,
            &mut self.v,
            self.diff,
            dt,
        );
    }
    pub fn new() -> fluidobj {
        let mut x = fluidobj {
            p_den: [0_f64; Nex as usize],
            den: [0_f64; Nex as usize],
            u: [0_f64; Nex as usize],
            v: [0_f64; Nex as usize],
            p_u: [0_f64; Nex as usize],
            p_v: [0_f64; Nex as usize],
            diff: 0.1_f64,
            visc: 0.1_f64,
        };
        //x.den[IX(N/2, N/2)] = 1_f64;
        return x;
    }
}
/*
fn thingy() {
    println!("Hello, world!");

    /*

    let mut a = [0_f64, 0_f64, 1_f64, 3_f64];

    let mut b = [3_f64, 7_f64, 7_f64, 2_f64];

    println!("a: {:?}, b: {:?}", a, b);

    std::mem::swap(&mut a, &mut b);

    println!("a: {:?}, b: {:?}", a, b);

    // */

    ///*
    let mut p_den = [0_f64; Nex as usize];
    let mut den = [0_f64; Nex as usize];
    let mut u = [100_f64; Nex as usize];
    let mut v = [0_f64; Nex as usize];
    let mut p_u = [0_f64; Nex as usize];
    let mut p_v = [0_f64; Nex as usize];
    let diff: f64 = 0.1;
    let visc: f64 = 0.1;
    let dt: f64 = 0.1;

    den[8] = 1.2;

    for i in 0..10 {
        println!("DEN:");
        print_matrix(den);
        print!("u:");
        print_matrix(u);
        print!("v:");
        print_matrix(v);

        //diffuse(0, &mut d, s, 0.1_f64, 0.8_f64);
        vel_step(&mut u, &mut v, &mut p_u, &mut p_v, visc, dt);
        dens_step(&mut den, &mut p_den, &mut u, &mut v, diff, dt);
    }

    //println!("Hello {} is {number:.prec$}", "x", prec = 5, number = 0.01);

    // */
}
*/
pub fn IX(i: i32, j: i32) -> usize {
    return ((i) + (N + 2) * j) as usize;
}

fn add_source(destination: &mut [f64; Nex as usize], source: &mut [f64; Nex as usize], dt: f64) {
    let size = Nex;

    for i in 0..size {
        destination[i as usize] += dt * source[i as usize];
    }
}

fn diffuse(
    b: i32,
    den: &mut [f64; Nex as usize],
    p_den: &mut [f64; Nex as usize],
    diff: f64,
    dt: f64,
) {
    let a: f64 = dt * diff * ((N * N) as f64);
    for k in 0..20 {
        for i in 1..=N {
            for j in 1..=N {
                den[IX(i, j)] = (p_den[IX(i, j)]
                    + a * (den[IX(i - 1, j)]
                        + den[IX(i + 1, j)]
                        + den[IX(i, j - 1)]
                        + den[IX(i, j + 1)]))
                    / (1 as f64 + 4 as f64 * a);
            }
        }
        set_bnd(b, den);
    }
}

// /* note that it technically labels it as den, p_den, u_vel and v_vel but it's more generalized to arr, p_arr, p_u, p_v since you can advect the velocities as well
fn advect(
    b: i32,
    arr: &mut [f64; Nex as usize],
    p_arr: [f64; Nex as usize],
    u_vel: [f64; Nex as usize],
    v_vel: [f64; Nex as usize],
    dt: f64,
) {
    let dt0: f64 = dt * N as f64;

    for i in 1..=N {
        for j in 1..=N {
            let x = i as f64 - dt0 * u_vel[IX(i, j)];
            let y = j as f64 - dt0 * v_vel[IX(i, j)];

            let x = if (x < 0.5) {
                0.5
            } else if (x > N as f64 + 0.5) {
                N as f64 + 0.5
            } else {
                x
            };
            let y = if (x < 0.5) {
                0.5
            } else if (x > N as f64 + 0.5) {
                N as f64 + 0.5
            } else {
                x
            };

            let i0 = x as i32;
            let i1 = i0 + 1;
            let j0 = y as i32;
            let j1 = j0 + 1;

            let s1 = x - i0 as f64;
            let s0 = 1_f64 - s1;

            let t1 = y - j0 as f64;
            let t0 = 1_f64 - t1;

            arr[IX(i, j)] = s0 * (t0 * p_arr[IX(i0, j0)] + t1 * p_arr[IX(i0, j1)])
                + s1 * (t0 * p_arr[IX(i1, j0)] + t1 * p_arr[IX(i1, j1)]);
        }
    }
    set_bnd(b, arr);
}
// */
fn project(
    u: &mut [f64; Nex as usize],
    v: &mut [f64; Nex as usize],
    p: &mut [f64; Nex as usize],
    div: &mut [f64; Nex as usize],
) {
    let h: f64 = 1_f64 / N as f64;
    for i in 1..=N {
        for j in 1..=N {
            div[IX(i, j)] =
                -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]);
            p[IX(i, j)] = 0_f64;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);

    for k in 0..20 {
        for i in 1..=N {
            for j in 1..=N {
                p[IX(i, j)] = (div[IX(i, j)]
                    + p[IX(i - 1, j)]
                    + p[IX(i + 1, j)]
                    + p[IX(i, j - 1)]
                    + p[IX(i, j + 1)])
                    / 4_f64
            }
        }
        set_bnd(0, p);
    }

    for i in 1..=N {
        for j in 1..=N {
            u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
            v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
}

pub fn dens_step(
    den: &mut [f64; Nex as usize],
    p_den: &mut [f64; Nex as usize],
    u_vel: &mut [f64; Nex as usize],
    v_vel: &mut [f64; Nex as usize],
    diff: f64,
    dt: f64,
) {
    //add_source(den, p_den, dt);
    add_source(p_den, den, dt);
    std::mem::swap(den, p_den);
    diffuse(0, den, p_den, diff, dt);
    std::mem::swap(den, p_den);
    advect(0, den, *p_den, *u_vel, *v_vel, dt);
}

pub fn vel_step(
    u_vel: &mut [f64; Nex as usize],
    v_vel: &mut [f64; Nex as usize],
    p_u_vel: &mut [f64; Nex as usize],
    p_v_vel: &mut [f64; Nex as usize],
    visc: f64,
    dt: f64,
) {
    add_source(u_vel, p_u_vel, dt);
    add_source(v_vel, p_v_vel, dt);

    //
    std::mem::swap(v_vel, p_v_vel);
    std::mem::swap(u_vel, p_u_vel);

    diffuse(1, u_vel, p_u_vel, visc, dt);
    diffuse(2, v_vel, p_v_vel, visc, dt);
    //

    project(u_vel, v_vel, p_u_vel, p_v_vel);

    //
    std::mem::swap(v_vel, p_v_vel);
    std::mem::swap(u_vel, p_u_vel);

    advect(1, u_vel, *p_u_vel, *p_u_vel, *p_v_vel, dt);
    advect(2, v_vel, *p_v_vel, *p_u_vel, *p_v_vel, dt);
    //

    project(u_vel, v_vel, p_u_vel, p_v_vel);
}

fn set_bnd(b: i32, arr: &mut [f64; Nex as usize]) {
    for i in 1..=N {
        arr[IX(0, i)] = if (b == 1) {
            -arr[IX(1, i)]
        } else {
            arr[IX(1, i)]
        };
        arr[IX(N + 1, i)] = if (b == 1) {
            -arr[IX(N, i)]
        } else {
            arr[IX(N, i)]
        };
        arr[IX(i, 0)] = if (b == 2) {
            -arr[IX(i, 1)]
        } else {
            arr[IX(i, 1)]
        };
        arr[IX(i, N + 1)] = if (b == 2) {
            -arr[IX(i, N)]
        } else {
            arr[IX(i, N)]
        };
    }
    arr[IX(0, 0)] = 0.5 * (arr[IX(1, 0)] + arr[IX(0, 1)]);
    arr[IX(0, N + 1)] = 0.5 * (arr[IX(1, N + 1)] + arr[IX(0, N)]);
    arr[IX(N + 1, 0)] = 0.5 * (arr[IX(N, 0)] + arr[IX(N + 1, 1)]);
    arr[IX(N + 1, N + 1)] = 0.5 * (arr[IX(N, N + 1)] + arr[IX(N + 1, N)]);
}


/*
fn print_matrix(mat: [f64; Nex as usize]) {
    for i in 0..(N + 2) {
        let start_index = (i * (N + 2)) as usize;
        let end_index = ((i + 1) * (N + 2)) as usize;
        let thing = &mat[start_index..end_index];
        let colored = if (i == 0 || i == N + 1) { false } else { true };
        //println!("{:?}", thing);
        print_slice(thing, colored);
    }
    println!();
}

fn print_slice(s: &[f64], colored: bool) {
    print!("[ {number:.prec$}, ", number = s[0], prec = 5);
    for i in (1..s.len() - 1) {
        if (colored) {
            print!(
                "\x1b[94m{number:.prec$}, \x1b[0m",
                number = s[i as usize],
                prec = 5
            );
        } else {
            print!("{number:.prec$}, ", number = s[i as usize], prec = 5);
        }
    }
    print!(
        "{number:.prec$} ]\n",
        number = s[(s.len() - 1) as usize],
        prec = 5
    );
}
*/
