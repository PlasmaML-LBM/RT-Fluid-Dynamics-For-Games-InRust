use std::mem::swap;

pub const N: i32 = 10;

const Nex: i32 = (N + 2) * (N + 2);

pub fn IX(i: i32, j: i32) -> usize {
    (i + (N+2)*j) as usize
}

pub struct fluid_object {
    pub x: Vec<f64>, 
    x0: Vec<f64>, /* x and x0 are density and previous density */
    u: Vec<f64>,
    u0: Vec<f64>,
    v: Vec<f64>,
    v0: Vec<f64>,
    diff: f64,
    visc: f64,
    n: i32,
    size: i32,
}

fn vnew(n: i32) -> Vec<f64> {
    vec![0_f64; ((n+2) * (n+2)) as usize]
}

fn testvec(n: i32) -> Vec<f64> {
    let mut v: Vec<f64> = vec![0_f64; ((n+2) * (n+2)) as usize];
    v[28] = 10_f64;
    v
}

impl fluid_object {
    pub fn new(n: i32) -> Self {
        let mut f = fluid_object {
            x: vnew(n),
            //x: testvec(n), 
            // x0: testvec(n), /* x and x0 are density and previous density */
            x0: vnew(n),
            u: vnew(n),
            u0: vnew(n),
            v: vnew(n),
            v0: vnew(n),
            diff: 0.5_f64,
            visc: 0.5_f64,
            n: n,
            size: (n+2) * (n+2),
        };
        f.x[28] = 100_f64;
        return f;
    }
    pub fn clear_prev_vals(&mut self) {
        /*
        for i in 0..(self.size as usize) {
            self.x0[i as usize] = 0_f64;
            self.u0[i as usize] = 0_f64;
            self.v0[i as usize] = 0_f64;
        }
        // */

        // /*
        self.x0 = vec![0.000_f64; self.size as usize];
        self.u0 = vec![0.000_f64; self.size as usize];
        self.v0 = vec![0.000_f64; self.size as usize];
        // */
    }

    /*
    fn set_val(&mut self, i: i32, j: i32, val: f64) {

    }
    // */

    /*
    pub fn add_input(&mut self, x_s: &mut Vec<f64>, u_s: &mut Vec<f64>, v_s: &mut Vec<f64>) {
        add_source(self.n, &mut self.x0, x_s, 1_f64);
        add_source(self.n, &mut self.u0, u_s, 1_f64);
        add_source(self.n, &mut self.v0, v_s, 1_f64);
    }
    // */
    pub fn vel_step_self(&mut self, dt: f64) {
        //(n: i32, u: &mut Vec<f64>, v: &mut Vec<f64>, u0: &mut Vec<f64>, v0: &mut Vec<f64>, visc: f64, dt: f64)
        vel_step(self.n, &mut self.u, &mut self.v, &mut self.u0, &mut self.v0, self.visc, dt);
    }
    pub fn dens_step_self(&mut self, dt: f64) {
        //(n: i32, x: &mut Vec<f64>, x0: &mut Vec<f64>, u: &mut Vec<f64>, v: &mut Vec<f64>, diff: f64, dt: f64)
        dens_step(self.n, &mut self.x, &mut self.x0, &mut self.u, &mut self.v, self.diff, dt);
    }

    pub fn full_step(&mut self, dt: f64) {
        //self.add_input(&mut vec![0_f64; self.size as usize], &mut vec![0_f64; self.size as usize], &mut vec![0_f64; self.size as usize]);
        self.clear_prev_vals();
        self.vel_step_self(dt);
        self.dens_step_self(dt);
        //self.clear_prev_vals();
    }
}

fn add_source(n: i32, x: &mut Vec<f64>, s: &mut Vec<f64>, dt: f64) {
    let size = (n+2) * (n+2);
    for i in 0..size {
        x[i as usize] += dt * s[i as usize];
    }
}

fn set_bnd(n: i32, b: i32, x: &mut Vec<f64>) {
    for i in 1..=n {
        x[IX(0, i)] = if (b == 1) {-x[IX(1, i)]} else {x[IX(1, i)]};
        x[IX(N+1, i)] = if (b == 1) {-x[IX(N, i)]} else {x[IX(N, i)]};
        x[IX(i, 0)] = if (b == 2) {-x[IX(i, 1)]} else {x[IX(i, 1)]};
        x[IX(i, N+1)] = if (b == 2) {-x[IX(i, N)]} else {x[IX(i, N)]};    
    }
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N+1)] = 0.5 * (x[IX(1, N+1)] + x[IX(0, N)]);
    x[IX(N+1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N+1, 1)]);
    x[IX(N+1, N+1)] = 0.5 * (x[IX(N, N+1)] + x[IX(N+1, N)]);
}

fn lin_solve(n: i32, b: i32, x: &mut Vec<f64>, x0: &mut Vec<f64>, a: f64, c: f64) {
    for k in 0..20 {
        for i in 1..=n {
            for j in 1..=n {
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
            }
        }
        set_bnd(n, b, x)
    }
}

fn diffuse(n: i32, b: i32, x: &mut Vec<f64>, x0: &mut Vec<f64>, diff: f64, dt: f64) {
    let a: f64 = dt * diff * (n as f64) * (n as f64);
    lin_solve(n, b, x, x0, a, 1_f64 + 4_f64*a);
}

fn advect(n: i32, b: i32, d: &mut Vec<f64>, d0: Vec<f64>, u: Vec<f64>, v: Vec<f64>, dt: f64) {

    let mut i0: i32; let mut j0: i32; let mut i1: i32; let mut j1: i32;

    let mut x: f64; let mut y: f64; let mut s0: f64; let mut t0: f64; let mut s1: f64; let mut t1: f64;

    let dt0 = dt * n as f64;

    for i in 1..=n {
        for j in 1..=n {
            x = i as f64 - dt0*u[IX(i, j)];
            y = j as f64 - dt0*v[IX(i, j)];
            
            //x = if (x < 0.5_f64) {0.5_f64} else if (x > N as f64 + 0.5_f64) {N as f64 + 0.5_f64} else {x};
            
            if (x < 0.5_f64) {x = 0.5_f64} else if (x > n as f64 + 0.5_f64) {x = n as f64 + 0.5_f64};
            i0 = x as i32;
            i1 = i0 + 1;

            if (y < 0.5_f64) {y = 0.5_f64} else if (y > n as f64 + 0.5_f64) {y = n as f64 + 0.5_f64};
            j0 = y as i32;
            j1 = j0 + 1;

            s1 = x - i0 as f64;
            s0 = 1_f64 - s1;
            t1 = y - j0 as f64;
            t0 = 1_f64 - t1;

            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
					 s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
        }
    }
    set_bnd(n, b, d);
}

fn project(n: i32, u: &mut Vec<f64>, v: &mut Vec<f64>, p: &mut Vec<f64>, div: &mut Vec<f64>) {
    for i in 1..=n {
        for j in 1..=n {
            div[IX(i,j)] = -0.5_f64*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/(n as f64);
		    p[IX(i,j)] = 0.5_f64;
        }
    }
    set_bnd(n, 0, div); set_bnd(n, 0, p);

    lin_solve(n, 0, p, div, 1_f64, 4_f64);

    for i in 1..=n {
        for j in 1..=n {
            u[IX(i,j)] -= 0.5_f64*(n as f64)*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		    v[IX(i,j)] -= 0.5_f64*(n as f64)*(p[IX(i,j+1)]-p[IX(i,j-1)]);
        }
    }

    set_bnd(n, 1, u); set_bnd(n, 2, v);
}

fn dens_step(n: i32, x: &mut Vec<f64>, x0: &mut Vec<f64>, u: &mut Vec<f64>, v: &mut Vec<f64>, diff: f64, dt: f64) {
    add_source(n, x, x0, dt);
    swap(x0, x);
    diffuse(n, 0, x, x0, diff, dt);
    swap(x0, x);
    advect(n, 0, x, x0.clone(), u.clone(), v.clone(), dt);
}

fn vel_step(n: i32, u: &mut Vec<f64>, v: &mut Vec<f64>, u0: &mut Vec<f64>, v0: &mut Vec<f64>, visc: f64, dt: f64) {
    add_source(n, u, u0, dt); add_source(n, v, v0, dt);
    swap(u0, u); diffuse(n, 1, u, u0, visc, dt);
    swap(v0, v); diffuse(n, 2, v, v0, visc, dt);
    project(n, u, v, u0, v0);
    swap(u0, u); swap(v0, v);
    advect(n, 1, u, u0.clone(), u0.clone(), v0.clone(), dt); advect(n, 2, v, v0.clone(), u0.clone(), v0.clone(), dt);
    project(n, u, v, u0, v0);
}

pub fn print_matrix(mat: Vec<f64>) {
    for i in 0..(N+2) {
        let start_index = (i * (N+2)) as usize;
        let end_index = ((i+1) * (N+2)) as usize;
        let thing = &mat[start_index..end_index];
        let colored = if (i == 0 || i == N+1) {false} else {true};
        //println!("{:?}", thing);
        print_slice(thing, colored);
    }
    println!();
}

fn print_slice(s: &[f64], colored: bool) {
    print!("[ {number:.prec$}, ", number = s[0], prec = 5);
    for i in (1..s.len()-1) {
        if (colored) {
            print!("\x1b[94m{number:.prec$}, \x1b[0m", number = s[i as usize], prec = 5);
        }
        else {
            print!("{number:.prec$}, ", number = s[i as usize], prec = 5);
        }
    }
    print!("{number:.prec$} ]\n", number = s[(s.len() - 1) as usize], prec = 5);
}