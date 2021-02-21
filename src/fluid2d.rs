use std::ops::{Index, IndexMut};

pub struct Array2DF {
    pub M: usize,
    pub N: usize,
    pub data: Vec<f32>,
}

impl Array2DF {
    pub fn new(m: usize, n: usize) -> Self {
        let mut data: Vec<f32> = Vec::new();
        for _ in 0..(m * n) {
            data.push(0.0);
        }
        Self { M: m, N: n, data }
    }
}

impl Index<[usize; 2]> for Array2DF {
    type Output = f32;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.data[index[0] + index[1] * self.N]
    }
}

impl IndexMut<[usize; 2]> for Array2DF {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.data[index[0] + index[1] * self.N]
    }
}

pub struct Fluid2D {
    pub size: usize,
    pub dt: f32,
    pub diff: f32,
    pub visc: f32,

    pub s: Array2DF,
    pub density: Array2DF,

    pub Vx: Array2DF,
    pub Vy: Array2DF,

    pub Vx0: Array2DF,
    pub Vy0: Array2DF,
}

impl Fluid2D {
    pub fn new(size: usize, dt: f32, diff: f32, visc: f32) -> Self {
        Self {
            size,
            diff,
            dt,
            visc,
            s: Array2DF::new(size, size),
            density: Array2DF::new(size, size),
            Vx: Array2DF::new(size, size),
            Vy: Array2DF::new(size, size),
            Vx0: Array2DF::new(size, size),
            Vy0: Array2DF::new(size, size),
        }
    }

    pub fn add_density(&mut self, x: usize, y: usize, radius: usize, amount: f32) {
        for i in x - radius..x + radius {
            for j in y - radius..y + radius {
                if i <= self.size
                    && j <= self.size
                    && (((i - x) * (i - x) + (j - y) * (j - y)) as f32).sqrt() <= radius as f32
                {
                    self.density[[i, j]] += amount;
                }
            }
        }
    }

    pub fn add_velocity(
        &mut self,
        x: usize,
        y: usize,
        radius: usize,
        amount_x: f32,
        amount_y: f32,
    ) {
        for i in x - radius..x + radius {
            for j in y - radius..y + radius {
                if i <= self.size
                    && j <= self.size
                    && (((i - x) * (i - x) + (j - y) * (j - y)) as f32).sqrt() <= radius as f32
                {
                    self.Vx[[i, j]] += amount_x;
                    self.Vy[[i, j]] += amount_y;
                }
            }
        }
    }

    pub fn step(&mut self, iter: usize) {
        diffuse(
            1,
            &mut self.Vx0,
            &mut self.Vx,
            self.visc,
            self.dt,
            iter,
            self.size,
        );
        diffuse(
            2,
            &mut self.Vy0,
            &mut self.Vy,
            self.visc,
            self.dt,
            iter,
            self.size,
        );

        project(
            &mut self.Vx0,
            &mut self.Vy0,
            &mut self.Vx,
            &mut self.Vy,
            iter,
            self.size,
        );

        advect(1, &mut self.Vx, &self.Vx0, &self.Vx0, &self.Vy0, self.dt);
        advect(2, &mut self.Vy, &self.Vy0, &self.Vx0, &self.Vy0, self.dt);

        project(
            &mut self.Vx,
            &mut self.Vy,
            &mut self.Vx0,
            &mut self.Vy0,
            iter,
            self.size,
        );

        diffuse(
            0,
            &mut self.s,
            &mut self.density,
            self.diff,
            self.dt,
            iter,
            self.size,
        );
        advect(
            0,
            &mut self.density,
            &mut self.s,
            &mut self.Vx,
            &mut self.Vy,
            self.dt,
        );
    }
}

pub fn lin_solve(
    b: i32,
    x: &mut Array2DF,
    x0: &mut Array2DF,
    a: f32,
    c: f32,
    iter: usize,
    N: usize,
) {
    let c_recip = 1.0 / c;
    for _ in 0..iter {
        for j in 1..(N - 1) {
            for i in 1..(N - 1) {
                x[[i, j]] = (x0[[i, j]]
                    + a * (x[[i + 1, j]] + x[[i - 1, j]] + x[[i, j + 1]] + x[[i, j - 1]]))
                    * c_recip;
            }
        }
        set_bnd(b, x);
    }
}

pub fn diffuse(
    b: i32,
    x: &mut Array2DF,
    x0: &mut Array2DF,
    diff: f32,
    dt: f32,
    iter: usize,
    N: usize,
) {
    let a = dt * diff * (N - 2) as f32 * (N - 2) as f32;
    lin_solve(b, x, x0, a, 1.0 + 6.0 * a, iter, N);
}

pub fn project(
    veloc_x: &mut Array2DF,
    veloc_y: &mut Array2DF,
    p: &mut Array2DF,
    div: &mut Array2DF,
    iter: usize,
    N: usize,
) {
    for j in 1..(veloc_x.N - 1) {
        for i in 1..(veloc_x.N - 1) {
            div[[i, j]] = -0.5
                * (veloc_x[[i + 1, j]] - veloc_x[[i - 1, j]] + veloc_y[[i, j + 1]]
                    - veloc_y[[i, j - 1]])
                / veloc_x.N as f32;
            p[[i, j]] = 0.0;
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1.0, 6.0, iter, N);

    for j in 1..(veloc_x.N - 1) {
        for i in 1..(veloc_x.N - 1) {
            veloc_x[[i, j]] -= 0.5 * (p[[i + 1, j]] - p[[i - 1, j]]) * N as f32;
            veloc_y[[i, j]] -= 0.5 * (p[[i, j + 1]] - p[[i, j - 1]]) * N as f32;
        }
    }
    set_bnd(1, veloc_x);
    set_bnd(2, veloc_y);
}

pub fn advect(
    b: i32,
    d: &mut Array2DF,
    d0: &Array2DF,
    veloc_x: &Array2DF,
    veloc_y: &Array2DF,
    dt: f32,
) {
    let mut i0;
    let mut i1;
    let mut j0;
    let mut j1;

    let dtx = dt * (d.N - 2) as f32;
    let dty = dt * (d.N - 2) as f32;

    let mut s0;
    let mut s1;
    let mut t0;
    let mut t1;
    let mut tmp1;
    let mut tmp2;
    let mut x;
    let mut y;

    let Nfloat: f32 = veloc_x.N as f32 - 2.0;
    let mut ifloat;
    let mut jfloat;

    for j in 1..(d.N - 1) {
        jfloat = j;
        for i in 1..(d.N - 1) {
            ifloat = i;
            tmp1 = dtx * veloc_x[[i, j]];
            tmp2 = dty * veloc_y[[i, j]];

            x = ifloat as f32 - tmp1;
            y = jfloat as f32 - tmp2;

            if x < 0.5 {
                x = 0.5;
            }
            if x > Nfloat + 0.5 {
                x = Nfloat + 0.5;
            }
            i0 = x.floor();
            i1 = i0 + 1.0;
            if y < 0.5 {
                y = 0.5;
            }
            if y > Nfloat + 0.5 {
                y = Nfloat + 0.5;
            }
            j0 = y.floor();
            j1 = j0 + 1.0;

            s1 = x - i0;
            s0 = 1.0 - s1;
            t1 = y - j0;
            t0 = 1.0 - t1;

            let i0i = i0;
            let i1i = i1;
            let j0i = j0;
            let j1i = j1;

            d[[i, j]] = s0
                * (t0 * d0[[i0i as usize, j0i as usize]]
                    + t1 * d0[[i0i as usize, j1i as usize]] as f32)
                + s1 * (t0 * d0[[i1i as usize, j0i as usize]]
                    + t1 * d0[[i1i as usize, j1i as usize]] as f32);
        }
    }
    set_bnd(b, d);
}

fn set_bnd(b: i32, x: &mut Array2DF) {
    let N = x.N;
    for i in 1..(N - 1) {
        x[[i, 0]] = if b == 2 { -x[[i, 1]] } else { x[[i, 1]] };
        x[[i, N - 1]] = if b == 2 {
            -x[[i, x.M - 2]]
        } else {
            x[[i, x.M - 2]]
        };
    }

    for j in 1..(x.N - 1) {
        x[[0, j]] = if b == 1 { -x[[1, j]] } else { x[[1, j]] };
        x[[N - 1, j]] = if b == 1 {
            -x[[x.N - 2, j]]
        } else {
            x[[x.N - 2, j]]
        };
    }

    x[[0, 0]] = 0.5 * (x[[1, 0]] + x[[0, 1]]);
    x[[0, N - 1]] = 0.5 * (x[[1, x.N - 1]] + x[[0, x.N - 2]]);
    x[[N - 1, 0]] = 0.5 * (x[[x.N - 2, 0]] + x[[x.N - 1, 1]]);
    x[[N - 1, N - 1]] = 0.5 * (x[[x.N - 2, x.N - 1]] + x[[x.N - 1, x.N - 2]]);
}
