use crate::field::{Field};
use crate::finit_diff::poisson_relax;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use crate::PEREODIC;
#[derive(Serialize, Deserialize, Copy, Clone,Debug)]
pub struct ParamsConc {
    pub le: f64,
    pub l: f64,
    #[serde(default ="zero")]
    pub sor:f64,
    #[serde(default ="zero")]
    pub psi:f64
}
fn zero() ->f64{
    0.0
}

#[derive(Serialize, Deserialize, Copy, Clone,Debug)]
pub struct Params {
    pub rel: f64,
    pub bm: f64,
    pub bm2: f64,
    #[serde(default ="zero")]
    pub rel_c1:f64,
    #[serde(default ="zero")]
    pub rel_c2:f64,

    pub pr: f64,
    pub time: f64,
    pub params_c: ParamsConc,
    pub params_c2: ParamsConc,
}
use crate::{H, NX, NY};
pub struct Concentration {
    pub f: Array2<f64>,
    qew: Array2<f64>,
    qsn: Array2<f64>,
    vx: Array2<f64>,
    vy: Array2<f64>,
    pub params: ParamsConc,
}

pub struct Temperatura {
    pub f: Array2<f64>,
    delta: Array2<f64>,
}
pub struct Phi {
    pub f: Array2<f64>,
    delta: Array2<f64>,
    pub params: Params,
}

pub struct Psi {
    pub f: Array2<f64>,
}

impl Field for Temperatura {
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}
impl Field for Phi {
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}
impl Field for Psi {
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}
impl Field for Concentration {
    fn get_f(&self) -> &Array2<f64> {
        return &self.f;
    }
}

impl Field for Array2<f64> {
    fn get_f(&self) -> &Array2<f64> {
        return &self;
    }
}
impl Concentration {
    pub fn new(nx: usize, ny: usize, params: ParamsConc) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
            qew: Array2::<f64>::zeros((nx, ny)),
            qsn: Array2::<f64>::zeros((nx, ny)),
            vx: Array2::<f64>::zeros((nx, ny)),
            vy: Array2::<f64>::zeros((nx, ny)),
            params: params,
        };
    }
    fn calc_v(&mut self, psi: &Psi) {
        unsafe {
            for i in 1..NX {
                let mut k = 0;
                self.vx[[i, k]] = psi.f[[i, k + 1]] + psi.f[[i - 1, k + 1]];

                for k in 1..NY - 1 {
                    self.vx[[i, k]] = psi.dy((i, k)) + psi.dy((i - 1, k));
                }
                k = NY - 1;
                self.vx[[i, k]] = -psi.f[[i, k - 1]] - psi.f[[i - 1, k - 1]];
            }
            self.vx /= 2.0;

            for i in 1..NX - 1 {
                for k in 1..NY {
                    self.vy[[i, k]] = psi.dx((i, k)) + psi.dx((i, k - 1));
                }
            }

            for k in 1..NY {
                {
                    let i = 0;
                    self.vy[[i, k]] = psi.f[[i + 1, k]] + psi.f[[i + 1, k - 1]];
                }
                let i = NX - 1;
                self.vy[[i, k]] = -psi.f[[i - 1, k]] - psi.f[[i - 1, k - 1]];
            }
            self.vy /= 2.0;
           
        }
    }
    pub fn step(&mut self, psi: &Psi,temp:&Temperatura, dt: f64) {
        let le = self.params.le;
        let sor = self.params.sor;
        let l = self.params.l;

        unsafe {
            self.calc_v(&psi);

            for i in 1..NX {
                for k in 0..NY {
                    let mut tmp = (self.vx[[i, k]] ) * self.mx_b((i, k));
                    tmp += -le * self.dx_b((i, k));
                    tmp += le * sor*temp.dx_b((i, k));
                    self.qew[[i, k]] = tmp / H;
                }
            }
            for i in 0..NX {
                for k in 1..NY {
                    let mut tmp = (-self.vy[[i, k]]-le/l) * self.my_b((i, k));
                    tmp += -le * (self.dy_b((i, k)));
                    tmp += le * sor*temp.dy_b((i, k));
                    self.qsn[[i, k]] = tmp / H;
                }
            }

            for i in 1..NX - 1 {
                for k in 1..NY - 1 {
                    let q = -(self.qew.dx_f((i, k)) + self.qsn.dy_f((i, k)));
                    self.f[[i, k]] += dt / H * q;
                }
            }
            for i in 1..NX - 1 {
                let k = 0;
                self.f[[i, k]] += dt / H * (-self.qew.dx_f((i, k)) - 2.0 * self.qsn[[i, k + 1]]);

                let k = NY - 1;
                self.f[[i, k]] += dt / H * (-self.qew.dx_f((i, k)) + 2.0 * self.qsn[[i, k]]);
            }
            if PEREODIC{
                for k in 0..NY {
                    self.f[[0, k]] = self.f[[NX - 2, k]];
                    self.f[[NX - 1, k]] = self.f[[1, k]];
                }
            }else{
                for k in 1..NY - 1 {
                    let i= 0;
                    self.f[[i, k]] += dt / H * (-self.qsn.dy_f((i, k)) - 2.0 * self.qew[[i+1, k]]);
    
                    let i = NX - 1;
                    self.f[[i, k]] += dt / H * (-self.qsn.dy_f((i, k)) +2.0 * self.qew[[i, k]]);
                }
                self.f[[0, 0]] += dt / H * 2.0 *( -self.qew[[1, 0]]- self.qsn[[0, 1]]);
                self.f[[0, NY-1]] +=  dt / H * 2.0 *( -self.qew[[1, NY-1]]+ self.qsn[[0, NY-1]]);
                self.f[[NX-1, 0]] += dt / H * 2.0 *( self.qew[[NX-1, 0]]- self.qsn[[NX-1, 1]]);
                self.f[[NX-1, NY-1]] += dt / H * 2.0 *( self.qew[[NX-1, NY-1]]+ self.qsn[[NX-1, NY-1]]);
            }
        }
    }
}

impl Temperatura {
    pub fn new(nx: usize, ny: usize) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
            delta: Array2::<f64>::zeros((nx, ny)),
        };
    }
    pub fn step(&mut self, psi: &Psi, dt: f64) {
        unsafe {
            for i in 1..NX - 1 {
                for j in 1..NY - 1 {
                    let mut tmp = self.lap((i, j));
                    tmp -= psi.dy((i, j)) * (self.dx((i, j)));
                    tmp += psi.dx((i, j)) * (self.dy((i, j)));
                    tmp /= H * H;
                    *self.delta.uget_mut((i, j)) = tmp * dt;
                }
            }
        }
        self.f += &self.delta;
    }
}
impl Phi {
    pub fn new(nx: usize, ny: usize, params: Params) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
            delta: Array2::<f64>::zeros((nx, ny)),
            params: params,
        };
    }
    pub fn step(&mut self, psi: &Psi, t: &Temperatura, conc: &Concentration, conc2: &Concentration, dt: f64) {
        //let mut delta = Array2::<f64>::zeros((NX,NY));
        //let shape= delta.dim();
        unsafe {
            for i in 1..NX - 1 {
                for j in 1..NY - 1 {
                    let mut tmp = self.lap((i, j)) * self.params.pr;
                    tmp -= psi.dy((i, j)) * (self.dx((i, j)));
                    tmp += psi.dx((i, j)) * (self.dy((i, j)));
                    tmp /= H * H;
                    tmp += self.params.pr * self.params.rel * (t.dx((i, j))) / H;
                    tmp -= self.params.pr * self.params.rel_c1 * (conc.dx((i, j))) / H;
                    tmp -= self.params.pr * self.params.rel_c2 * (conc2.dx((i, j))) / H;
                    *self.delta.uget_mut((i, j)) = tmp * dt;
                }
            }
        }
        self.f += &self.delta;
    }
}
impl Psi {
    pub fn new(nx: usize, ny: usize) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
        };
    }
    pub fn step(&mut self, phi: &Phi) {
        poisson_relax(&phi.f, &mut self.f, H);
    }
}
