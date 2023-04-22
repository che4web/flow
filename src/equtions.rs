use crate::field::{Field};
use crate::finit_diff::poisson_relax;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use crate::PEREODIC;

use crate::io::write_mat;
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
trait Equation{
    fn step(&mut self, psi: &Psi, t: &Temperatura,  dt: f64) ;
    fn write_map(&mut self, time:f64) ;

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


impl Field for Array2<f64> {
    fn get_f(&self) -> &Array2<f64> {
        return &self;
    }
}

unsafe fn conv_mul(v:&Array2<f64>,t:&Array2<f64>,index:(usize,usize))->f64{
    return v.dx(index)*t.dy(index)-v.dy(index)*t.dx(index)
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
                    let mut tmp = self.f.lap((i, j));
                    tmp += conv_mul(&psi.f,&self.f,(i, j));
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
impl Equation for Phi{
    fn step(&mut self, psi: &Psi, t: &Temperatura,  dt: f64) {
        //let mut delta = Array2::<f64>::zeros((NX,NY));
        //let shape= delta.dim();
        unsafe {
            for i in 1..NX - 1 {
                for j in 1..NY - 1 {
                    let mut tmp = self.f.lap((i, j)) * self.params.pr;
                    tmp += conv_mul(&psi.f,&self.f,(i, j));
                    tmp /= H * H;
                    tmp += self.params.pr * self.params.rel * (t.f.dx((i, j))) / H;
                    //tmp -= self.params.pr * self.params.rel_c1 * (conc.dx((i, j))) / H;
                    //tmp -= self.params.pr * self.params.rel_c2 * (conc2.dx((i, j))) / H;
                    *self.delta.uget_mut((i, j)) = tmp * dt;
                }
            }
        }
        self.f += &self.delta;
    }

    fn write_map(&mut self,time:f64){
            write_mat(&self.f,String::from(format!("res_psi/psi_t={:05}",time)),H);

    } 
    


}
