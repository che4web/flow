const NX:usize=32;
const NY:usize=32;
const L:f64=1.0;
const H:f64=L/((NY-1) as f64);
const DT:f64=H*H/8.0;

//pub mod io;
pub mod finit_diff;
pub mod io;
//use io::{writeMat};

use ndarray::{Array2};
use ndarray_stats::QuantileExt;

use csv::Writer;
use serde::{Serialize,Deserialize};


use finit_diff::{dx,dy,laplace,poisson_relax,dx_f,dy_f,dx_b,dy_b,meanCx,meanCy};
use io::writeMat;
use std::fs;

use toml;

#[derive(Serialize,Deserialize)]
struct Params{
    R:f64
}
struct System{
    T:Array2<f64>,
    psi:Array2<f64>,
    phi:Array2<f64>,
    C:Array2<f64>,
    params:Params,
}
impl System{
    fn step_t(&self)->Array2<f64> {
        let mut delta = Array2::<f64>::zeros((NX,NY));
        let mut jx = Array2::<f64>::zeros((NX,NY));
        let mut jy = Array2::<f64>::zeros((NX,NY));
        delta+=&(laplace(&self.T));
        delta-=&(dy(&self.psi)*dx(&self.T));
        delta+=&(dx(&self.psi)*dy(&self.T));
        delta/=(H*H);
        return delta
    }
    fn step_phi(&self)->Array2<f64> {
        let mut delta = Array2::<f64>::zeros((NX,NY));
        delta+=&(laplace(&self.phi));
        delta-=&(dy(&self.psi)*dx(&self.phi));
        delta+=&(dx(&self.psi)*dy(&self.phi));
        delta/=(H*H);
        delta+=&(self.params.R*dx(&self.T)/H);
        return delta
    }
    fn step_c(&self)->Array2<f64> {
        let mut delta = Array2::<f64>::zeros((NX,NY));
        let mut jx = Array2::<f64>::zeros((NX,NY));
        let mut jy = Array2::<f64>::zeros((NX,NY));
        let vx=dy(&self.psi);
        let vy=dx(&self.psi);

        let cx = dx_b(&self.C);
        let cy = dy_b(&self.C);
        let tx = dx_b(&self.T);
        let ty = dy_b(&self.T);
        let mut c_loc= meanCx(&self.C);

        jx+=&tx;
        jx+=&cx;
        jx+=&(dx_b(&vx)*&c_loc);

        c_loc= meanCy(&self.C);
        jy+=&ty;
        jy+=&cy;
        jy+=&(dy_b(&vx)*&c_loc);

        for i in 1..NX{
            let c_l = C[[i,0]]+C[[i-1,0]];
            jx[[i,0]] = (psi[[i,1]]+psi[[i-1,1]])*2*c_l;
            jx[[i,0]] += cx[[i,0]];
            jx[[i,0]] += tx[[i,0]]*c_l;

            c_l = C[[i,NY-1]]+C[[i-1,NY-1]];

            jx[[i,NY-1]] = (-psi[[i,NY-2]]-psi[[i-1,NY-2]])*2*c_l;
            jx[[i,NY-1]] += cx[[i,NY-1]];
            jx[[i,NY-1]] += tx[[i,NY-1]]*c_l;
        
        }

        delta+=&(&dx_f(&jx)+&dy_f(&jy));
        return delta;
    }

    fn step_psi(&self)->Array2<f64>{
        poisson_relax(&self.phi,&self.psi,H)
    }

    fn next_step(&mut self,_DT:f64){
        self.phi+= &(self.step_phi()*_DT);
        self.psi = self.step_psi();
        self.T+=&(self.step_t()*_DT);
        self.boundary_condition();
    }
    fn boundary_condition(&mut self){

        for i in 1..NY-1{
            self.phi[[0,i]]=-self.psi[[1,i]]/(H*H)*2.0;
            self.phi[[NX-1,i]]=-self.psi[[NX-2,i]]/(H*H)*2.0;
        }

        for i in 1..NX-1{
            self.phi[[i,0]]=-self.psi[[i,1]]/(H*H)*2.0;
            self.phi[[i,NY-1]]=-self.psi[[i,NY-2]]/(H*H)*2.0;
        }

    }
}
fn log_params(time:f64,system: &System)->Row{
    let mut nu  =0.0;
    let shape = system.T.dim();
    for i in 1..shape.0-1{
        nu+= system.T[[i,shape.1-1]]-system.T[[i,shape.1-2]]
    }
    nu/=H*(shape.0 as f64)-2.0;

   return Row{
        t:time,
        psi_m:*(system.psi.max().unwrap()),
        nu:nu,
   } 
}


#[derive(Serialize)]
pub struct Row{
    t:f64,
    psi_m:f64,
    nu:f64,
}


fn main() {
    println!("Hello, world!");
    const R:f64 = 6000.0;
    //let mut dela =Array2::<f64>::zeros((NX,NY));
    let contents = fs::read_to_string("config.toml").unwrap() ;
    let params:Params = toml::from_str(&contents).unwrap();
    println!("{:?}",params.R);

    let mut system = System{
        psi:Array2::<f64>::zeros((NX,NY)),
        phi:Array2::<f64>::zeros((NX,NY)),
        T:Array2::<f64>::zeros((NX,NY)),
        C:Array2::<f64>::zeros((NX,NY)),
        params:params,
    };
    //set initial 
    system.phi[[NX/2,NY/2]]=1.0;
    for i in 0..NX{
        for j in 0..NY{
            system.T[[i,j]] = 1.0-(j as f64)*H;
        }
    }
    // log open
    let mut wtr = Writer::from_path("foo.csv").unwrap();

    let mut time= 0.0;
    for _x in 0..10{
        system.next_step(DT);
        time+=DT;

        if _x % 10==0{

            let res =  log_params(time,&system);
            wtr.serialize(res).unwrap();
        }

    }
    println!("==============");
    println!("{:?}",system.psi);
    println!("{:?}",system.phi);

    writeMat(&system.phi,String::from("phi"),H);
    writeMat(&system.psi,String::from("psi"),H);
    writeMat(&system.T,String::from("t"),H);
}
