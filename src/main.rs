const NY:usize=33;
const NX:usize=(NY-1)*4+1;
const L:f64=4.0;
const H:f64=L/((NX-1) as f64);
const DT:f64=H*H/8.0;


//pub mod io;
pub mod finit_diff;
pub mod io;
//use io::{writeMat};

use ndarray::{Array2};
use ndarray_stats::QuantileExt;

use csv::Writer;
use serde::{Serialize,Deserialize};


use finit_diff::{dx,dy,laplace,poisson_relax,dx_f,dy_f,dx_b,dy_b,mean_cx,mean_cy};
use io::write_mat;
use std::fs;

use std::time::{Duration, Instant};
use toml;

#[derive(Serialize,Deserialize,Copy,Clone)]
struct ParamsConc{
    Le:f64,
    l:f64,
}

#[derive(Serialize,Deserialize,Copy,Clone)]
struct Params{
    R:f64,
    B:f64,
    P:f64,
    time:f64,
    paramsC:ParamsConc,

}
struct System{
    T:Array2<f64>,
    psi:Array2<f64>,
    phi:Array2<f64>,
    C:Array2<f64>,
    vx:Array2<f64>,
    vy:Array2<f64>,
    params:Params,
    solverConc:SolverConc,
}
struct SolverConc{
    vx:Array2<f64>,
    vy:Array2<f64>,
    jx:Array2<f64>,
    jy:Array2<f64>,
    djx:Array2<f64>,
    djy:Array2<f64>,
    delta:Array2<f64>,
    c_loc:Array2<f64>,

}
impl SolverConc {
    pub fn new()->Self{
        return Self{
            vx:Array2::<f64>::zeros((NX,NY)),
            vy:Array2::<f64>::zeros((NX,NY)),
            jx:Array2::<f64>::zeros((NX,NY)),
            jy:Array2::<f64>::zeros((NX,NY)),
            djx:Array2::<f64>::zeros((NX,NY)),
            djy:Array2::<f64>::zeros((NX,NY)),
            delta:Array2::<f64>::zeros((NX,NY)),
            c_loc:Array2::<f64>::zeros((NX,NY)),
        }
    }
    fn calc_v(&mut self,psi:&Array2<f64>){
        self.vx=-dy(&psi);
        self.vy=dx(&psi);
        for i in 1..NX{
            self.vx[[i,0]] = -psi[[i,1]];
            self.vx[[i,NY-1]] = psi[[i,NY-2]];
        }

        for j in 1..NY{
            self.vy[[0,j]] = psi[[1,j]];
            self.vy[[NX-1,j]] = -psi[[NX-2,j]];
        }
        self.vx = mean_cx(&self.vx);
        self.vy = mean_cy(&self.vy);
    }
    fn solve(&mut self,psi:&Array2<f64>,t:&Array2<f64>,c:&Array2<f64>,params:&ParamsConc)->&Array2<f64> {

        let le = params.Le;
        let l = params.l;
        self.delta.fill(0.0);
        self.jx.fill(0.0);
        self.jy.fill(0.0);
        self.calc_v(&psi);


        self.c_loc= mean_cx(&c);

        self.jx+=&(dx_b(&t)*0.0);
        self.jx+=&(dx_b(&c)*le);
        self.jx+=&(&self.vx*&self.c_loc);

        self.c_loc= mean_cy(&c);
        self.jy+=&(dy_b(&t)*0.0);
        self.jy+=&(dy_b(&c)*le);
        self.jy+=&(&self.vy*&self.c_loc);
        self.jy+=&(le/l*&self.c_loc*H);

        dx_f(&self.jx,&mut self.djx);
        dy_f(&self.jy,&mut self.djy);

        for j in 1..NY{
            self.djx[[0,j]] = 2.0*self.jx[[1,j]];
            self.djx[[NX-1,j]] = -2.0*self.jx[[NX-1,j]];
            //djy[[NX-1,j]] = -0.0;
        }

        for i in 1..NX{
            self.djy[[i,0]] = 2.0*(self.jy[[i,1]]);
            self.djy[[i,NY-1]] = -2.0*(self.jy[[i,NY-1]]);
        }

        // djx = dx_f(&jx);
        // djy = dy_f(&jy);


        self.delta+=&(&self.djx+&self.djy);
        self.delta[[0,0]]=-2.0*(-self.jx[[1,0]]-self.jy[[0,1]]);
        self.delta[[0,NY-1]]=-2.0*(-self.jx[[1,NY-1]]+self.jy[[0,NY-1]]);
        self.delta[[NX-1,0]]=-2.0*(self.jx[[NX-1,0]]-self.jy[[NX-1,1]]);
        self.delta[[NX-1,NY-1]]=-2.0*(self.jx[[NX-1,NY-1]]+self.jy[[NX-1,NY-1]]);
        
        self.delta/=H*H;
        return &self.delta;
    }

}
impl System{
    fn step_t(&self)->Array2<f64> {
        let mut delta = Array2::<f64>::zeros((NX,NY));
        delta+=&(laplace(&self.T));
        //delta-=&(dy(&self.psi)*dx(&self.T));
       // delta+=&(dx(&self.psi)*dy(&self.T));
       //
        delta+=&(&self.vx*dx(&self.T));
        delta+=&(&self.vy*dy(&self.T));
        delta/=H*H;
        return delta
    }
    fn step_phi(&self)->Array2<f64> {
        let mut delta = Array2::<f64>::zeros((NX,NY));
        delta+=&(laplace(&self.phi));
        delta+=&(&self.vx*dx(&self.phi));
        delta+=&(&self.vy*dy(&self.phi));
        delta/=H*H;
        delta+=&(self.params.R*dx(&self.T)/H);
        delta-=&(self.params.B*dx(&self.C)/H);
        return delta
    }
    fn step_psi(&mut self){
        poisson_relax(&self.phi,&mut self.psi,H)
    }

    fn next_step(&mut self,_dt:f64){
        self.phi+= &(self.step_phi()*_dt);
        self.step_psi();
        self.vx=-dy(&self.psi);
        self.vy=dx(&self.psi);
        self.T+=&(self.step_t()*_dt);
        self.C+=&(self.solverConc.solve(&self.psi,&self.T,&self.C,&self.params.paramsC)*_dt);
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
        psi_l:(system.psi[[NX/4,NY/2]]),
        nu:nu,
   } 
}


#[derive(Serialize)]
pub struct Row{
    t:f64,
    psi_m:f64,
    psi_l:f64,
    nu:f64,
}


fn main() {
    println!("Hello, world!");
    //let mut dela =Array2::<f64>::zeros((NX,NY));
    let contents = fs::read_to_string("config.toml").unwrap() ;
    let params:Params = toml::from_str(&contents).unwrap();
    println!("{:?}",params.R);

    let mut system = System{
        psi:Array2::<f64>::zeros((NX,NY)),
        phi:Array2::<f64>::zeros((NX,NY)),
        T:Array2::<f64>::zeros((NX,NY)),
        C:Array2::<f64>::zeros((NX,NY)),
        vx:Array2::<f64>::zeros((NX,NY)),
        vy:Array2::<f64>::zeros((NX,NY)),
        params:params,
        solverConc:SolverConc::new(),
    };
    //set initial 
    system.phi[[NX/2,NY/2]]=1.0;
    for i in 0..NX{
        for j in 0..NY{
            system.T[[i,j]] = 1.0-(j as f64)*H;
            system.C[[i,j]] = 30.5-(j as f64)*H;
        }
    }
    // log open
    let mut wtr = Writer::from_path("foo.csv").unwrap();

    let mut time= 0.0;
    let time_step = (params.time/DT)as usize;
    let start = Instant::now();
    for _x in 0..time_step{
        system.next_step(DT);
        time+=DT;
/*
        if _x % 1000==0{

            let res =  log_params(time,&system);
            wtr.serialize(res).unwrap();
            //write_mat(&system.C,String::from(format!("res/c_{:05}",_x/1000)),H);
        }
*/
    }
    println!("==============");

     let duration = start.elapsed();
    println!("Time elapsed in expensive_function() is: {:?}", duration);

    write_mat(&system.phi,String::from("phi"),H);
    write_mat(&system.psi,String::from("psi"),H);
    write_mat(&system.T,String::from("t"),H);
    write_mat(&system.C,String::from("c"),H);
}
