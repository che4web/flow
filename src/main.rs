const NY:usize=65;
const NX:usize=(NY-1)*2+1;
const L:f64=2.0;
const H:f64=L/((NX-1) as f64);
const DT:f64=H*H/5.0/10.;
const pereodic:bool = true;


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
        delta*=self.params.P;
        delta-=&(dy(&self.psi)*dx(&self.phi));
        delta+=&(dx(&self.psi)*dy(&self.phi));
        delta/=H*H;
        delta+=&(self.params.P*self.params.R*dx(&self.T)/H);
        delta-=&(self.params.P*self.params.B*dx(&self.C)/H);
        return delta
    }
    fn step_c(&self,psi:&Array2<f64>,t:&Array2<f64>,c:&Array2<f64>,params:&ParamsConc)->Array2<f64> {

        let le = params.Le;
        let l = params.l;

        let mut delta = Array2::<f64>::zeros((NX,NY));
        let mut vx=-dy(&psi);
        let mut vy=dx(&psi);
        for i in 1..NX{
            vx[[i,0]] = -psi[[i,1]];
            vx[[i,NY-1]] = psi[[i,NY-2]];
        }

        for j in 1..NY{
            vy[[0,j]] = psi[[1,j]];
            vy[[NX-1,j]] = -psi[[NX-2,j]];
        }
        vx = mean_cx(&vx);
        vy = mean_cy(&vy);

        let cx = dx_b(&c)*le;
        let cy = dy_b(&c)*le;
        let tx = dx_b(&t)*0.0;
        let ty = dy_b(&t)*0.0;

        let mut c_loc= mean_cx(&c);

        let mut jx = Array2::<f64>::zeros((NX,NY));
        let mut jy = Array2::<f64>::zeros((NX,NY));
        jx+=&tx;
        jx+=&cx;
        jx+=&(&vx*&c_loc);

        c_loc= mean_cy(&c);
        jy+=&ty;
        jy+=&cy;
        jy+=&(&vy*&c_loc);
        jy+=&(le/l*&c_loc*H);

        let mut djx = dx_f(&jx);
        let mut djy = dy_f(&jy);
        for j in 1..NY{
            djx[[0,j]] = 2.0*jx[[1,j]];
            djx[[NX-1,j]] = -2.0*jx[[NX-1,j]];
            //djy[[NX-1,j]] = -0.0;
        }

        for i in 1..NX{
            djy[[i,0]] = 2.0*(jy[[i,1]]);
            djy[[i,NY-1]] = -2.0*(jy[[i,NY-1]]);
        }

        // djx = dx_f(&jx);
        // djy = dy_f(&jy);


        delta+=&(&djx+&djy);
        delta[[0,0]]=-2.0*(-jx[[1,0]]-jy[[0,1]]);
        delta[[0,NY-1]]=-2.0*(-jx[[1,NY-1]]+jy[[0,NY-1]]);
        delta[[NX-1,0]]=-2.0*(jx[[NX-1,0]]-jy[[NX-1,1]]);
        delta[[NX-1,NY-1]]=-2.0*(jx[[NX-1,NY-1]]+jy[[NX-1,NY-1]]);
        
        delta/=H*H;
        return delta;
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
        self.C+=&(self.step_c(&self.psi,&self.T,&self.C,&self.params.paramsC)*_dt);
        self.boundary_condition();
    }
    fn boundary_condition(&mut self){
        if pereodic{

            for i in 0..NY{
                self.psi[[0,i]]=self.psi[[NX-2,i]];
                self.psi[[NX-1,i]]=self.psi[[2,i]];

                self.phi[[0,i]]=self.phi[[NX-2,i]];
                self.phi[[NX-1,i]]=self.phi[[2,i]];

                self.T[[0,i]]=self.T[[NX-2,i]];
                self.T[[NX-1,i]]=self.T[[2,i]];

                self.C[[0,i]]=self.C[[NX-2,i]];
                self.C[[NX-1,i]]=self.C[[2,i]];
            }
        }else{
            for i in 1..NY-1{
                self.phi[[0,i]]=-self.psi[[1,i]]/(H*H)*2.0;
                self.phi[[NX-1,i]]=-self.psi[[NX-2,i]]/(H*H)*2.0;
            }
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
    let mut  time_i=0;
    while time < params.time{
        for i in 0..1000{
            system.next_step(DT);
            time+=DT;
        }
        time_i+=1;
        let res =  log_params(time,&system);
        wtr.serialize(res).unwrap();
        if time_i %100 == 0{
            write_mat(&system.C,String::from(format!("res/c_{:05}",time)),H);
        }
    }
    println!("==============");

     let duration = start.elapsed();
    println!("Time elapsed in expensive_function() is: {:?}", duration);

    write_mat(&system.phi,String::from("phi"),H);
    write_mat(&system.psi,String::from("psi"),H);
    write_mat(&system.T,String::from("t"),H);
    write_mat(&system.C,String::from("c"),H);
}
