const NY:usize=65;
const NX:usize=(NY-1)*2+1;
const L:f64=2.0;
const H:f64=L/((NX-1) as f64);
const DT:f64=H*H/5.0/10.;
const PEREODIC:bool = true;

use ndarray_stats::QuantileExt;
//pub mod io;
pub mod finit_diff;
pub mod io;
pub mod field;
//use io::{writeMat};
pub mod equtions;
use equtions::{Concentration,Temperatura,Psi,Phi,Params,ParamsConc};
use csv::Writer;
use serde::{Serialize,Deserialize};



use io::write_mat;
use std::fs;

use std::time::{Instant};
use toml;


struct System{
    temp:Temperatura,
    psi:Psi,
    phi:Phi,
    conc:Concentration,
}



impl System{


    fn next_step(&mut self,_dt:f64){
        self.phi.step(&self.psi,&self.temp,&self.conc,_dt);
        self.psi.step(&self.phi);
        self.temp.step(&self.psi,_dt);
        self.conc.step(&self.psi,&self.temp.f,_dt);
        //self.temp.diff();

        //self.conc=self.step_c(&self.psi.f,&self.temp.f,&self.conc,&self.params.params_c,_dt);
        self.boundary_condition();
    }
    fn boundary_condition(&mut self){
        if PEREODIC{

            for i in 0..NY{
                self.psi.f[[0,i]]=self.psi.f[[NX-2,i]];
                self.psi.f[[NX-1,i]]=self.psi.f[[1,i]];

                self.phi.f[[0,i]]=self.phi.f[[NX-2,i]];
                self.phi.f[[NX-1,i]]=self.phi.f[[1,i]];

                self.temp.f[[0,i]]=self.temp.f[[NX-2,i]];
                self.temp.f[[NX-1,i]]=self.temp.f[[1,i]];

                //self.C[[0,i]]=self.C[[NX-2,i]];
                //self.C[[NX-1,i]]=self.C[[2,i]];
            }
        }else{
            for i in 1..NY-1{
                self.phi.f[[0,i]]=-self.psi.f[[1,i]]/(H*H)*2.0;
                self.phi.f[[NX-1,i]]=-self.psi.f[[NX-2,i]]/(H*H)*2.0;
            }
        }

        for i in 1..NX-1{
            self.phi.f[[i,0]]=-self.psi.f[[i,1]]/(H*H)*2.0;
            self.phi.f[[i,NY-1]]=-self.psi.f[[i,NY-2]]/(H*H)*2.0;
        }

    }
}
fn log_params(time:f64,system: &System)->Row{
    let mut nu  =0.0;
    let shape = system.temp.f.dim();
    for i in 1..shape.0-1{
        nu+= system.temp.f[[i,shape.1-1]]-system.temp.f[[i,shape.1-2]]
    }
    nu/=H*(shape.0 as f64)-2.0;

   return Row{
        t:time,
        psi_m:*(system.psi.f.max().unwrap()),
        psi_l:(system.psi.f[[NX/4,NY/2]]),
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
    println!("{:?}",params.rel);

    let mut system = System{
        temp:Temperatura::new(NX,NY),
        phi:Phi::new(NX,NY,params),
        psi:Psi::new(NX,NY),
        conc:Concentration::new(NX,NY,params.params_c),
    };
    //set initial 
    system.phi.f[[NX/2,NY/2]]=1.0;
    for i in 0..NX{
        for j in 0..NY{
            system.temp.f[[i,j]] = 1.0-(j as f64)*H;
            system.conc.f[[i,j]] = 30.5-(j as f64)*H;
        }
    }
    system.phi.f[[NX/2,NY/2]] = 1.0;
    // log open
    let mut wtr = Writer::from_path("foo.csv").unwrap();

    let mut time= 0.0;
    let start = Instant::now();
    let mut  time_i=0;
    while time < params.time{
        for _i in 0..1000{
            system.next_step(DT);
            time+=DT;
        }
        time_i+=1;
        let res =  log_params(time,&system);
        wtr.serialize(res).unwrap();
        if time_i %10 == 0{
            let duration = start.elapsed();
            println!("{:0.2},{:?}",time/params.time,duration) ;
            write_mat(&system.conc.f,String::from(format!("res/c_{:05}",time)),H);
        }
    }
    println!("==============");

     let duration = start.elapsed();
    println!("Time elapsed in expensive_function() is: {:?}", duration);

    write_mat(&system.phi.f,String::from("phi"),H);
    write_mat(&system.psi.f,String::from("psi"),H);
    write_mat(&system.temp.f,String::from("t"),H);
    write_mat(&system.conc.f,String::from("c"),H);
}
