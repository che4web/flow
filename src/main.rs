const NY:usize=65;
const NX:usize=(NY-1)*2+1;
const L:f64=2.0;
const H:f64=L/((NX-1) as f64);
const DT:f64=H*H/5.0/10.;
const PEREODIC:bool = true;
extern crate queues;
use queues::*;


//pub mod io;
pub mod finit_diff;
pub mod io;
//use io::{writeMat};


use ndarray::{Array2};
use ndarray_stats::QuantileExt;

use csv::Writer;
use serde::{Serialize,Deserialize};


use finit_diff::{dx,dy,laplace,poisson_relax,dx_f,dy_f,dx_b,dy_b,mean_cx,mean_cy,dx_mut,dy_mut,laplace_mut};
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
    params_c:ParamsConc,

}
trait Field{
    fn get_f(&self)-> &Array2<f64>;
    unsafe fn dx(&self,index:(usize,usize))->f64{
        return (*self.get_f().uget((index.0+1,index.1))-*self.get_f().uget((index.0-1,index.1)))/2.0
    }

    unsafe fn dy(&self,i:(usize,usize))->f64{
        return (*self.get_f().uget((i.0,i.1+1))-*self.get_f().uget((i.0,i.1-1)))/2.0
    }

    unsafe fn lap(&self,i:(usize,usize))->f64{
        return *self.get_f().uget((i.0+1,i.1  ))+
               *self.get_f().uget((i.0-1,i.1  ))+
               *self.get_f().uget((i.0  ,i.1+1))+
               *self.get_f().uget((i.0  ,i.1-1))-
               *self.get_f().uget(i)*4.0;
    }
}

struct System{
    temp:Temperatura,
    psi:Psi,
    phi:Phi,
    conc:Concentration,
    params:Params,
}
struct Concentration{
    f:Array2<f64>,
    delta:Array2<f64>,
    qew:Array2<f64>,
    qsn:Array2<f64>,
    params:ParamsConc,
}

struct Temperatura{
    f:Array2<f64>,
    delta:Array2<f64>,
}
struct Phi{
    f:Array2<f64>,
    delta:Array2<f64>,
    params:Params,
}

struct Psi{
    f:Array2<f64>,
}



impl Field for Temperatura{
    fn get_f(&self)-> &Array2<f64>{
        return &self.f
    }
}
impl Field for Phi{
    fn get_f(&self)-> &Array2<f64>{
        return &self.f
    }
}
impl Field for Psi{
    fn get_f(&self)-> &Array2<f64>{
        return &self.f
    }
}
impl Field for Concentration{
    fn get_f(&self)-> &Array2<f64>{
        return &self.f
    }
}





impl Concentration{
    fn new(nx:usize,ny:usize,params:ParamsConc)-> Self{
        return Self{
        f:Array2::<f64>::zeros((nx,ny)),
        delta:Array2::<f64>::zeros((nx,ny)),
        qsn:Array2::<f64>::zeros((nx,ny)),
        qew:Array2::<f64>::zeros((nx,ny)),
        params:params
        }
    }
    fn step(&mut self,psi:&Array2<f64>,t:&Array2<f64>,dt:f64) {

        let le = self.params.Le;
        let l = self.params.l;

        let a = vec![
            0.25/H*0.5,
            -0.25/H*0.5,
            le/l*0.5,
            -le/H,
            -le/H,
            dt/H,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0
        ];
        for i in 1..NX{
            for k in 1..NY-1{
                let c_local=self.f[[i,k]]+self.f[[i-1,k]];
                let  mut tmp=a[0]*(psi[[i,k+1]]-psi[[i,k-1]]+psi[[i-1,k+1]]-psi[[i-1,k-1]])*c_local;
                tmp+=a[3]*(self.f[[i,k]]-self.f[[i-1,k]]);
                self.qew[[i,k]] = tmp;
            }
        }
        for i in 1..NX-1{
            for k in 1..NY{
                let c_local=self.f[[i,k]]+self.f[[i,k-1]];
                let  mut tmp=a[1]*(psi[[i+1,k]]-psi[[i-1,k]]+psi[[i+1,k-1]]-psi[[i-1,k-1]])*c_local;
                tmp+=a[4]*(self.f[[i,k]]-self.f[[i,k-1]]);
                self.qsn[[i,k]] = tmp;
            }
        }

        for i in 1..NX{
            let mut  k=0;
            let mut c_local=self.f[[i,k]]+self.f[[i-1,k]];
            let  mut tmp=a[0]*(psi[[i,k+1]]+psi[[i-1,k+1]])*c_local*2.0;
            tmp+=a[3]*(self.f[[i,k]]-self.f[[i-1,k]]);
            self.qew[[i,k]] = tmp;

            k=NY-1;
            c_local=self.f[[i,k]]+self.f[[i-1,k]];
            tmp=a[0]*(-psi[[i,k-1]]-psi[[i-1,k-1]])*c_local*2.0;
            tmp+=a[3]*(self.f[[i,k]]-self.f[[i-1,k]]);
            self.qew[[i,k]] = tmp;

        }
        for k in 1..NY{
            let mut  i=0;
            let mut c_local=self.f[[i,k]]+self.f[[i,k-1]];
            let  mut tmp=a[1]*(psi[[i+1,k]]+psi[[i+1,k-1]])*c_local*2.0;
            tmp+=a[3]*(self.f[[i,k]]-self.f[[i,k-1]]);
            self.qsn[[i,k]] = tmp;

            i=NX-1;
            c_local=self.f[[i,k]]+self.f[[i,k-1]];
            tmp=a[1]*(-psi[[i-1,k]]-psi[[i-1,k-1]])*c_local*2.0;
            tmp+=a[3]*(self.f[[i,k]]-self.f[[i,k-1]]);
            self.qsn[[i,k]] = tmp;
        }
        for i in 1..NX-1{
            for k in 1..NY-1{
                let q = self.qew[[i,k]]-self.qew[[i+1,k]]+self.qsn[[i,k]]-self.qsn[[i,k+1]];
                self.f[[i,k]]=self.f[[i,k]]+a[5]*q;
            }
        }
        for i in 1..NX-1{
            let k=0;
            self.f[[i,k]]=self.f[[i,k]]+a[5]*((self.qew[[i,k]]-self.qew[[i+1,k]])-2.0*self.qsn[[i,k+1]]);

            let k=NY-1;
            self.f[[i,k]]=self.f[[i,k]]+a[5]*((self.qew[[i,k]]-self.qew[[i+1,k]])+2.0*self.qsn[[i,k]]);
        }
        for k in 0..NY{
            self.f[[0,k]]=self.f[[NX-2,k]];

            self.f[[NX-1,k]]=self.f[[1,k]];
            //delta[[0,k]]=c[[NX-2,k]];
            //delta[[NX-1,k]]=c[[1,k]];
        }

    }

}



impl Temperatura{
    fn new(nx:usize,ny:usize)-> Self{
        return Self{
        f:Array2::<f64>::zeros((nx,ny)),
        delta:Array2::<f64>::zeros((nx,ny)),
        }
    }
    fn step(&mut self,psi:&Psi,dt:f64){
        unsafe{
        for i in 1..NX-1{
            for j in 1..NY-1{
                let mut tmp = self.lap((i,j));
                tmp-=psi.dy((i,j))*(self.dx((i,j)));
                tmp+=psi.dx((i,j))*(self.dy((i,j)));
                tmp/=H*H;
                *self.delta.uget_mut((i,j))=tmp*dt;
            }
        }
        }
        self.f+=&self.delta;
    }
}
impl Phi{
    fn new(nx:usize,ny:usize,params:Params)-> Self{
        return Self{
        f:Array2::<f64>::zeros((nx,ny)),
        delta:Array2::<f64>::zeros((nx,ny)),
        params:params,
        }
    }
    fn step(&mut self,psi:&Psi,t:&Temperatura,conc:&Concentration,dt:f64) {
        //let mut delta = Array2::<f64>::zeros((NX,NY));
        //let shape= delta.dim();
        unsafe{
        for i in 1..NX-1{
            for j in 1..NY-1{
                let mut tmp = self.lap((i,j))*self.params.P;
                tmp-=psi.dy((i,j))*(self.dx((i,j)));
                tmp+=psi.dx((i,j))*(self.dy((i,j)));
                tmp/=H*H;
                tmp+=self.params.P*self.params.R*(t.dx((i,j)))/H;
                tmp-=self.params.P*self.params.B*(conc.dx((i,j)))/H;
                *self.delta.uget_mut((i,j))=tmp*dt;
            
            }
        }
        }
        self.f+=&self.delta;
    }
}
impl Psi{
    fn new(nx:usize,ny:usize)-> Self{
        return Self{
        f:Array2::<f64>::zeros((nx,ny)),
        }
    }
    fn step(&mut self,phi:&Phi,) {
        poisson_relax(&phi.f,&mut self.f,H);
    }
}






impl System{


    fn next_step(&mut self,_dt:f64){
        self.phi.step(&self.psi,&self.temp,&self.conc,_dt);
        self.psi.step(&self.phi);
        self.temp.step(&self.psi,_dt);
        self.conc.step(&self.psi.f,&self.temp.f,_dt);
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
    println!("{:?}",params.R);

    let mut system = System{
        temp:Temperatura::new(NX,NY),
        phi:Phi::new(NX,NY,params),
        psi:Psi::new(NX,NY),
        conc:Concentration::new(NX,NY,params.params_c),
        params:params,
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
