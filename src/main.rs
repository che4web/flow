const NY:usize=33;
const NX:usize=(NY-1)*2+1;
const L:f64=2.0;
const H:f64=L/((NX-1) as f64);
const DT:f64=H*H/5.0/10.;
const PEREODIC:bool = true;


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
struct System{
    temp:Field,
    psi:Field,
    phi:Field,
    conc:Array2<f64>,
    params:Params,
}
struct Field{
    f:Array2<f64>,
    dx:Array2<f64>,
    dy:Array2<f64>,
    lap:Array2<f64>,
    //buf:Array2<f64>,
}
impl Field{
    fn new(nx:usize,ny:usize)-> Self{
        return Self{
        f:Array2::<f64>::zeros((nx,ny)),
        dx:Array2::<f64>::zeros((nx,ny)),
        dy:Array2::<f64>::zeros((nx,ny)),
        lap:Array2::<f64>::zeros((nx,ny)),
      //  buf:Array2::<f64>::zeros((nx,ny)),
        }
    }
    unsafe fn dx_f(&self,index:(usize,usize))->f64{
        return (*self.f.uget((index.0+1,index.1))-*self.f.uget((index.0-1,index.1)))/2.0
    }

    unsafe fn dy_f(&self,i:(usize,usize))->f64{
        return (*self.f.uget((i.0,i.1+1))-*self.f.uget((i.0,i.1-1)))/2.0
    }

    unsafe fn lap_f(&self,i:(usize,usize))->f64{
        return *self.f.uget((i.0+1,i.1  ))+
               *self.f.uget((i.0-1,i.1  ))+
               *self.f.uget((i.0  ,i.1+1))+
               *self.f.uget((i.0  ,i.1-1))-
               *self.f.uget(i)*4.0;
    }
    fn diff(&mut self){
        dx_mut(&self.f,&mut self.dx);
        dy_mut(&self.f,&mut self.dy);
        //laplace_mut(&self.f,&mut self.lap);
    }
    
}
impl System{
    fn step_t(&mut self,dt:f64){
        //let mut delta = Array2::<f64>::zeros((NX,NY));

        unsafe{
        for i in 1..NX-1{
            for j in 1..NY-1{
                let mut tmp = self.temp.lap_f((i,j));
                tmp-=*self.psi.dy.uget((i,j))*(self.temp.dx_f((i,j)));
                tmp+=*self.psi.dx.uget((i,j))*(self.temp.dy_f((i,j)));
                tmp/=H*H;
                *self.temp.lap.uget_mut((i,j))=tmp*dt;
            }
        }
        }
        self.temp.f+=&self.temp.lap;
      //  delta+=&self.temp.lap;
       // delta-=&(&self.psi.dy*&self.temp.dx);
        //delta+=&(&self.psi.dx*&self.temp.dy);
        //delta/=H*H;
    }
    fn step_phi(&mut self,dt:f64) {
        //let mut delta = Array2::<f64>::zeros((NX,NY));
        //let shape= delta.dim();
        unsafe{
        for i in 1..NX-1{
            for j in 1..NY-1{
                let dc = (*self.conc.uget((i+1,j))-*self.conc.uget((i-1,j)))/2.0;
                let mut tmp = self.phi.lap_f((i,j))*self.params.P;
                tmp-=*self.psi.dy.uget((i,j))*(self.phi.dx_f((i,j)));
                tmp+=*self.psi.dx.uget((i,j))*(self.phi.dy_f((i,j)));
                tmp/=H*H;
                tmp+=self.params.P*self.params.R*(self.temp.dx_f((i,j)))/H;
                tmp-=self.params.P*self.params.B*(dc)/H;
                *self.phi.lap.uget_mut((i,j))=tmp*dt;
            
            }
        }
        }
        //delta+=&self.phi.lap;
        //delta*=self.params.P;
        //delta-=&(&self.psi.dy*&self.phi.dx);
        //delta+=&(&self.psi.dx*&self.phi.dy);
        //delta/=H*H;
        //delta+=&(self.params.P*self.params.R*&self.temp.dx/H);
        //delta*=dt;
        self.phi.f+=&self.phi.lap;
   //     delta-=&(self.params.P*self.params.B*dx(&self.conc)/H);
       // return delta
    }
    fn step_c(&self,psi:&Array2<f64>,t:&Array2<f64>,c:&Array2<f64>,params:&ParamsConc,dt:f64)->Array2<f64> {

        let le = params.Le;
        let l = params.l;

        let mut delta = Array2::<f64>::zeros((NX,NY));
        let mut qew = Array2::<f64>::zeros((NX,NY));
        let mut qsn = Array2::<f64>::zeros((NX,NY));
        
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
                let c_local=c[[i,k]]+c[[i-1,k]];
                let  mut tmp=a[0]*(psi[[i,k+1]]-psi[[i,k-1]]+psi[[i-1,k+1]]-psi[[i-1,k-1]])*c_local;
                tmp+=a[3]*(c[[i,k]]-c[[i-1,k]]);
                qew[[i,k]] = tmp;
            }
        }
        for i in 1..NX-1{
            for k in 1..NY{
                let c_local=c[[i,k]]+c[[i,k-1]];
                let  mut tmp=a[1]*(psi[[i+1,k]]-psi[[i-1,k]]+psi[[i+1,k-1]]-psi[[i-1,k-1]])*c_local;
                tmp+=a[4]*(c[[i,k]]-c[[i,k-1]]);
                qsn[[i,k]] = tmp;
            }
        }

        for i in 1..NX{
            let mut  k=0;
            let mut c_local=c[[i,k]]+c[[i-1,k]];
            let  mut tmp=a[0]*(psi[[i,k+1]]+psi[[i-1,k+1]])*c_local*2.0;
            tmp+=a[3]*(c[[i,k]]-c[[i-1,k]]);
            qew[[i,k]] = tmp;

            k=NY-1;
            c_local=c[[i,k]]+c[[i-1,k]];
            tmp=a[0]*(-psi[[i,k-1]]-psi[[i-1,k-1]])*c_local*2.0;
            tmp+=a[3]*(c[[i,k]]-c[[i-1,k]]);
            qew[[i,k]] = tmp;

        }
        for k in 1..NY{
            let mut  i=0;
            let mut c_local=c[[i,k]]+c[[i,k-1]];
            let  mut tmp=a[1]*(psi[[i+1,k]]+psi[[i+1,k-1]])*c_local*2.0;
            tmp+=a[3]*(c[[i,k]]-c[[i,k-1]]);
            qsn[[i,k]] = tmp;

            i=NX-1;
            c_local=c[[i,k]]+c[[i,k-1]];
            tmp=a[1]*(-psi[[i-1,k]]-psi[[i-1,k-1]])*c_local*2.0;
            tmp+=a[3]*(c[[i,k]]-c[[i,k-1]]);
            qsn[[i,k]] = tmp;
        }
        for i in 1..NX-1{
            for k in 1..NY-1{
                let q = qew[[i,k]]-qew[[i+1,k]]+qsn[[i,k]]-qsn[[i,k+1]];
                delta[[i,k]]=c[[i,k]]+a[5]*q;
            }
        }
        for i in 1..NX-1{
            let k=0;
            delta[[i,k]]=c[[i,k]]+a[5]*((qew[[i,k]]-qew[[i+1,k]])-2.0*qsn[[i,k+1]]);

            let k=NY-1;
            delta[[i,k]]=c[[i,k]]+a[5]*((qew[[i,k]]-qew[[i+1,k]])+2.0*qsn[[i,k]]);
        }
        for k in 0..NY{
            delta[[0,k]]=delta[[NX-2,k]];
            delta[[NX-1,k]]=delta[[1,k]];
            //delta[[0,k]]=c[[NX-2,k]];
            //delta[[NX-1,k]]=c[[1,k]];
        }



        return delta
    }

 
    fn step_psi(&mut self){
        poisson_relax(&self.phi.f,&mut self.psi.f,H)
    }

    fn next_step(&mut self,_dt:f64){
        self.step_phi(_dt);
        //self.phi.diff();

        self.step_psi();
        self.psi.diff();

        self.step_t(_dt);
        //self.temp.diff();

        self.conc=self.step_c(&self.psi.f,&self.temp.f,&self.conc,&self.params.params_c,_dt);
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
        temp:Field::new(NX,NY),
        phi:Field::new(NX,NY),
        psi:Field::new(NX,NY),
        conc:Array2::<f64>::zeros((NX,NY)),
        params:params,
    };
    //set initial 
    system.phi.f[[NX/2,NY/2]]=1.0;
    for i in 0..NX{
        for j in 0..NY{
            system.temp.f[[i,j]] = 1.0-(j as f64)*H;
            system.conc[[i,j]] = 30.5-(j as f64)*H;
        }
    }
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
            write_mat(&system.conc,String::from(format!("res/c_{:05}",time)),H);
        }
    }
    println!("==============");

     let duration = start.elapsed();
    println!("Time elapsed in expensive_function() is: {:?}", duration);

    write_mat(&system.phi.f,String::from("phi"),H);
    write_mat(&system.psi.f,String::from("psi"),H);
    write_mat(&system.temp.f,String::from("t"),H);
    write_mat(&system.conc,String::from("c"),H);
}
