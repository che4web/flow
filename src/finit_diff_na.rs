use ndarray::{Array2};
use crate::{FMatrix};
const PEREODIC:bool  = crate::PEREODIC; 

pub fn poisson_relax(phi: &FMatrix,psi: &mut FMatrix,h:f64){
    const OMEGA:f64 =1.92; 
    //let mut delta = Array2::<f64>::zeros((N,N));
    let shape= psi.shape();
    let nx = shape.0;
    let ny = shape.1;
    let h2=h*h;
    for _k in  0..1000{
        //let lap = laplace(&psi_new);
        let mut nev =0.0;
        for i in 1..shape.0-1{
            for j in 1..shape.1-1{
                let tmp =-OMEGA*psi[(i,j)]+ OMEGA/4.0*(psi[(i+1,j)] +psi[(i-1,j)]+psi[(i,j+1)] +psi[(i,j-1)]+h2*phi[(i,j)]);

                psi[(i,j)] +=tmp;
                if nev<tmp.abs(){
                    nev =tmp.abs();
                }
            }
        }
        if PEREODIC{
            for i in 0..ny{
                psi[(0,i)]=psi[(nx-2,i)];
                psi[(nx-1,i)]=psi[(2,i)];
            }
        }

        //println!("{:?}| {:?} | {:?} {:?}",nev,k,i_max,j_max);
        if nev < (1e-5 as f64){
            break;
        }

    }
}


pub fn laplace(v: &Array2<f64>)->Array2<f64> {
    let shape = v.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
           delta[[j,i]]= v[[j+1,i]]+v[[j-1,i]]+v[[j,i+1]]+v[[j,i-1]]-4.0*v[[j,i]];
       }
    }
    /*
    let tmp = -4. * &v.slice(s![1..-1, 1..-1])
        + v.slice(s![ ..-2, 1..-1])
        + v.slice(s![1..-1,  ..-2])
        + v.slice(s![1..-1, 2..  ])
        + v.slice(s![2..  , 1..-1]);
   delta.slice_mut(s![1..-1, 1..-1]).assign(&tmp);
   */ 
   return delta
}

pub fn dx(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
            delta[[j,i]]=( arr[[j+1,i]]-arr[[j-1,i]])/(2.0);
        }
    }
    return delta
}
pub fn dx_f(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 0..shape.0-1{
        for i in 0..shape.1{
            delta[[j,i]]= arr[[j+1,i]]-arr[[j,i]];
        }
    }
    return delta
}
pub fn dy_f(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 0..shape.0{
        for i in 0..shape.1-1{
            delta[[j,i]]= arr[[j,i+1]]-arr[[j,i]];
        }
    }
    return delta
}

pub fn dx_b(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 1..shape.0{
        for i in 0..shape.1{
            delta[[j,i]]= arr[[j,i]]-arr[[j-1,i]];
        }
    }
    return delta
}
pub fn dy_b(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 0..shape.0{
        for i in 1..shape.1{
            delta[[j,i]]= arr[[j,i]]-arr[[j,i-1]];
        }
    }
    return delta
}
pub fn mean_cx(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 1..shape.0{
        for i in 0..shape.1{
            delta[[j,i]]= (arr[[j,i]]+arr[[j-1,i]])/2.0;
        }
    }
    return delta
}
pub fn mean_cy(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 0..shape.0{
        for i in 1..shape.1{
            delta[[j,i]]= (arr[[j,i]]+arr[[j,i-1]])/2.0;
        }
    }
    return delta
}




pub fn dy(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
           delta[[j,i]]=( arr[[j,i+1]]-arr[[j,i-1]])/(2.0);
       }
   }

    //let arr = array![0, 1, 2, 3];
    //assert_eq!(arr.slice(s![1..3;-1]), array![2, 1]);
//     let tmp=(&arr.slice(s![1..-1, 2..])-&arr.slice(s![1..-1, 0..-2]))/2.0;

  //   delta.slice_mut(s![1..-1, 1..-1]).assign(&tmp);
    //println!("{:?}",tmp.shape());
    //delta= &x.slice(s![1..]) - &x.slice(s![..-1]);
    return delta
}


pub fn dy_mut(arr: &FMatrix,delta:&mut FMatrix){
    let shape = arr.shape();
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
           delta[(j,i)]=( arr[(j,i+1)]-arr[(j,i-1)])/(2.0);
       }
   }
}
pub fn dx_mut(arr: &FMatrix,delta:&mut FMatrix){
    let shape = arr.shape();
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
           delta[(j,i)]=( arr[(j+1,i)]-arr[(j-1,i)])/(2.0);
       }
   }
}
pub fn laplace_mut(v: &FMatrix,delta:&mut FMatrix) {
    let shape = v.shape();
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
           delta[(j,i)]= v[(j+1,i)]+v[(j-1,i)]+v[(j,i+1)]+v[(j,i-1)]-4.0*v[(j,i)];
       }
    }
}
 









