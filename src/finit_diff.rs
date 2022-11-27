use ndarray::{ Array2};

pub fn poisson_relax(phi: &Array2<f64>,psi0: &Array2<f64>,H:f64)->Array2<f64>{
    const OMEGA:f64 =1.8; 
    //let mut delta = Array2::<f64>::zeros((N,N));
    let mut psi= (psi0).clone();
    let shape= psi0.shape();
    for _k in  0..1000{
        //let lap = laplace(&psi_new);
        let mut nev =0.0;
        for i in 1..shape[0]-1{
            for j in 1..shape[1]-1{
                let tmp =-OMEGA*psi[[i,j]]+ OMEGA/4.0*(psi[[i+1,j]] +psi[[i-1,j]]+psi[[i,j+1]] +psi[[i,j-1]]+H*H*phi[[i,j]]);

                psi[[i,j]] +=tmp;
                if nev<tmp.abs(){
                    nev =tmp.abs();
                }
            }
        }
        //println!("{:?}| {:?} | {:?} {:?}",nev,k,i_max,j_max);
        if nev < (1e-7 as f64){
            break;
        }

    }
    return psi
}


pub fn laplace(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
            delta[[j,i]]= arr[[j+1,i]]+arr[[j-1,i]]+arr[[j,i+1]]+arr[[j,i-1]]-4.0*arr[[j,i]];
        }
    }
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
pub fn dy(arr: &Array2<f64>)->Array2<f64> {
    let shape = arr.dim();
    let mut delta = Array2::<f64>::zeros(shape);
    for j in 1..shape.0-1{
        for i in 1..shape.1-1{
            delta[[j,i]]=( arr[[j,i+1]]-arr[[j,i-1]])/(2.0);
        }
    }
    return delta
}


