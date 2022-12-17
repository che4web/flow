use std::fs::File;
use std::io::prelude::*;
use ndarray::{ Array2};
use crate::{FMatrix};
pub fn write_mat(phi:&FMatrix,name:String,h:f64){
    let mut file = File::create(name+".txt").unwrap();
    for j in  0..phi.shape().1{
        for i in  0..phi.shape().0{
            writeln!(file, "{:+1.6} {:+1.6} {:+1.6}",(i as f64)*h,(j as f64)*h,phi[(i,j)]).unwrap();
        }
    }
}


