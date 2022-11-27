use std::fs::File;
use std::io::prelude::*;
use ndarray::{ Array2};
pub fn writeMat(phi: &Array2<f64>,name:String,H:f64){
    let mut file = File::create(name+".txt").unwrap();
    for j in  0..phi.shape()[1]{
        for i in  0..phi.shape()[0]{
            writeln!(file, "{:+1.6} {:+1.6} {:+1.6}",(i as f64)*H,(j as f64)*H,phi[[i,j]]).unwrap();
        }
    }
}


