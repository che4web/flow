use ndarray::{Array2};
pub trait Field{
    fn get_f(&self)-> &Array2<f64>;
    unsafe fn mx_b(&self,index:(usize,usize))->f64{
        return (*self.get_f().uget((index.0,index.1))+*self.get_f().uget((index.0-1,index.1)))/2.0
    }

    unsafe fn my_b(&self,index:(usize,usize))->f64{
        return (*self.get_f().uget((index.0,index.1))+*self.get_f().uget((index.0,index.1-1)))/2.0
    }

    unsafe fn dx_f(&self,index:(usize,usize))->f64{
        return *self.get_f().uget((index.0+1,index.1))-*self.get_f().uget((index.0,index.1))
    }

    unsafe fn dx_b(&self,index:(usize,usize))->f64{
        return *self.get_f().uget((index.0,index.1))-*self.get_f().uget((index.0-1,index.1))
    }

    unsafe fn dy_b(&self,index:(usize,usize))->f64{
        return *self.get_f().uget((index.0,index.1))-*self.get_f().uget((index.0,index.1-1))
    }

    unsafe fn dy_f(&self,index:(usize,usize))->f64{
        return *self.get_f().uget((index.0,index.1+1))-*self.get_f().uget((index.0,index.1))
    }

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


pub unsafe fn dx_b(f:&Array2<f64>,index:(usize,usize))->f64{
    return *f.uget((index.0,index.1))-*f.uget((index.0-1,index.1))
}

pub unsafe fn dy_b(f:&Array2<f64>,index:(usize,usize))->f64{
    return *f.uget((index.0,index.1))-*f.uget((index.0,index.1-1))
}