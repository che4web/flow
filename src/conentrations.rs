pub struct Concentration {
    pub f: Array2<f64>,
    qew: Array2<f64>,
    qsn: Array2<f64>,
    vx: Array2<f64>,
    vy: Array2<f64>,
    pub params: ParamsConc,
}

impl Concentration {
    pub fn new(nx: usize, ny: usize, params: ParamsConc) -> Self {
        return Self {
            f: Array2::<f64>::zeros((nx, ny)),
            qew: Array2::<f64>::zeros((nx, ny)),
            qsn: Array2::<f64>::zeros((nx, ny)),
            vx: Array2::<f64>::zeros((nx, ny)),
            vy: Array2::<f64>::zeros((nx, ny)),
            params: params,
        };
    }
    fn calc_v(&mut self, psi: &Psi) {
        unsafe {
            for i in 1..NX {
                let mut k = 0;
                self.vx[[i, k]] = psi.f[[i, k + 1]] + psi.f[[i - 1, k + 1]];

                for k in 1..NY - 1 {
                    self.vx[[i, k]] = psi.f.dy((i, k)) + psi.f.dy((i - 1, k));
                }
                k = NY - 1;
                self.vx[[i, k]] = -psi.f[[i, k - 1]] - psi.f[[i - 1, k - 1]];
            }
            self.vx /= 2.0;

            for i in 1..NX - 1 {
                for k in 1..NY {
                    self.vy[[i, k]] = psi.dx((i, k)) + psi.dx((i, k - 1));
                }
            }

            for k in 1..NY {
                {
                    let i = 0;
                    self.vy[[i, k]] = psi.f[[i + 1, k]] + psi.f[[i + 1, k - 1]];
                }
                let i = NX - 1;
                self.vy[[i, k]] = -psi.f[[i - 1, k]] - psi.f[[i - 1, k - 1]];
            }
            self.vy /= 2.0;
           
        }
    }
    pub fn step(&mut self, psi: &Psi,temp:&Temperatura, dt: f64) {
        let le = self.params.le;
        let sor = self.params.sor;
        let l_sed= 30.0/H;
       //println!("{:?}",le/l);
        unsafe {
            self.calc_v(&psi);

            for i in 1..NX {
                for k in 0..NY {
                    let mut tmp = (self.vx[[i, k]] ) * self.mx_b((i, k));
                    tmp += -le * self.dx_b((i, k));
                    tmp += -le * sor*temp.dx_b((i, k))* self.my_b((i, k));
                    self.qew[[i, k]] = tmp / H;
                }
            }
            for i in 0..NX {
                for k in 1..NY {
                    let mut tmp = (-self.vy[[i, k]]-le/l_sed) * self.my_b((i, k));
                    tmp += -le * (self.dy_b((i, k)));
                    tmp += -le * sor*temp.dy_b((i, k)) *self.my_b((i, k));
                    self.qsn[[i, k]] = tmp / H;
                }
            }

            for i in 1..NX - 1 {
                for k in 1..NY - 1 {
                    let q = -(self.qew.dx_f((i, k)) + self.qsn.dy_f((i, k)));
                    self.f[[i, k]] += dt / H * q;
                }
            }
            for i in 1..NX - 1 {
                let k = 0;
                self.f[[i, k]] += dt / H * (-self.qew.dx_f((i, k)) - 2.0 * self.qsn[[i, k + 1]]);

                let k = NY - 1;
                self.f[[i, k]] += dt / H * (-self.qew.dx_f((i, k)) + 2.0 * self.qsn[[i, k]]);
            }
            if PEREODIC{
                for k in 0..NY {
                    self.f[[0, k]] = self.f[[NX - 2, k]];
                    self.f[[NX - 1, k]] = self.f[[1, k]];
                }
            }else{
                for k in 1..NY - 1 {
                    let i= 0;
                    self.f[[i, k]] += dt / H * (-self.qsn.dy_f((i, k)) - 2.0 * self.qew[[i+1, k]]);
    
                    let i = NX - 1;
                    self.f[[i, k]] += dt / H * (-self.qsn.dy_f((i, k)) +2.0 * self.qew[[i, k]]);
                }
                self.f[[0, 0]] += dt / H * 2.0 *( -self.qew[[1, 0]]- self.qsn[[0, 1]]);
                self.f[[0, NY-1]] +=  dt / H * 2.0 *( -self.qew[[1, NY-1]]+ self.qsn[[0, NY-1]]);
                self.f[[NX-1, 0]] += dt / H * 2.0 *( self.qew[[NX-1, 0]]- self.qsn[[NX-1, 1]]);
                self.f[[NX-1, NY-1]] += dt / H * 2.0 *( self.qew[[NX-1, NY-1]]+ self.qsn[[NX-1, NY-1]]);
            }
        }
    }
}


