//! Module to simply create, add and multiply complex numbers in cartesian form, both the real and
//! imaginary part of these complexes will be f64

pub struct Complex {
    pub a: f64,
    pub b: f64,
}

impl Complex {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            a,
            b,
        }
    }

    pub fn add(&mut self, other: &Complex) {
        self.a += other.a;
        self.b += other.b;
    }

    pub fn mult(&mut self, other: &Complex) {
        self.a = self.a*other.a - self.b*other.b;
        self.b = self.a*other.b + self.b*other.a;
    }

    ///multiplying by an int and returning the results
    pub fn mult_u64(&self, other: u64) -> Self {
        Self{
            a: self.a*other as f64,
            b: self.b*other as f64
        }
    }

    ///adds a u64 "in place"
    pub fn add_u64(&mut self, n: u64) {
        self.a += n as f64;
    }

    ///substracts a u64 "in place"
    pub fn sub_u64(&mut self, n: u64) {
        self.a -= n as f64;
    }

    //function that takes a "center of mass" as self and a "rotation complex", and checks if the 
    //kmer represented by this mass is a decycler or not
    pub fn check_decycler(&self, rot: &Complex) -> bool {
        let epsilon: f64 = 0.00001;
        if self.b > epsilon {
            if self.a*rot.b + self.b*rot.a < -epsilon {
                return true
            }
        }
        false
    }
}
