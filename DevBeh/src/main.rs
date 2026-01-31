use rand::Rng;
use rand::thread_rng;
use rand::seq::SliceRandom;
use rand_distr::{Normal, Distribution};
use statrs::distribution::{Normal as nm};
use weighted_rand::builder::*;
use std::error::Error;
use std::fs::OpenOptions;
use rayon::prelude::*;
use std::env;
use std::fs::{self, File};
use std::path::Path;
use csv::Writer;
use std::io::BufWriter;
use std::io;

// R stuff
use std::process::{Command, Stdio};
use std::io::Write;
use std::sync::Mutex;

pub struct RSession {
    child: std::process::Child,
}

impl RSession {
    pub fn new() -> std::io::Result<Self> {
        let child = Command::new("R")
            .args(["--vanilla", "--quiet", "--slave"])
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()?;

        Ok(Self { child })
    }

    pub fn exec(&mut self, code: &str) -> std::io::Result<()> {
        let stdin = self.child.stdin.as_mut().unwrap();
        stdin.write_all(code.as_bytes())?;
        stdin.write_all(b"\n")?;
        stdin.flush()?;
        Ok(())
    }
}

//Functions
fn write_csv_header(path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file_path1 = format!("{}summaries.csv", path);

    let file = match File::create(&file_path1) {
        Ok(f) => f,
        Err(e) => return Err(Box::new(e)),
    };

    let buf_writer = BufWriter::new(file);
    let mut wtr = Writer::from_writer(buf_writer);

    wtr.write_record(&[
        "i",
        "pop_len", 
        "b_f", // Benefit of additive parental effort to fecundity
        "b_s", // mort
        "b_p", // brood predation risk
        "c_q", // survival cost of fast POL in winter
        "c_v", // surival cost of making observations in winter
        "c_psi", // opportunity cost of making observations on provisioning
        "c_u", // survival cost of parental effort in winter
        "mu", // mutation rate of loci
        "mut_size",
        "sigma0", // Locus determining prior alpha
        "sigma_cue",
        "div_rate",
        "varsigma", // the maximum sample size an individual is able to gather from cues
        "h", // hawkishness slope
        "theta",
        "mean_fitness",
        "m",
        "c",
        "c2",
        "u_base",
        "rho",
        "nu",
        "gamma",
        "lambda",
        "mean_faster_effort",
        "mean_slower_effort",
        "mean_fast_h",
        "mean_slow_h",
        "mean_f_s_w",
        "mean_f_s_d",
        "f_w",
        "s_w",
        "effort_dif",
        "mean_fast_obs",
        "mean_slow_obs"
    ])?;

    wtr.flush()?;
    Ok(())
}

fn write_csv_header2(path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file_path1 = format!("{}", path);

    let file = match File::create(&file_path1) {
        Ok(f) => f,
        Err(e) => return Err(Box::new(e)),
    };

    let buf_writer = BufWriter::new(file);
    let mut wtr = Writer::from_writer(buf_writer);

    wtr.write_record(&[
        "i",
        "pop_len", 
        "b_f", // Benefit of additive parental effort to fecundity
        "b_s", // mort
        "b_p", // brood predation risk
        "c_q", // survival cost of fast POL in winter
        "c_v", // surival cost of making observations in winter
        "c_psi", // opportunity cost of making observations on provisioning
        "c_u", // survival cost of parental effort in winter
        "mu", // mutation rate of loci
        "mut_size",
        "sigma0", // Locus determining prior alpha
        "sigma_cue",
        "div_rate",
        "varsigma", // the maximum sample size an individual is able to gather from cues
        "h", // hawkishness slope
        "theta",
        "mean_fitness",
        "m",
        "c",
        "c2",
        "u_base",
        "rho",
        "nu",
        "gamma",
        "lambda",
        "mean_faster_effort",
        "mean_slower_effort",
        "mean_fast_h",
        "mean_slow_h",
        "mean_f_s_w",
        "mean_f_s_d",
        "f_w",
        "s_w",
        "effort_dif",
        "mean_fast_obs",
        "mean_slow_obs"
    ])?;

    wtr.flush()?;
    Ok(())
}

fn write_csv_header3(path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file_path = format!("{path}summaries.csv");

    let file = OpenOptions::new()
        .write(true)
        .create_new(true) // prevents overwrite
        .open(&file_path)?;

    let buf_writer = BufWriter::new(file);
    let mut wtr = Writer::from_writer(buf_writer);

    wtr.write_record(["qm", "qf", "u1", "u2", "ben", "pred"])?;
    wtr.flush()?;

    Ok(())
}

fn run_r_sim_script(path: String) -> io::Result<()> {
    // Create a new command to run "Rscript" (must be installed and in PATH)
    let status = Command::new("Rscript")
        // Pass the path to your R script as the first argument
        .args(&["--vanilla", "--quiet", "./src/sim_plot.r", &path])
        .status()?;

    // Check if the script ran successfully (i.e., exited with code 0)
    if status.success() {
    } else {
        // If not, print the exit status
        eprintln!("❌ R script failed with status: {}", status);
    }

    Ok(()) // Return Ok if all went well
}

fn init_pop(n: usize, agent: Agent, sigma0:f64) -> Vec<Agent> {
    // returns vector with new agents
    let normal = Normal::new(0.0, sigma0).unwrap();
    let mut pop = Vec::new();
    for i in 0..n {
        pop.push(agent.clone());
        pop[i].q = normal.sample(&mut thread_rng());
        
    }
    return pop;
}

fn try_print(vec:Vec<f64>, file:&str) -> Result<(), Box<dyn Error>> {
    let strings: Vec<String> = vec.iter().map(|n| n.to_string()).collect();

    let outfile = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(file)
        .unwrap();
    let mut wtr = csv::Writer::from_writer(outfile);
    wtr.write_record(strings)?;
    wtr.flush()?;
    Ok(())
}

fn dice() -> f64 {
    let mut rng = rand::thread_rng();
    let dice: f64 = rng.gen(); 
    return dice
}

fn sigmoid(q:f64) -> f64 {
    return 1. / (1. + (-q).exp())
}

// Structs
#[derive(Clone, Copy, Debug)]
struct Agent {
    sigma: f64,
    pi: f64,
    q: f64,
    fitness: f64,
    u:f64, // negotiated parental effort
    partner:Option<usize>, // partner id
    m: f64, // Loci determining gradient of observation effort against pace-of-life
    c: f64, // Locus determining baseline obs rate
    c2: f64, // Locus determining max obs rate
    u_base: f64, // Locus - baseline parental effort
    rho: f64, // locus - influence of self POL information on negotiation
    nu: f64, // locus - influence of self POL information on negotiation
    lambda:f64, // locus - influence of prior partner bids on negotiation
    gamma:f64, // locus - influence of partner POL information on negotiation
}
#[derive(Clone, Debug)]
struct Environment {
    pop: Vec<Agent>,
    singles: Vec<usize>,
    dead: Vec<usize>,
    mean_fitness:f32,
    tol:f64,
    b_f:f64, // Benefit of additive parental effort to fecundity
    b_p:f64, // Offspring survival benefit of parental effort
    b_s:f64, // baseline survival rate of parents in winter
    c_q:f64, // survival cost of fast POL in winter
    c_v:f64, // surival cost of making observations in winter
    c_psi:f64, // surival cost of making observations in winter
    c_u:f64, // survival cost of parental effort in winter
    mu:f64, // mutation rate of loci
    sigma0: f64, // Prior sigma
    sigma_cue:f64, // variation of social cue
    varsigma:i64, // the maximum sample size an individual is able to gather from cues
    theta:f64, // slope of mortality against q-value
    mut_size:f64, // standard deviation of mutations
    div_rate:f64, // divorce rate
    h:f64, // hawkishness slope
    mean_faster_effort:f64,
    mean_slower_effort:f64,
    mean_fast_h:f64,
    mean_slow_h:f64,
    mean_f_s_w:f64,
    mean_f_s_d:f64,
    f_w:f64,
    s_w:f64,
    effort_dif:f64,
}

// Implementations
impl Agent {
    fn mutate(&mut self, mu:f64, mut_size:f64) {

        let mut rng = rand::thread_rng();

        let mut dice: f64 = rng.gen(); //m
        if dice <= mu { //mutates
            // mutate loci
            // loc0: f64, between 0.5 and 1. (proportion)
            let n = nm::new(self.m, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < -1. {
                self.m = -1.;
            } else if loc_mut > 1.{
                self.m = 1.;
            } else {
                self.m = loc_mut;
            }
        }
            
        dice = rng.gen(); //c
        if dice <= mu/5. { //mutates
            // loc2: f64 (c)
            // observation effort intercept
            let n = nm::new(self.c, mut_size).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < 0. {
                self.c = 0.;
            } else if loc_mut > 1.{
                self.c = 1.;
            } else {
                self.c = loc_mut;
            }
        }
        dice = rng.gen(); //c2
        if dice <= mu/5. { //mutates
            // loc2: f64 (c2)
            // observation effort intercept
            let n = nm::new(self.c2, mut_size).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < 0. {
                self.c2 = 0.;
            } else if loc_mut > 1.{
                self.c2 = 1.;
            } else {
                self.c2 = loc_mut;
            }
        }
        
        dice = rng.gen(); //ubase
        if dice <= mu { //mutates
            // loc2: f64 (u)
            // parental effort baseline
            let n = nm::new(self.u_base, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < -1. {
                self.u_base = -1.;
            } else if loc_mut > 1.{
                self.u_base = 1.;
            } else {
                self.u_base = loc_mut;
            }
        }

        dice = rng.gen(); //rho
        if dice <= mu { //mutates
            // loc2: f64 (rho)
            // parental effort selfishness bias
            let n = nm::new(self.rho, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < -1. {
                self.rho = -1.;
            } else if loc_mut > 1.{
                self.rho = 1.;
            } else {
                self.rho = loc_mut;
            }
        }

        dice = rng.gen(); //nu
        if dice <= mu { //mutates
            // loc2: f64 (nu)
            // pace-of-life influence effort baseline
            let n = nm::new(self.nu, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < -1. {
                self.nu = -1.;
            } else if loc_mut > 1.{
                self.nu = 1.;
            } else {
                self.nu = loc_mut;
            }
        }

        
        dice = rng.gen(); //lambda
        if dice <= mu { //mutates
            // loc3: f64, between 0. and 1. (lambda)
            // scales the perceived informativeness of partner gathering outcomes.
            let n = nm::new(self.lambda, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < -1.0 {
                self.lambda = -1.0;
            } else if loc_mut > 1.0 {
                self.lambda = 1.0;
            } else {
                self.lambda = loc_mut;
            }
        }
        
        dice = rng.gen(); //gamma
        if dice <= mu { //mutates
            // loc3: f64, between 0. and 1. (lambda)
            // scales the perceived informativeness of partner gathering outcomes.
            let n = nm::new(self.gamma, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < -1.0 {
                self.gamma = -1.0;
            } else if loc_mut > 1.0 {
                self.gamma = 1.0;
            } else {
                self.gamma = loc_mut;
            }
        }
    }

    fn pol_v(&self) -> f64 {
        // return self.c2*(1. - self.c) + self.c / (1. + (self.m * (self.q) ).exp());
        return self.c + self.c2*(1. - self.c)*sigmoid(self.m * (self.q));
    }

    fn social_cue(&mut self, sigma_cue:f64, phen:f64) {
        let b = 1. / (sigma_cue * sigma_cue);

        let a = 1. / (self.sigma * self.sigma);

        self.updating(a, b, phen);
    }

    fn updating(&mut self, a:f64, b:f64, x:f64) {
        self.pi = (a * self.pi + b * x) / (a + b);
        self.sigma = (1. / (a + b)).sqrt();
    }

    fn uncertainty(&self, sigma0:f64) -> f64 {
        // let hx = 0.5  * (2.0 * std::f64::consts::PI * std::f64::consts::E * self.sigma * self.sigma).log2();
        // return hx / (0.5 * (2.0 * std::f64::consts::PI * std::f64::consts::E * sigma0*sigma0).log2())
        return self.sigma / sigma0
    }

    fn r(&self, u2:f64, h:f64) -> f64 {
        let r:f64
            = self.rho
            + self.nu*self.q
            + self.gamma*(1.-h)*(self.pi-self.q) 
            - self.lambda*(u2 - self.u_base);
        return r
    }

    fn surivorship(&self, b_s:f64, c_q:f64, c_v:f64, c_u:f64, theta:f64) -> f64 {
        let s: f64 =  b_s*(1. - sigmoid(self.q * theta) * c_q)*(1. - self.pol_v() * c_v)*(1. - self.u * c_u);
        return s
    }
}

impl Environment {
    fn observations(&mut self) {
        let mut v_f: f64;
        let mut q_m: f64;
        for i in 0..self.pop.len(){
            let female = i;
            let male = self.pop[i].partner.unwrap();
            v_f = self.pop[female].pol_v();

            q_m = self.pop[male].q;
            let normal = Normal::new(q_m, self.sigma_cue).unwrap();
            
            for _j in 0..self.varsigma {
                if dice() < v_f {
                    let phen:f64 = normal.sample(&mut thread_rng());
                    self.pop[female].social_cue(self.sigma_cue, phen);
                }
            }
        } 
    }

    fn negotiations(&mut self) {  
        // inits
            let mut u1:f64;
            let mut u2:f64;
            let mut h1:f64;
            let mut h2:f64;
            // let mut old:f64;
            let mut male:usize; 
            let mut female:usize;
            let mut ben:f64;
            let mut q1:f64;
            let mut q2:f64;
            let mut pred:f64;
            self.mean_faster_effort = 0.;
            self.mean_slower_effort = 0.;
            self.mean_fast_h = 0.;
            self.mean_slow_h = 0.;
            let mut fast_count: f64 = 0.;
            let mut slow_count: f64 = 0.;
            let mut fitness: f64;
            let mut dif: f64 = 0.;
            self.mean_f_s_w = 0.;
            self.mean_f_s_d = 0.;
            let mut f_count: f64 = 0.;
            let mut s_count: f64 = 0.;
            self.effort_dif = 0.;
            let mut count: f64 = 0.;

        for i in 0..self.pop.len(){
            female = i;
            male = self.pop[i].partner.unwrap();
            if male > female { // so that there aren't repeats
                h1 = self.pop[male].uncertainty(self.sigma0);
                h2 = self.pop[female].uncertainty(self.sigma0);
                u1=0.;
                u2=0.;
                
                // Parental egotiations
                for _j in 0..100 {
                    u1 = self.pop[male].r(u2,h1);
                    u2 = self.pop[female].r(u1,h2);
                }
                let alpha1 = sigmoid(2.*u1)*(1. - self.c_psi*self.pop[male].pol_v());
                let alpha2 = sigmoid(2.*u2)*(1. - self.c_psi*self.pop[female].pol_v());
                let beta1 = (1.-sigmoid(2.*u1))*(1. - self.c_psi*self.pop[male].pol_v());
                let beta2 = (1.-sigmoid(2.*u2))*(1. - self.c_psi*self.pop[female].pol_v());
                if (alpha1 < 0.) | (alpha2 < 0.) {
                    print!("\n{}, {}", u1, u2)
                }
                // Provisioning
                ben = ((alpha1.max(0.0) + alpha2.max(0.0))/2.).powf(0.5);

                // Predation
                self.pop[male].u = alpha1;
                self.pop[female].u = alpha2;
                q1 = self.pop[male].q;
                q2 = self.pop[female].q;
                pred = self.predation(beta1, beta2, q1, q2);

                // Record data
                {
                    if q1 > self.sigma0 {
                        self.mean_faster_effort += alpha1;
                        self.mean_fast_h += h1;
                        fast_count += 1.;
                    } else if q1 < -self.sigma0 {
                        self.mean_slower_effort += alpha1;
                        self.mean_slow_h += h1;
                        slow_count += 1.;
                    }
                    if q2 > self.sigma0 {
                        self.mean_faster_effort += alpha2;
                        self.mean_fast_h += h2;
                        fast_count += 1.;
                    } else if q2 < -self.sigma0 {
                        self.mean_slower_effort += alpha2;
                        self.mean_slow_h += h2;
                        slow_count += 1.;
                    }
                }

                // Fledging
                fitness = (ben * (1.-pred) / 2.).max(0.0);
                self.pop[male].fitness = fitness;
                self.pop[female].fitness = fitness;

                // Record more data
                {
                    if (q1 > 0.) & (q2 < 0.) {
                        dif += 1.;
                        self.mean_f_s_w += fitness;
                        self.mean_f_s_d += alpha1 / (alpha1+alpha2);
                    }
                    else if (q1 < 0.) & (q2 > 0.) {
                        dif += 1.;
                        self.mean_f_s_w += fitness;
                        self.mean_f_s_d += alpha2 / (alpha1+alpha2);
                    } else if (q1 > 0.) & (q2 > 0.) {
                        f_count += 1.;
                        self.f_w += fitness;
                    } else if (q1 < 0.) & (q2 <0.) {
                        s_count += 1.;
                        self.s_w += fitness;
                    }
                    if alpha1 > alpha2 {
                        self.effort_dif += alpha1 - alpha2;
                        count+=1.;
                    } else {
                        self.effort_dif += alpha2 - alpha1;
                        count+=1.;
                    }
                }
            }
        }

        self.f_w /= f_count;
        self.s_w /= s_count;
        self.mean_f_s_w /= dif;
        self.mean_f_s_d /= dif;
        self.mean_faster_effort /= fast_count;
        self.mean_slower_effort /= slow_count;
        self.mean_fast_h /= fast_count;
        self.mean_slow_h /= slow_count;
        self.effort_dif /= count;
    }
    
    fn negotiations_baseline(&mut self, file:&str) {  
        // determine male and female parental effort and then update fitness
        // let mut dif:f64;
        let mut u1:f64;
        let mut u2:f64;
        let mut h1:f64;
        let mut h2:f64;
        // let mut old:f64;
        let mut male:usize; 
        let mut female:usize;
        let mut ben:f64;
        let mut q1:f64;
        let mut q2:f64;
        let mut pred:f64;

        for i in 0..self.pop.len(){
            let mut data = Vec::new();
            female = i;
            male = self.pop[i].partner.unwrap();
            if male > female { // so that there aren't repeats
                // let mut gen: f64 = 0.;
                // dif = 9999.;
                h1 = self.pop[male].uncertainty(self.sigma0);
                h2 = self.pop[female].uncertainty(self.sigma0);
                u1=0.;
                u2=0.;
                // old=u1+u2;
                // negotiations:
                // while dif > self.tol {
                for _j in 0..100 {
                    
                    u1 = self.pop[male].r(u2,h1);
                    u2 = self.pop[female].r(u1,h2);
                    
                }
                u1 = sigmoid(2.*u1);
                u2 = sigmoid(2.*u2);
                let alpha1 = u1*(1. - self.c_psi*self.pop[male].pol_v());
                let alpha2 = u2*(1. - self.c_psi*self.pop[female].pol_v());
                let beta1 = (1.-u1)*(1. - self.c_psi*self.pop[male].pol_v());
                let beta2 = (1.-u2)*(1. - self.c_psi*self.pop[female].pol_v());
                if (alpha1 < 0.) | (alpha2 < 0.) {
                    print!("\n{}, {}", u1, u2)
                }
                // Provisioning
                ben = ((alpha1.max(0.0) + alpha2.max(0.0))/2.).powf(0.5);

                // Predation
                self.pop[male].u = alpha1;
                self.pop[female].u = alpha2;
                q1 = self.pop[male].q;
                q2 = self.pop[female].q;
                pred = self.predation(beta1, beta2, q1, q2);

                data.push(q1);
                data.push(q2);
                data.push(u1);
                data.push(u2);
                data.push(ben);
                data.push(pred);
                let _ = try_print(data, file);
                // Record data
            }
        }
    }

    fn predation(&mut self, beta1:f64,beta2:f64,q1:f64,q2:f64) -> f64 {
        let g1 = sigmoid(self.h*q1);
        let g2 = sigmoid(self.h*q2);
        return self.b_p * (-(beta1*g1 + beta2*g2)).exp()
    }

    fn kills(&mut self, i: usize) {
        // individual i is killed by the environment
        if self.dead.contains(&i){
            return ();
        }
        self.dead.push(i);
    }

    fn widowed(&mut self) {
        for i in 0..self.dead.len() {
            let partner_id: Option<usize> = self.pop[self.dead[i]].partner;
            match partner_id {
                None => (),
                _ => {
                    self.pop[self.dead[i]].partner = None;
                    self.pop[partner_id.unwrap()].partner = None;
                }
            }
        }
    }

    fn winter_deaths(&mut self) {
        // winter deaths given phenotype environment match of each individual
        let mut rng = rand::thread_rng();
        for i in 0..self.pop.len() {
            if rng.gen::<f64>() > self.pop[i].surivorship(self.b_s, self.c_q, self.c_v, self.c_u, self.theta) {
                self.kills(i);
            }
        }
    }

    fn divorce(&mut self) { // divorce with prob div_rate
        
        let mut rng = rand::thread_rng();

        for i in 0..self.pop.len() {
            
            let dice: f64 = rng.gen(); 
            if dice < self.div_rate {
                match self.pop[i].partner {
                    None => (),
                    _ => {
                        let partner_id: usize = self.pop[i].partner.unwrap();
                        self.pop[partner_id].partner = None;
                        self.pop[i].partner = None;
                    }
                }
            }
        }
    }

    fn pairing_pool(&mut self) { // divorces and then all attempt to pair
        // all pairs divorce
        self.divorce();

        // Collect only those who are actually single
        let mut available: Vec<usize> = self.singles
            .iter()
            .copied()
            .filter(|i| self.pop[*i].partner.is_none())
            .collect();

        available.shuffle(&mut thread_rng()); // shuffle the true singles

        // Pair them in adjacent pairs
        for chunk in available.chunks(2) {
            if let [a, b] = chunk {
                // initialize states
                self.pop[*a].sigma = self.sigma0;
                self.pop[*a].pi = 0.0;
                self.pop[*a].partner = Some(*b);

                self.pop[*b].sigma = self.sigma0;
                self.pop[*b].pi = 0.0;
                self.pop[*b].partner = Some(*a);
            }
        }
    }

    fn natural_selection(&mut self) {
        let mut selected_male:usize;
        let mut selected_female:usize;
        let normal = Normal::new(0.0, self.sigma0).unwrap();
        // for each i in dead, replace strategy with that of a living individual (single or paired) biased by their fitness from either environment

        // weighted agent indices by fitness (dead or alive)
        let mut weights: Vec<f32> = self.pop.iter().map(|a| a.fitness as f32).collect();
        let tot: f32 = weights.iter().sum();
        self.mean_fitness = tot / weights.len() as f32;
        // println!("mean: {}, tot: {}, weights.len(): {}", self.mean_fitness, tot, weights.len());
        for w in 0..weights.len() {
            // println!("{}", weights[w]);
            weights[w] = weights[w]/tot;
            // println!("{}", weights[w]);
        }

        // table for random weighted choice of agent indices by fitness weights
        let builder = WalkerTableBuilder::new(&weights);
        let wa_table = builder.build();

        self.winter_deaths(); // over winter deaths as a result of degree of fit to the environment 
                             // and baseline mortality rate
        // seed
        let mut rng = rand::thread_rng();
        // for each dead agent, give new policy.
        for i in 0..self.dead.len() {
            // choose parent
            selected_male = wa_table.next_rng(&mut rng);
            // println!("{}, fit: {}",selected_male, self.pop[selected_male].fitness);
            selected_female = self.pop[selected_male].partner.unwrap();

            let parent = if rng.gen::<f64>() < 0.5 { selected_male } else { selected_female };
            self.pop[self.dead[i]].m= self.pop[parent].m;
            self.pop[self.dead[i]].c = self.pop[parent].c;
            self.pop[self.dead[i]].c2 = self.pop[parent].c2;
            self.pop[self.dead[i]].u_base = self.pop[parent].u_base;
            self.pop[self.dead[i]].rho = self.pop[parent].rho;
            self.pop[self.dead[i]].nu = self.pop[parent].nu;
            self.pop[self.dead[i]].lambda = self.pop[parent].lambda;
            self.pop[self.dead[i]].gamma = self.pop[parent].gamma;
            self.pop[self.dead[i]].fitness = 0.;
            self.pop[self.dead[i]].q = normal.sample(&mut thread_rng());

            // apply mutations to policies
            self.pop[self.dead[i]].mutate(self.mu, self.mut_size);   
        }
        self.widowed(); // unpair widows and the dead
        self.dead.clear(); // all dead should be replaced so clear
    }

    fn csvgen(&mut self, pop_len:f64, i:f64) -> Vec<f64> {
        let mut data: Vec<f64> = Vec::new();
        let mut mean_loc_m: f64 = 0.;
        let mut mean_loc_c: f64 = 0.;
        let mut mean_loc_c2: f64 = 0.;
        let mut mean_loc_u_base: f64 = 0.;
        let mut mean_loc_rho: f64 = 0.;
        let mut mean_loc_nu: f64 = 0.;
        let mut mean_loc_gamma: f64 = 0.;
        let mut mean_loc_lambda: f64 = 0.;
        let mut mean_fast_obs:f64 = 0.;
        let mut mean_slow_obs:f64 = 0.;
        let mut fast_count: f64 = 0.;
        let mut slow_count: f64 = 0.;
        // println!("Generation: {}, Mean_fitness: {:.4}, prior_wet: {:.4}, prior_general: {:.4}, prior_dry: {:.4}, prop_obs: {:.4}", 
            // i, env.mean_fitness, env.prop_obs);
        for a in &self.pop {
            mean_loc_m += a.m;
            mean_loc_c += a.c;
            mean_loc_c2 += a.c2;
            mean_loc_u_base += a.u_base;
            mean_loc_rho += a.rho;
            mean_loc_nu += a.nu;
            mean_loc_lambda += a.lambda;
            mean_loc_gamma += a.gamma;
            if a.q > 0. {
                mean_fast_obs += a.pol_v();
                fast_count += 1.;
            } else {
                mean_slow_obs += a.pol_v();
                slow_count += 1.;
            }
        }
        mean_loc_m = mean_loc_m/pop_len;
        mean_loc_c = mean_loc_c/pop_len;
        mean_loc_c2 = mean_loc_c2/pop_len;
        mean_loc_u_base = mean_loc_u_base/pop_len;
        mean_loc_rho = mean_loc_rho/pop_len;
        mean_loc_nu = mean_loc_nu/pop_len;
        mean_loc_lambda = mean_loc_lambda/pop_len;
        mean_loc_gamma = mean_loc_gamma/pop_len;
        mean_fast_obs = mean_fast_obs/fast_count;
        mean_slow_obs = mean_slow_obs/slow_count;
            
        data.push(i);
        data.push(pop_len);
        data.push(self.b_f); // Baseline fecundity
        data.push(self.b_s); // parental mort
        data.push(self.b_p); // predation risk
        data.push(self.c_q); // survival cost of fast POL in winter
        data.push(self.c_v); // surival cost of making observations in winter
        data.push(self.c_psi); // surival cost of making observations in winter
        data.push(self.c_u); // survival cost of parental effort in winter
        data.push(self.mu); // mutation rate of loci
        data.push(self.mut_size);
        data.push(self.sigma0); // Locus determining prior alpha
        data.push(self.sigma_cue);
        data.push(self.div_rate);
        data.push(self.varsigma as f64); // the maximum sample size an individual is able to gather from cues
        data.push(self.h);
        data.push(self.theta);
        data.push(self.mean_fitness as f64);
        data.push(mean_loc_m);
        data.push(mean_loc_c);
        data.push(mean_loc_c2);
        data.push(mean_loc_u_base);
        data.push(mean_loc_rho);
        data.push(mean_loc_nu);
        data.push(mean_loc_gamma);
        data.push(mean_loc_lambda);
        data.push(self.mean_faster_effort);
        data.push(self.mean_slower_effort);
        data.push(self.mean_fast_h);
        data.push(self.mean_slow_h);
        data.push(self.mean_f_s_w / self.mean_fitness as f64);
        data.push(self.mean_f_s_d);
        data.push(self.f_w / self.mean_fitness as f64);
        data.push(self.s_w / self.mean_fitness as f64);
        data.push(self.effort_dif);
        data.push(mean_fast_obs);
        data.push(mean_slow_obs);
        return data
    }

}

// main Functions
fn run(
    generations:i32, 
    res_path: &String, 
    trial: Option<&str>, 
    gen: Option<&i32>,
    mut env:Environment) {

    let trial = trial.unwrap_or("");
    let file_path2 = format!("{}_{}_{}policies.csv", res_path,trial,gen.unwrap_or(&0));
    let _ = write_csv_header2(&file_path2);
    
    // let mut r = RSession::new()?;
    // r.exec("source('src/sim_plot.r')")?;

    /////////////////////////////////////////// Run simulation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    for i in 1..=generations {
        // env.evolve_priors(); // updates priors for beliefs about partner phenotype (prior for climate is const)
        env.pairing_pool(); // available individuals pair up randomly
        env.observations(); // Individuals decide observation effort
        env.negotiations(); // individuals negotiate their parental effort
        
        if i%1000 == 0 {
            let data: Vec<f64> = env.csvgen(env.pop.len() as f64, i as f64);
            let _ = try_print(data, &file_path2);
            // let _ = run_r_sim_script(file_path2.clone());
        }

        if i==generations-1 {
            let data: Vec<f64> = env.csvgen(env.pop.len() as f64, i as f64);
            let file_path1 = format!("{}summaries.csv", res_path);
            let _ = try_print(data, &file_path1);
        }
       
        env.natural_selection(); // selection, reproduction and mutation of the fittest
    }
}

fn run2(
    generations:i32, 
    res_path: &String,
    mut env:Environment) {

    /////////////////////////////////////////// Run simulation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    for i in 1..=generations {
        // env.evolve_priors(); // updates priors for beliefs about partner phenotype (prior for climate is const)
        env.pairing_pool(); // available individuals pair up randomly
        env.observations(); // Individuals decide observation effort
        env.negotiations(); // individuals negotiate their parental effort

        if i==generations-1 {
            let file_path1 = format!("{}summaries.csv", res_path);
            env.negotiations_baseline(&file_path1); // print care outcomes of negotiations
        }
       
        env.natural_selection(); // selection, reproduction and mutation of the fittest
    }
}

fn main() -> std::io::Result<()>  {
////////////////////////////////////////////////////// Comand-line arguments
        let args: Vec<String> = env::args().collect();
        let project_id = &args[1];

/////////////////////////////////////////// Initialise baseline parameters \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        let generations: i32 = 30000; // number of generations the sim lasts
        let sigma0:f64 = 2.0;
        let iterations:i32 = 500;
        let pop_size = 1000;

        let agent:Agent = Agent {
            sigma: sigma0,
            pi: 0.,
            q: 0.5,
            fitness: 0.,
            u:0., 
            partner:None,
            m: 0., // Loci determining gradient of observation effort against pace-of-life
            c: 0., // Locus determining effort min of observation effort against POL
            c2: 0., // Locus determining effort max of observation effort against POL
            u_base: 1., // locus - perceived optimal partner effort
            rho: 1., // locus - starting point of parental negotiation
            nu: 1., // locus - influence of self POL information on negotiation
            lambda: 0., // locus - influence of prior partner bids on negotiation
            gamma: -0.5, // locus - influence of partner POL information on negotiation
        };
        
        let env = Environment {
            pop: init_pop(pop_size, agent,sigma0),
            singles: (0..pop_size as usize).collect(),
            dead: Vec::new(),
            mean_fitness:0.,
            tol: 0.001,
            mu:0.1, // mutation rate of loci
            mut_size:0.2,

            b_f:1.0, // Baseline fecundity (offspring independence)
            b_s:0.95, // Baseline survival rate of adults (5 years lifespan)
            b_p:0., // Brood predation risk (assuming twice the mortality of adults)
            c_q:0.4, // survival cost of fast POL in winter 
            c_v:0.0, // surival cost of making observations in winter
            c_psi:1.0, // opportunity cost of making observations on provisioning
            c_u:0.25, // survival cost of parental effort in winter (highest energy expenditure)
            sigma0, // Locus determining prior alpha
            sigma_cue:0.5, // sd of the cue
            varsigma:100, // the maximum sample size an individual is able to gather from cues
            h:1., // slow-fast slope of nest-defence/size/aggression (behavioural phenotypes more extreme therefore h>theta)
            theta: 1.0, // slow-fast sigmoid slope of mortality risk against q-value (physiological
            div_rate:0.1, // divorce rate
            mean_faster_effort:0.,
            mean_slower_effort:0.,
            mean_fast_h:0.,
            mean_slow_h:0.,
            mean_f_s_d:0.,
            mean_f_s_w:0.,
            f_w:0.,
            s_w:0.,
            effort_dif:0.,
        };

////////////////////////////////////////////////////// start r session
    let mut r = RSession::new()?;
    r.exec("source('src/plots.r')")?;
    let r_mutex = Mutex::new(r);

////////////////////////////// Cue detectability (cue_sigma) sims
    // Construct the full path
    let path = format!("./Results/{}/sigma_cue/", project_id);
    let path_construct = Path::new(&path);
     // Ensure parent directory exists
    let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
    let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.001..50.0); // uniform sample from parameter space
        env.sigma_cue = x;
        
        println!("Simulation started: sigmacue: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: sigmacue: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 5 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_sigmacue_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_sigmacue_plot('{}')", path)).unwrap();

////////////////////////////// Brood predation risk (b_p) sims
    // Construct the full path
        let path = format!("./Results/{}/b_p/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.b_p = x;
        
        println!("Simulation started: predation risk: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: predation risk: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_b_p_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_b_p_plot('{}')", path)).unwrap();


////////////////////////////// Adult mortality rate (b_s) simulations
    // Construct the full path
    let path = format!("./Results/{}/b_s/", project_id);
    let path_construct = Path::new(&path);
    // Ensure parent directory exists
    let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
    let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.0..100.0); // uniform sample from parameter space
        env.b_s = 1.-1./x;
        
        println!("Simulation started: mort: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: mort: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_b_s_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_b_s_plot('{}')", path)).unwrap();


////////////////////////////// Observation opportunity cost c_psi
    // Construct the full path
    let path = format!("./Results/{}/c_psi/", project_id);
    let path_construct = Path::new(&path);
    // Ensure parent directory exists
    let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
    let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.c_v = 0.;
        env.c_psi = x;
        
        println!("Simulation started: c_psi: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: c_psi: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_c_psi_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_c_psi_plot('{}')", path)).unwrap();
    
////////////////////////////// Slope of mort for POL theta
        // Construct the full path
        let path = format!("./Results/{}/theta/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.0..2.0); // uniform sample from parameter space
        env.theta = x;
        
        println!("Simulation started: theta: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: theta: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_theta_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_theta_plot('{}')", path)).unwrap();

////////////////////////////// Slope of defence for POL h
        // Construct the full path
        let path = format!("./Results/{}/h/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.0..2.0); // uniform sample from parameter space
        env.h = x;
        
        println!("Simulation started: h: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: h: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_h_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_h_plot('{}')", path)).unwrap();


////////////////////////////// Baseline high
    // Construct the full path
        let path = format!("./Results/{}/Baseline_high/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header3(&path);
    (0..50).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        env.b_p=1.;
        
        println!("Simulation started: gen: {}", g);
        run2(
            generations, 
            &path,  
            env
        );
        println!("Simulation done: gen: {}", g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_BaselineH_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_BaselineH_plot('{}')", path)).unwrap();

////////////////////////////// Baseline low
    // Construct the full path
        let path = format!("./Results/{}/Baseline_low/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header3(&path);
    (0..50).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        env.b_p=0.;
        
        println!("Simulation started: gen: {}", g);
        run2(
            generations, 
            &path,  
            env
        );
        println!("Simulation done: gen: {}", g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_BaselineL_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_BaselineL_plot('{}')", path)).unwrap();

////////////////////////////// Cost of Fast POLS c_q
    // Construct the full path
        let path = format!("./Results/{}/c_q/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.c_q = x;
        
        println!("Simulation started: c_q: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: c_q: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_c_q_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_c_q_plot('{}')", path)).unwrap();


////////////////////////////// Cue cost sims
        // Construct the full path
        let path = format!("./Results/{}/cue_cost/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.c_v = x;
        
        println!("Simulation started: cue cost: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: cue cost: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_cuecost_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_cuecost_plot('{}')", path)).unwrap();    

////////////////////////////// Pace-of-life variance (sigma) sims
        // Construct the full path
        let path = format!("./Results/{}/sigma/", project_id);
        let path_construct = Path::new(&path);
        // Ensure parent directory exists
        let _ = fs::create_dir_all(path_construct); // Create directory path if it doesn't exist
        let _ = write_csv_header(&path);
    (0..iterations).into_par_iter().for_each(|g|  {
        // Initialise stochastic variables
        let mut rng = rand::thread_rng();
        let mut env = env.clone();
        let mut agent = agent.clone();
        agent.mutate(1.0, 1.0); // randomize resident loci
        env.pop = init_pop(pop_size, agent, sigma0);
        let x = rng.gen_range(0.01..20.0); // uniform sample from parameter space
        env.sigma0 = x;
        
        println!("Simulation started: sigma: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: sigma: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_sigma_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_sigma_plot('{}')", path)).unwrap();
    

    Ok(())
}