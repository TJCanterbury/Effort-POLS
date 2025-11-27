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
        "c_u", // survival cost of parental effort in winter
        "mu", // mutation rate of loci
        "mut_size",
        "sigma0", // Locus determining prior alpha
        "mu_q",
        "sigma_cue",
        "div_rate",
        "varsigma", // the maximum sample size an individual is able to gather from cues
        "h", // hawkishness slope
        "theta",
        "mean_fitness",
        "m",
        "c",
        "u_base",
        "rho",
        "nu",
        "gamma",
        "lambda"
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
        "c_u", // survival cost of parental effort in winter
        "mu", // mutation rate of loci
        "mut_size",
        "sigma0", // Locus determining prior alpha
        "sigma_cue",
        "mu_q",
        "div_rate",
        "varsigma", // the maximum sample size an individual is able to gather from cues
        "h", // hawkishness slope
        "theta",
        "mean_fitness",
        "m",
        "c",
        "u_base",
        "rho",
        "nu",
        "gamma",
        "lambda"
    ])?;

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
    return 1. / (1. + (q).exp())
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
    c: f64, // Locus determining effort intercept of observation effort against POL
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
    c_u:f64, // survival cost of parental effort in winter
    mu:f64, // mutation rate of loci
    sigma0: f64, // Prior sigma
    mu_q:f64, // prior pi
    sigma_cue:f64, // variation of social cue
    varsigma:i64, // the maximum sample size an individual is able to gather from cues
    theta:f64, // slope of mortality against q-value
    mut_size:f64, // standard deviation of mutations
    div_rate:f64, // divorce rate
    h:f64, // hawkishness slope
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
        if dice <= mu { //mutates
            // loc2: f64 (c)
            // observation effort intercept
            let n = nm::new(self.c, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < 0. {
                self.c = 0.;
            } else if loc_mut > 1.{
                self.c = 1.;
            } else {
                self.c = loc_mut;
            }
        }
        
        dice = rng.gen(); //ubase
        if dice <= mu { //mutates
            // loc2: f64 (u)
            // parental effort baseline
            let n = nm::new(self.u_base, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            if loc_mut < 0. {
                self.u_base = 0.;
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
            if loc_mut < 0. {
                self.rho = 0.;
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
            self.nu = loc_mut;
        }

        
        dice = rng.gen(); //lambda
        if dice <= mu { //mutates
            // loc3: f64, between 0. and 1. (lambda)
            // scales the perceived informativeness of partner gathering outcomes.
            let n = nm::new(self.lambda, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            self.lambda = loc_mut;
        }
        
        dice = rng.gen(); //gamma
        if dice <= mu { //mutates
            // loc3: f64, between 0. and 1. (lambda)
            // scales the perceived informativeness of partner gathering outcomes.
            let n = nm::new(self.gamma, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            self.gamma = loc_mut;
        }
    }

    fn pol_v(&self) -> f64 {
        // return self.c2*(1. - self.c) + self.c / (1. + (self.m * (self.q) ).exp());
        return self.m * self.q + self.c;
    }

    fn social_cue(&mut self, partner_q:f64, sigma_cue:f64) {
        let normal = Normal::new(partner_q, sigma_cue).unwrap();
        let phen:f64;
        let b = 1. / (sigma_cue * sigma_cue);

        // println!("prior pi: {}, sig: {}, q: {}, sig_cue: {}", self.pi, self.sigma, partner_q, sigma_cue);
        let a = 1. / (self.sigma * self.sigma);
        phen = normal.sample(&mut thread_rng());
        self.updating(a, b, phen);
            // println!("n: {}, post pi: {}, sig: {}, phen: {}, q: {}, H(q): {}\n",i, self.pi, self.sigma, phen, partner_q, self.uncertainty());
    }

    fn updating(&mut self, a:f64, b:f64, x:f64) {
        self.pi = (a * self.pi + b * x) / (a + b);
        self.sigma = (1. / (a + b)).sqrt();
    }

    fn uncertainty(&self, sigma0:f64) -> f64 {
        let hx = 0.5  * (2.0 * std::f64::consts::PI * std::f64::consts::E * self.sigma * self.sigma).log2();
        return hx / (0.5 * (2.0 * std::f64::consts::PI * std::f64::consts::E * sigma0*sigma0).log2())
    }

    fn r(&self, u2:f64, h:f64) -> f64 {
        let r:f64
            = self.rho
            + self.nu*self.q
            + self.gamma*(self.q - self.pi) 
            - self.lambda*(u2 - self.u_base);
        return r.clamp(0., 1.)
    }

    fn surivorship(&self, b_s:f64, c_q:f64, c_v:f64, c_u:f64, theta:f64) -> f64 {
        let s: f64 =  b_s*(1. - sigmoid(-theta*(self.q-0.5)) * c_q)*(1. - self.pol_v() *c_v)*(1. - self.u * c_u); 
        return s
    }
}

impl Environment {
    fn predation(&mut self, u1:f64,u2:f64,q1:f64,q2:f64) -> f64 {
        return self.b_p * (1. - sigmoid(self.h*(q1-0.5))*(1.-u1)) * (1. - sigmoid(self.h*(q2-0.5))*(1.-u2))
    }

    fn negotiations(&mut self) {  
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
                    // Observation preference
                    // gen+=1.;
                    u1 = self.pop[male].r(u2,h1);
                    u2 = self.pop[female].r(u1,h2);
                }

                ben = self.b_f + u1 + u2;
                self.pop[male].u = u1;
                self.pop[female].u = u2;
                q1 = self.pop[male].q;
                q2 = self.pop[female].q;
                pred = self.predation(u1, u2, q1, q2);
                
                self.pop[male].fitness = (ben * (1.-pred) / 2.).max(0.0);
                self.pop[female].fitness = (ben * (1.-pred) / 2.).max(0.0);
            }
        }
    }
    
    fn observations(&mut self) {
        let mut v_f: f64;
        let mut q_m: f64;
        let mut n: f64;
        let mut xbar:f64;
        for i in 0..self.pop.len(){
            n=0.;
            let female = i;
            let male = self.pop[i].partner.unwrap();
            v_f = self.pop[female].pol_v();
            q_m = self.pop[male].q;
            let normal = Normal::new(q_m, self.sigma_cue).unwrap();
            xbar = 0.;
            
            // Number of observations made
            for _j in 0..self.varsigma {
                if dice() < v_f {
                    n += 1.;
                    
                    xbar += normal.sample(&mut thread_rng());
                }
            }

            // observation mean xbar
            xbar = xbar / n;

            // Bayesian updating
            self.pop[female].sigma = (1./self.sigma0.powf(2.) + n/self.sigma_cue.powf(2.)).powf(-1.).powf(1./2.);
            self.pop[female].pi = self.pop[female].sigma * (self.mu_q/self.sigma0.powf(2.) + (n*xbar)/self.sigma_cue.powf(2.));
        } 
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
                self.pop[*a].pi = self.mu_q;
                self.pop[*a].partner = Some(*b);

                self.pop[*b].sigma = self.sigma0;
                self.pop[*b].pi = self.mu_q;
                self.pop[*b].partner = Some(*a);
            }
        }
    }

    fn natural_selection(&mut self) {
        let mut selected_male:usize;
        let mut selected_female:usize;
        let blend:f64 = 0.5;
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

            // Blending:
            self.pop[self.dead[i]].m = blend*self.pop[selected_male].m + (1.-blend)*self.pop[selected_female].m ;
            self.pop[self.dead[i]].c = blend*self.pop[selected_male].c + (1.-blend)*self.pop[selected_female].c;
            self.pop[self.dead[i]].u_base = blend*self.pop[selected_male].u_base + (1.-blend)*self.pop[selected_female].u_base;
            self.pop[self.dead[i]].rho = blend*self.pop[selected_male].rho + (1.-blend)*self.pop[selected_female].rho;
            self.pop[self.dead[i]].nu = blend*self.pop[selected_male].nu + (1.-blend)*self.pop[selected_female].nu ;
            self.pop[self.dead[i]].lambda = blend*self.pop[selected_male].lambda + (1.-blend)*self.pop[selected_female].lambda;
            self.pop[self.dead[i]].gamma = blend*self.pop[selected_male].gamma + (1.-blend)*self.pop[selected_female].gamma;
            self.pop[self.dead[i]].fitness = 0.;
            self.pop[self.dead[i]].q = normal.sample(&mut thread_rng());

            // apply mutations to policies
            self.pop[self.dead[i]].mutate(self.mu, self.mut_size);   
        }
        self.widowed(); // unpair widows and the dead
        self.dead.clear(); // all dead should be replaced so clear
    }

    fn csvgen(&mut self, pop_len:f64, i:f64) -> Vec<f64> {
        let mut data = Vec::new();
        let mut mean_loc_m = 0.;
        let mut mean_loc_c = 0.;
        let mut mean_loc_u_base = 0.;
        let mut mean_loc_rho = 0.;
        let mut mean_loc_nu = 0.;
        let mut mean_loc_gamma = 0.;
        let mut mean_loc_lambda = 0.;
        // println!("Generation: {}, Mean_fitness: {:.4}, prior_wet: {:.4}, prior_general: {:.4}, prior_dry: {:.4}, prop_obs: {:.4}", 
            // i, env.mean_fitness, env.prop_obs);
        for a in &self.pop {
            mean_loc_m += a.m;
            mean_loc_c += a.c;
            mean_loc_u_base += a.u_base;
            mean_loc_rho += a.rho;
            mean_loc_nu += a.nu;
            mean_loc_lambda += a.lambda;
            mean_loc_gamma += a.gamma;
        }
        mean_loc_m = mean_loc_m/pop_len;
        mean_loc_c = mean_loc_c/pop_len;
        mean_loc_u_base = mean_loc_u_base/pop_len;
        mean_loc_rho = mean_loc_rho/pop_len;
        mean_loc_nu = mean_loc_nu/pop_len;
        mean_loc_lambda = mean_loc_lambda/pop_len;
        mean_loc_gamma = mean_loc_gamma/pop_len;
            
        data.push(i);
        data.push(pop_len);
        data.push(self.b_f); // Baseline fecundity
        data.push(self.b_s); // parental mort
        data.push(self.b_p); // predation risk
        data.push(self.c_q); // survival cost of fast POL in winter
        data.push(self.c_v); // surival cost of making observations in winter
        data.push(self.c_u); // survival cost of parental effort in winter
        data.push(self.mu); // mutation rate of loci
        data.push(self.mut_size);
        data.push(self.sigma0); // Locus determining prior alpha
        data.push(self.mu_q);
        data.push(self.sigma_cue);
        data.push(self.div_rate);
        data.push(self.varsigma as f64); // the maximum sample size an individual is able to gather from cues
        data.push(self.h);
        data.push(self.theta);
        data.push(self.mean_fitness as f64);
        data.push(mean_loc_m);
        data.push(mean_loc_c);
        data.push(mean_loc_u_base);
        data.push(mean_loc_rho);
        data.push(mean_loc_nu);
        data.push(mean_loc_gamma);
        data.push(mean_loc_lambda);
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

fn main() -> std::io::Result<()>  {
////////////////////////////////////////////////////// Comand-line arguments
        let args: Vec<String> = env::args().collect();
        let project_id = &args[1];

/////////////////////////////////////////// Initialise baseline parameters \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        let generations: i32 = 20000; // number of generations the sim lasts
        let sigma0:f64 = 1./20.;
        let iterations:i32 = 200;
        let pop_size = 1000;

        let agent:Agent = Agent {
            sigma: sigma0,
            pi: 100.,
            q: 0.5,
            fitness: 0.,
            u:0., 
            partner:None,
            m: 0., // Loci determining gradient of observation effort against pace-of-life
            c: 0., // Locus determining effort intercept of observation effort against POL
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
            mu:0.2, // mutation rate of loci
            mut_size:0.2,

            b_f:1.0, // Baseline fecundity (offspring independence)
            b_s:0.99, // Baseline survival rate of adults (5 years lifespan)
            b_p:1.0, // Brood predation risk (assuming twice the mortality of adults)
            c_q:0.25, // survival cost of fast POL in winter 
            c_v:0.1, // surival cost of making observations in winter
            c_u:0.25, // survival cost of parental effort in winter (highest energy expenditure)
            sigma0, // Locus determining prior alpha
            mu_q:0.5,
            sigma_cue:4.0, // sd of the cue
            varsigma:10, // the maximum sample size an individual is able to gather from cues
            h:10., // slow-fast slope of nest-defence/size/aggression (behavioural phenotypes more extreme therefore h>theta)
            theta: 10.0, // slow-fast sigmoid slope of mortality risk against q-value (physiological
            div_rate:1.0, // divorce rate
        };

////////////////////////////////////////////////////// start r session
    let mut r = RSession::new()?;
    r.exec("source('src/b_s.r')")?;
    let r_mutex = Mutex::new(r);
////////////////////////////////////////////////////// hawkishness of fast individuals
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
        let x = rng.gen_range(0.0..50.0); // uniform sample from parameter space
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

////////////////////////////////////////////////////// mort (theta) of fast individuals
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
        let x = rng.gen_range(0.0..50.0); // uniform sample from parameter space
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

////////////////////////////////////////////////////// Adult mortality rate (b_s) simulations
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
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.b_s = x;
        
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

////////////////////////////////////////////////////// Cue cost sims
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

////////////////////////////////////////////////////// Cue detectability (cue_sigma) sims
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
        let x = rng.gen_range(0.001..10.0); // uniform sample from parameter space
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
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_sigmacue_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_sigmacue_plot('{}')", path)).unwrap();

////////////////////////////////////////////////////// Offspring independence (b_f) simulations
    // Construct the full path
    let path = format!("./Results/{}/b_f/", project_id);
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
        let x = rng.gen_range(0.0..5.0); // uniform sample from parameter space
        env.b_f = x;
        
        println!("Simulation started: baseline fecundity: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: baseline fecundity: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_b_f_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_b_f_plot('{}')", path)).unwrap();

////////////////////////////////////////////////////// Brood predation risk (b_p) sims
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

////////////////////////////////////////////////////// hawkishness of fast individuals
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
        let x = rng.gen_range(0.0..50.0); // uniform sample from parameter space
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

////////////////////////////////////////////////////// Divorce rate: 
    // Construct the full path
    let path = format!("./Results/{}/divorce/", project_id);
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
        env.div_rate = x;
        
        println!("Simulation started: divorce: {}, trial: {}", x, g);
        run(generations, &path, Some(&(x.to_string())), Some(&g),env);
        println!("Simulation done: divorce: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        if g % 10 == 0 {
            let mut r_guard = r_mutex.lock().unwrap(); 
            r_guard.exec(&format!("run_divorce_plot('{}')", path)).unwrap();
        }
    });
    let mut r = RSession::new()?;
    r.exec(&format!("run_divorce_plot('{}')", path)).unwrap();

    
Ok(())
}