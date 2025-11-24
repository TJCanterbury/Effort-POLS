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
        "b_s", // Offspring survival benefit of parental effort
        "s_p", // baseline survival rate of parents in winter
        "c_q", // survival cost of fast POL in winter
        "c_v", // surival cost of making observations in winter
        "c_u", // survival cost of parental effort in winter
        "mu", // mutation rate of loci
        "sigma0", // Locus determining prior alpha
        "sigma_cue",
        "div_rate",
        "varsigma", // the maximum sample size an individual is able to gather from cues
        "mean_fitness",
        "m",
        "c",
        "u_base",
        "rho",
        "nu",
        "gamma",
        "lambda",
        "x"
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
        "b_s", // Offspring survival benefit of parental effort
        "s_p", // baseline survival rate of parents in winter
        "c_q", // survival cost of fast POL in winter
        "c_v", // surival cost of making observations in winter
        "c_u", // survival cost of parental effort in winter
        "mu", // mutation rate of loci
        "sigma0", // Locus determining prior alpha
        "sigma_cue",
        "div_rate",
        "varsigma", // the maximum sample size an individual is able to gather from cues
        "mean_fitness",
        "m",
        "c",
        "u_base",
        "rho",
        "nu",
        "gamma",
        "lambda",
        "x"
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

fn init_pop(n: u32, agent: Agent) -> Vec<Agent> {
    // returns vector with new agents
    let mut pop = Vec::new();
    for _i in 0..n {
        pop.push(agent.clone());
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
    x:f64, // Expected POL
}
#[derive(Clone, Debug)]
struct Environment {
    pop: Vec<Agent>,
    singles: Vec<usize>,
    dead: Vec<usize>,
    mean_fitness:f32,
    tol:f64,
    b_f:f64, // Benefit of additive parental effort to fecundity
    b_s:f64, // Offspring survival benefit of parental effort
    s_p:f64, // baseline survival rate of parents in winter
    c_q:f64, // survival cost of fast POL in winter
    c_v:f64, // surival cost of making observations in winter
    c_u:f64, // survival cost of parental effort in winter
    mu:f64, // mutation rate of loci
    sigma0: f64, // Prior sigma
    sigma_cue:f64, // variation of social cue
    varsigma:f64, // the maximum sample size an individual is able to gather from cues
    theta:f64, // slope of mortality against q-value
    mut_size:f64, // standard deviation of mutations
    div_rate:f64, // divorce rate
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

        dice = rng.gen(); //x
        if dice <= mu { //mutates
            // x
            let n = nm::new(self.gamma, mut_size as f64).unwrap(); // generates normal distribution 
            let loc_mut = n.sample(&mut rng);
            self.x = loc_mut;
        }
    }

    fn pol_v(&self) -> f64 {
        return self.c / (1. + (self.m * (self.q) ).exp());
    }

    fn social_cue(&mut self, partner_q:f64, sigma_cue:f64, n:i64) {
        let normal = Normal::new(partner_q, sigma_cue).unwrap();
        let mut phen:f64;
        let b = 1. / (sigma_cue * sigma_cue);

        // println!("n: {}, prior pi: {}, sig: {}, phen: {}, q: {}",i, self.pi, self.sigma, phen, partner_q);
        for _i in 0..n {
            let a = 1. / (self.sigma * self.sigma);
            phen = normal.sample(&mut thread_rng());
            self.updating(a, b, phen);
            // println!("n: {}, post pi: {}, sig: {}, phen: {}, q: {}\n",i, self.pi, self.sigma, phen, partner_q);

        }
    }

    fn updating(&mut self, a:f64, b:f64, x:f64) {
        self.pi = (a * self.pi + b * x) / (a + b);
        self.sigma = (1. / (a + b)).sqrt();
    }

    fn r(&self, u2:f64) -> f64 {
        let r:f64 = self.u_base + self.rho + self.nu * self.q + self.gamma * self.pi - self.lambda * (self.u_base  - u2);
        let zero:f64 = 0.;
        return r.max(zero).min(1.0)
    }

    fn surivorship(&self, s_p:f64, c_q:f64, c_v:f64, c_u:f64, theta:f64) -> f64 {
        let s: f64 =  s_p*(1. - (1. / (1.+(-self.q * theta).exp())) * c_q)*(1. - self.pol_v() * c_v)*(1. - self.u * c_u);
        // println!("s: {}, c: {}, u: {}, q..: {}",s, (1. - self.pol_v() * c_v), (1. - self.u * c_u), (1. - (1. / (1.+(-self.q * theta).exp())) * c_q));
        return s
    }
}

impl Environment { 
    fn development(&mut self) {
        let mut phen:f64;
        for i in 0..self.pop.len() {
            let normal = Normal::new(self.pop[i].x, self.sigma0).unwrap();
            phen = normal.sample(&mut thread_rng());
            self.pop[i].q = phen;
        }
    }
    
    fn observations(&mut self) {
        let mut v_f: f64;
        let mut n_f: i64;
        let mut q_m: f64;
        for i in 0..self.pop.len(){
            let female = i;
            let male = self.pop[i].partner.unwrap();
            v_f = self.pop[female].pol_v();

            n_f = (v_f * self.varsigma) as i64;

            q_m = self.pop[male].q;
            
            self.pop[female].social_cue(q_m, self.sigma_cue, n_f);
        } 
    }

    fn negotiations(&mut self) {  
        // determine male and female parental effort and then update fitness
        // let mut dif:f64;
        let mut u1:f64;
        let mut u2:f64;
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
                u1=0.;
                u2=0.;
                // old=u1+u2;

                // negotiations:
                // while dif > self.tol {
                for _j in 0..100 {
                    // Observation preference
                    // gen+=1.;
                    u1 = self.pop[male].r(u2);
                    u2 = self.pop[female].r(u1);
                    // dif = (u1 + u2)-old;
                    // old = u1 + u2;
                    // println!("gen: {}, dif: {}, old: {}, u1: {}, u2: {}", gen, dif, old, u1, u2)
                }
                
                self.pop[male].u = u1;
                self.pop[female].u = u2;

                ben = self.b_f * (u1 + u2);
                q1 = self.pop[male].q;
                q2 = self.pop[female].q;
                pred = self.b_s * (u1 * (1. / (1.+(-q1 * self.theta).exp())) + u2 * (1. / (1.+(-q2 * self.theta).exp())));
                
                self.pop[male].fitness = (0. as f64).max(ben * pred / 2.);
                self.pop[female].fitness = (0. as f64).max(ben * pred / 2.);
            }
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
            if rng.gen::<f64>() > self.pop[i].surivorship(self.s_p, self.c_q, self.c_v, self.c_u, self.theta) {
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
                self.pop[*a].pi = self.pop[*a].x;
                self.pop[*a].partner = Some(*b);

                self.pop[*b].sigma = self.sigma0;
                self.pop[*b].pi = self.pop[*b].x;
                self.pop[*b].partner = Some(*a);
            }
        }
    }

    fn natural_selection(&mut self) {
        let mut selected_male:usize;
        let mut selected_female:usize;
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

            // println!("{}", self.pop[selected].fitness);
            // inherit policies
            
            // let parent = if rng.gen::<f64>() < 0.5 { selected_male } else { selected_female };
            self.pop[self.dead[i]].m = (self.pop[selected_male].m + self.pop[selected_female].m) / 2. ;
            self.pop[self.dead[i]].c = (self.pop[selected_male].c + self.pop[selected_female].c) / 2. ;
            self.pop[self.dead[i]].u_base = (self.pop[selected_male].u_base + self.pop[selected_female].u_base) / 2. ;
            self.pop[self.dead[i]].rho = (self.pop[selected_male].rho + self.pop[selected_female].rho) / 2. ;
            self.pop[self.dead[i]].nu = (self.pop[selected_male].nu + self.pop[selected_female].nu) / 2. ;
            self.pop[self.dead[i]].lambda = (self.pop[selected_male].lambda + self.pop[selected_female].lambda) / 2. ;
            self.pop[self.dead[i]].gamma = (self.pop[selected_male].gamma + self.pop[selected_female].gamma) / 2. ;
            self.pop[self.dead[i]].x = (self.pop[selected_male].x + self.pop[selected_female].x) / 2. ;
            // self.pop[self.dead[i]].m = self.pop[parent].m;
            // self.pop[self.dead[i]].c = self.pop[parent].c;
            // self.pop[self.dead[i]].u_base = self.pop[parent].u_base;
            // self.pop[self.dead[i]].rho = self.pop[parent].rho;
            // self.pop[self.dead[i]].nu = self.pop[parent].nu;
            // self.pop[self.dead[i]].lambda = self.pop[parent].lambda;
            // self.pop[self.dead[i]].gamma = self.pop[parent].gamma;
            // self.pop[self.dead[i]].x = self.pop[parent].x;
            self.pop[self.dead[i]].fitness = 0.;
            // self.pop[self.dead[i]].partner = None;

            // apply mutations to policies
            self.pop[self.dead[i]].mutate(self.mu, self.mut_size);   
        }
        self.widowed();
        self.dead.clear();
        // // move dead into singles
        // while self.dead.len() > 0 {
        //     let index = self.dead[self.dead.len() - 1];
        //     self.pop[index].fitness=0.;
        //     self.pop[index].sigma=self.sigma0;
        //     self.pop[index].pi=0.;

        //     // add developmental noise

        //     self.singles.push(self.dead.pop().unwrap());
        // }
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
        let mut mean_loc_x:f64=0.;
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
            mean_loc_x += a.x;
        }
        mean_loc_m = mean_loc_m/pop_len;
        mean_loc_c = mean_loc_c/pop_len;
        mean_loc_u_base = mean_loc_u_base/pop_len;
        mean_loc_rho = mean_loc_rho/pop_len;
        mean_loc_nu = mean_loc_nu/pop_len;
        mean_loc_lambda = mean_loc_lambda/pop_len;
        mean_loc_gamma = mean_loc_gamma/pop_len;
        mean_loc_x = mean_loc_x/pop_len;
            
        data.push(i);
        data.push(pop_len);
        data.push(self.b_f); // Benefit of additive parental effort to fecundity
        data.push(self.b_s); // Offspring survival benefit of parental effort
        data.push(self.s_p); // baseline survival rate of parents in winter
        data.push(self.c_q); // survival cost of fast POL in winter
        data.push(self.c_v); // surival cost of making observations in winter
        data.push(self.c_u); // survival cost of parental effort in winter
        data.push(self.mu); // mutation rate of loci
        data.push(self.sigma0); // Locus determining prior alpha
        data.push(self.sigma_cue);
        data.push(self.div_rate);
        data.push(self.varsigma); // the maximum sample size an individual is able to gather from cues
        data.push(self.mean_fitness as f64);
        data.push(mean_loc_m);
        data.push(mean_loc_c);
        data.push(mean_loc_u_base);
        data.push(mean_loc_rho);
        data.push(mean_loc_nu);
        data.push(mean_loc_gamma);
        data.push(mean_loc_lambda);
        data.push(mean_loc_x);
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
        env.development(); // Individuals develop a random pace-of=life phenotype q
        env.pairing_pool(); // all individuals pair up randomly
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
        let generations: i32 = 30000; // number of generations the sim lasts
        let sigma0:f64 = 6.0;
        let iterations:i32 = 100;

        let agent:Agent = Agent {
            sigma: sigma0,
            pi: 0.,
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
            x:0.,
        };
        
        let env = Environment {
            pop: init_pop(1000, agent),
            singles: (0..1000 as usize).collect(),
            dead: Vec::new(),
            mean_fitness:0.,
            tol: 0.001,
            b_f:2.0, // Benefit of additive parental effort to fecundity
            b_s:0.8, // Offspring survival benefit of parental effort
            s_p:0.8, // baseline survival rate of parents in winter
            c_q:0.25, // survival cost of fast POL in winter
            c_v:0.25, // surival cost of making observations in winter
            c_u:0.25, // survival cost of parental effort in winter
            mu:1.0, // mutation rate of loci
            sigma0, // Locus determining prior alpha
            sigma_cue:6.0, // sd of the cue
            varsigma:10.0, // the maximum sample size an individual is able to gather from cues
            theta: 2.0, // sigmoid slope of mortality risk against q-value
            mut_size:0.25,
            div_rate:1.0, // divorce rate
        };

////////////////////////////////////////////////////// start r session
    let mut r = RSession::new()?;
    r.exec("source('src/b_s.r')")?;
    let r_mutex = Mutex::new(r);

////////////////////////////////////////////////////// Run b_f simulations
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
        env.pop = init_pop(1000, agent);
        let x = rng.gen_range(0.0..10.0); // uniform sample from parameter space
        env.b_f = x;
        
        println!("Simulation started: benefit f: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: benefit f: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        let mut r_guard = r_mutex.lock().unwrap(); 
        r_guard.exec(&format!("run_b_f_plot('{}')", path)).unwrap();
    });

////////////////////////////////////////////////////// Run b_s simulations
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
        env.pop = init_pop(1000, agent);
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.b_s = x;
        
        println!("Simulation started: benefit s: {}, trial: {}", x, g);
        run(
            generations, 
            &path, 
            Some(&(x.to_string())), 
            Some(&g),
            env
        );
        println!("Simulation done: benefit s: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        let mut r_guard = r_mutex.lock().unwrap(); 
        r_guard.exec(&format!("run_b_s_plot('{}')", path)).unwrap();
    });

////////////////////////////////////////////////////// Run cue cost sims
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
        env.pop = init_pop(1000, agent);
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
        let mut r_guard = r_mutex.lock().unwrap(); 
        r_guard.exec(&format!("run_cuecost_plot('{}')", path)).unwrap();
    });

////////////////////////////////////////////////////// divorce: 
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
        env.pop = init_pop(1000, agent);
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.div_rate = x;
        
        println!("Simulation started: divorce: {}, trial: {}", x, g);
        run(generations, &path, Some(&(x.to_string())), Some(&g),env);
        println!("Simulation done: divorce: {}, trial: {}", x, g);

        // ensure only one thread runs r at a time (plotting):
        let mut r_guard = r_mutex.lock().unwrap(); 
        r_guard.exec(&format!("run_divorce_plot('{}')", path)).unwrap();

    });

////////////////////////////////////////////////////// Run mort sims
    // Construct the full path
        let path = format!("./Results/{}/mort/", project_id);
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
        env.pop = init_pop(1000, agent);
        let x = rng.gen_range(0.0..1.0); // uniform sample from parameter space
        env.s_p = x;
        
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
        let mut r_guard = r_mutex.lock().unwrap(); 
        r_guard.exec(&format!("run_mort_plot('{}')", path)).unwrap();
    });

////////////////////////////////////////////////////// Run POL sigma sims
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
        env.pop = init_pop(1000, agent);
        let x = rng.gen_range(0.01..10.0); // uniform sample from parameter space
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
        let mut r_guard = r_mutex.lock().unwrap(); 
        r_guard.exec(&format!("run_sigma_plot('{}')", path)).unwrap();
    });
        
////////////////////////////////////////////////////// Run sigma_cue sims
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
        env.pop = init_pop(1000, agent);
        let x = rng.gen_range(0.01..10.0); // uniform sample from parameter space
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
        let mut r_guard = r_mutex.lock().unwrap(); 
        r_guard.exec(&format!("run_sigmacue_plot('{}')", path)).unwrap();
    });



        

    
Ok(())
}