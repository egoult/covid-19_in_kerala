// Implementation of the coupled covid model with stringecy as a covariate

//start_globs
// Expit transform for parameters constrained in interval [a,b]
static double expitCons(double x, double a, double b) {
	double out = (a + b * exp(x)) / (1.0 + exp(x));
	if(ISNAN(out)) out = (b + a * exp(-x)) / (1.0 + exp(-x)); // If x=+Inf, must return b
	return out;
}

// Logit transform for parameters constrained in interval [a,b]
static double logitCons(double x, double a, double b) {
	x = (x <= a) ? a : (x >= b ? b : x);
	double out = log((x - a) / (b - x));
	return out;
}

// Probability of transition for event of rate r during time step delta_t
// p = 1 - exp(-r * delta_t)
static double pTrans(double r, double delta_t) {
// r: event (instantaneous rate)
// delta_t: time step
// Returns: probability of transition
	double out = 1.0 - exp(-r * delta_t);
	return out;
}
//end_globs

//start_toest
int i; 
double *misc_vec = (double *) &misc1;
double *T_misc_vec = (double *) &T_misc1;

T_N = N;
T_lambda1 = log(lambda1);
T_lambda2 = log(lambda2);
T_lambda3 = log(lambda3);
T_sigma = log(sigma);
T_omega   = log(omega);
T_sensitivity = logit(sensitivity);
T_specificity = logit(specificity);
T_d = logit(d); 
T_p = log(p);
T_r = log(r);
T_kc = log(kc);
T_kd = lod(kd);

for(i = 0; i < 9; i++) T_misc_vec[i] = misc_vec[i];
//end_toest

//start_fromest
int i;
double *misc_vec = (double *) &misc1;
double *T_misc_vec = (double *) &T_misc1;

N = T_N;
lambda1 = exp(T_lambda1);
lambda2 = exp(T_lambda2);
lambda3 = exp(T_lambda3);
sigma = exp(T_sigma);
omega   = exp(T_omega);
sensitivity = exp(T_sensitivity);
specificity = exp(T_specificity);
d = exp(T_d); 
p = exp(T_p);
r = exp(T_r);
kc = exp(T_kc);
kd = exp(T_kd);

for(i = 0; i < 9; i++) misc_vec[i] = T_misc_vec[i];
//end_fromest

//start_dmeas
double fD, fC;
fD = ISNA(deaths) ? 0.0 : dnbinom_mu(nearbyint(deaths), 1.0 / kd, DM, 1);
fC = ISNA(cases) ? 0.0 : dnbinom_mu(nearbyint(cases), 1.0 / kc, CM, 1);
lik = (give_log) ? (fD) : exp(fD);
//end_dmeas

//start_rmeas
deaths = rnbinom_mu(1.0 / tau, DM);
//end_rmeas


//start_rinit
double *Q_vec = (double *) &Q1;
int i; 

S =  nearbyint(N);
E = 0;
I = 0;
R = 0;

SH = 0;
EH = 0;
IH = 0;
RH = 0;

CM = 0
DM = 0;
//end_rinit

//start_skel
double prev_t; // Prevalence of infection
double lambda_t; // Force of infection
double beta_t; // transmission rate
double sigma2, gamma2; // calculation saving 
double fstring0, fstring1;// stringency function
int i; 

//Vectors
double *Q_vec = (double *) &Q1;
double *DQ_vec = (double *) &DQ1;


// Calculate transmission rate and force of infection
prev_t = I / N + iota; // Prevalence of infection

lambda_t = lambda * prev_t; // Force of infection, natural scale


sigma2 = 2.0 * sigma; // calculated multiple times
gamma2 = 2.0 * gamma; // calculated multiple times

fstring0 = b1 *  log(1- stringency_0/100); // function pf stringency index at time t
fstring1 = b1 *  log(1- stringency_1/100); // function pf stringency index at time t-1


// Derivatives
DS  = - lambda_t * S + omega * SH + specificity * deltas - total_delta;
DE = lambda_t * S - p * E + (1-sensitivity) * deltae;
DI = p * E - r * I - sigma * I + (1-sensitivity) * deltai ;
DR = r * I + omegaw * RH - specificity * deltar;


DSH = (1 - specificity) * deltas - omega * SH;
DEH = sensitivity * deltae - p * EH;
DIH = p * EH + sigma * I + sensitivity * deltai - r*IH - d * IH;
DRH = r * IH - omegaw * RQ + (1 - specificity) * deltar;

DCM = ;
DDM = d * IH;
//end_skel

//start_rsim
double prev_t; // Prevalence of infection
double beta_t; // Transmission rate
double lambda_t; // Force of infection
double sigma2, gamma2; // calculation saving 
double fstring0, fstring1;// stringency function
int i; 

//Vectors
double *Q_vec = (double *) &Q1;

// Calculate transmission rate and force of infection
double dW = rnorm(0.0, beta_sd * sqrt(dt)); // Brownian motion, effective reproduction number
double edW = exp(dW);
prev_t = (I1 + I2) / N + iota; // Prevalence of infection
beta_t = Re * gamma;  //transmission rate
lambda_t = beta_t * prev_t ; // Force of infection

sigma2 = 2.0 * sigma; // calculated multiple times
gamma2 = 2.0 * gamma; // calculated multiple times

fstring0 = b1 *  log(1- stringency_0/100); // function pf stringency index at time t
fstring1 = b1 *  log(1- stringency_1/100); // function pf stringency index at time t-1


// Calculate no of transitions 
double fromS, fromE1, fromE2, fromI1, fromI2, toC2, toQ1;
double fromQ[n_Q];

fromS = rbinom(S, pTrans( lambda_t, dt)); // S->E1
fromE1 = rbinom(E1, pTrans(sigma2, dt)); // E1->E2
fromE2 = rbinom(E2, pTrans( sigma2, dt)); //E2->I1
fromI1 = rbinom(I1, pTrans( gamma2, dt)); //I1->I2
fromI2 = rbinom(I2, pTrans( gamma2, dt)); // I2->R

for(i = 0; i < n_Q; i++) {
	fromQ[i] = rbinom(Q_vec[i], pTrans(n_Q * kappa, dt)); // Q_i->Q_{i+1}
}
toQ1 = rbinom(I1, pTrans(gamma2 * mu, dt)); //I1->Q1

// Balance equations
S += -fromS;
E1 += fromS - fromE1;
E2 += fromE1 - fromE2;
I1 += fromE2 - fromI1; 
I2 += fromI1 - fromI2; 
R += fromI2; 
Re *= exp( fstring0 ) * edW / exp( fstring1 ); // should use this one
if(debug) {
	Rprintf("t=%.2f, fstring = %.3f\n", t, fstring0);
}
// Re = (Re - fstring1) * edW + fstring0 > 0 ? (Re - fstring1) * edW + fstring0 : 0; 
//Re += (edW-1) * Re - fstring1 * edW + fstring0;
Q_vec[0] += toQ1 - fromQ[0];
for(i = 1; i < n_Q; i++) {
	Q_vec[i] += fromQ[i - 1] - fromQ[i]; 
} 
DM += fromQ[n_Q - 1];
//end_rsim
