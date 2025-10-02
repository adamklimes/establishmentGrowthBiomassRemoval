data{
  int N;
  int Nspec;
  vector[N] resp;
  array[N] int clon;
  array[N] int spec;
}
parameters {
  vector [2] B;
  vector[Nspec] randSpecStd;
  real<lower=0> sig;
  real<lower=0> sigR;
}
transformed parameters {
  vector[Nspec] randSpec;
  randSpec = randSpecStd * sigR;
}
model {
  vector[N] aux;
//priors
  B ~ normal(0,1);
  randSpec ~ normal(0,1);
  sig ~ normal(0,1);
  sigR ~ normal(0,1);
//model
  randSpecStd ~ std_normal();
  aux = B[clon] + randSpec[spec];
  resp ~ normal(aux, sig);
}
