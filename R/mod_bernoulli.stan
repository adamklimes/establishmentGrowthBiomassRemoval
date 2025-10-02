data{
  int N;
  int NspecTime;
  int NclonTime;
  array[N] int resp;
  array[N] int clonTime;
  array[N] int specTime;
}
parameters {
  vector[NclonTime] B;
  vector[NspecTime] randSpecTimeStd;
  real<lower=0> sigR;
}
transformed parameters {
  vector[NspecTime] randSpecTime;
  randSpecTime = randSpecTimeStd * sigR;
}
model {
  vector[N] aux;
//priors
  B ~ normal(0,1);
  randSpecTime ~ normal(0,1);
  sigR ~ normal(0,1);
//model
  randSpecTimeStd ~ std_normal();
  aux = B[clonTime] + randSpecTime[specTime];
  resp ~ bernoulli_logit(aux);
}
