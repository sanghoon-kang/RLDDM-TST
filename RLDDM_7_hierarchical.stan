functions{
  //random number generator
  vector wiener_rng(real a, real tau, real z, real d){
    vector[2] simul;
    real dt;
    real sigma;
    real p;
    real y;
    real i;
    real aa;
    int ch;
    real rt;
    dt=.0001;
    sigma=1;
    y=z*a;
    p=.5*(1+((d*sqrt(dt))/sigma));
    i=0;
    while (y<a && y>0){
      aa=uniform_rng(0,1);
      if (aa<=p){
        y=y+sigma*sqrt(dt);
        i=i+1;
      }else{
        y=y-sigma*sqrt(dt);
        i=i+1;
      }
    }
    ch=(y<=0)*1+1; //1 for left 2 for right
    rt=i*dt+tau;
    
    simul[1]=ch;
    simul[2]=rt;
    return simul;
  }
}

data {
  int<lower=1> N;  // Number of subjects
  int<lower=1> Tsubj[N];  // Number of trials for each subject
  int<lower=1> T; // Total number of trials
  int<lower=0, upper=1> transition[N,T];    // rare or common transition 
  int<lower=2, upper=3> secstagestate[N,T]; // second stage state
  int<lower=1, upper=2> choice1[N,T];
  int<lower=1, upper=2> choice2[N,T];
  real RT1[N,T];   // response times from first stage
  real RT2[N,T];   // response times from second stage
  real minRT[N];   // minimum RT for stage 1
  real RTbound;    // lower bound or RT (e.g., 0.1 second) for stage 1
  int<lower=0> reward[N,T];     // reward
}

parameters {
  vector[N] alpha_p; // first stage learning rate
  vector[N] p_p; //perseveration parameter
  vector[N] omega_p; // model-based weight
  vector[N] lambda_p; // decay parameter
  vector[N] a_p;  // first stage boundary separation
  vector[N] b_p; //scaling parameter for first stage
  vector[N] tau_p;  // nondecision time for first stage
  vector[7] mu_p;
  vector<lower=0>[7] sigma;
}

transformed parameters{
  vector<lower=0, upper=1>[N] alpha;
  vector[N] p;
  vector<lower=0, upper=1>[N] omega;
  vector<lower=0, upper=1>[N] lambda;
  vector<lower=0>[N] a;
  vector[N] b;
  vector<lower=RTbound,upper=max(minRT)>[N] tau;
  
  for (i in 1:N){
    alpha[i] = Phi_approx(mu_p[1]+sigma[1]*alpha_p[i]);
    p[i] = Phi_approx(mu_p[2] + sigma[2] * p_p[i]);
    omega[i] = Phi_approx(mu_p[3]+sigma[3]*omega_p[i]);
    lambda[i] = Phi_approx(mu_p[4]+sigma[4]*lambda_p[i]);
    tau[i] = Phi_approx(mu_p[7]+sigma[7]*tau_p[i])*(minRT[i]-RTbound) + RTbound;
  }
  a = exp(mu_p[5]+sigma[5]*a_p);
  b = mu_p[6] + sigma[6] * b_p;
} 

model {
  mu_p ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  
  alpha_p ~ normal(0.5,1);
  p_p ~ normal(0.5, 1);
  a_p ~ normal(0, 1);
  b_p ~ normal(0, 1);
  tau_p ~ normal(0, 1);
  omega_p ~ normal(0, 1);
  lambda_p ~ normal(0, 1);
  
  for (i in 1:N){
    real Qs[3,2];
    real per[2];
    real Qmb[2];
    real Qnet[2];
    real delta2;
    real PE1;
    real PE2;
    for (e in 1:3){
      for (r in 1:2){
        Qs[e,r] = 0;
      }
    }
    for (r in 1:2){
      Qmb[r] = 0;
      Qnet[r] = 0;
    }
    delta2 = 0;

    per[1] = 0;
    per[2] = 0;
    for (t in 1:Tsubj[i]){
      //print(Qs)
      //first stage choice
      Qmb[1]=.7*max(Qs[2,1:2])+.3*max(Qs[3,1:2]); //should we specify initial value?
      Qmb[2]=.3*max(Qs[2,1:2])+.7*max(Qs[3,1:2]);
      Qnet[1]=(1-omega[i])*Qs[1,1] + omega[i]*Qmb[1] + p[i]*per[1];
      Qnet[2]=(1-omega[i])*Qs[1,2] + omega[i]*Qmb[2] + p[i]*per[2];
      per[1] = 0;
      per[2] = 0;
      //choice
      if (choice1[i,t] == 1){
        RT1[i,t] ~ wiener(a[i], tau[i], 0.5, b[i]*(Qnet[1]-Qnet[2]));
        if (secstagestate[i,t] == 2){
          if(choice2[i,t] == 1){
            delta2=b[i]*(Qs[2,1]-Qs[2,2]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[2,1]-Qs[1,1]);
            PE2=(reward[i,t]-Qs[2,1]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,1]=Qs[2,1]+alpha[i]*PE2;
            per[1] = 1;
          }else{
            delta2=b[i]*(Qs[2,2]-Qs[2,1]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[2,2]-Qs[1,1]);
            PE2=(reward[i,t]-Qs[2,2]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,2]=Qs[2,2]+alpha[i]*PE2;
            per[1] = 1;
          }
        }else{
          if(choice2[i,t] == 1){
            delta2=b[i]*(Qs[3,1]-Qs[3,2]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[3,1]-Qs[1,1]);
            PE2=(reward[i,t]-Qs[3,1]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,1]=Qs[3,1]+alpha[i]*PE2;
            per[1] = 1;
          }else{
            delta2=b[i]*(Qs[3,2]-Qs[3,1]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[3,2]-Qs[1,1]);
            PE2=(reward[i,t]-Qs[3,2]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,2]=Qs[3,2]+alpha[i]*PE2;
            per[1] = 1;
          }
        }
      }else{
        RT1[i,t] ~ wiener(a[i], tau[i], 0.5, b[i]*(Qnet[2]-Qnet[1]));
        if (secstagestate[i,t] == 2){
          if(choice2[i,t] == 1){
            delta2=b[i]*(Qs[2,1]-Qs[2,2]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[2,1]-Qs[1,2]);
            PE2=(reward[i,t]-Qs[2,1]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,1]=Qs[2,1]+alpha[i]*PE2;
            per[2] = 1;
          }else{
            delta2=b[i]*(Qs[2,2]-Qs[2,1]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[2,2]-Qs[1,2]);
            PE2=(reward[i,t]-Qs[2,2]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,2]=Qs[2,2]+alpha[i]*PE2;
            per[2] = 1;
          }
        }else{
          if(choice2[i,t] == 1){
            delta2=b[i]*(Qs[3,1]-Qs[3,2]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[3,1]-Qs[1,2]);
            PE2=(reward[i,t]-Qs[3,1]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,1]=Qs[3,1]+alpha[i]*PE2;
            per[2] = 1;
          }else{
            delta2=b[i]*(Qs[3,2]-Qs[3,1]);
            RT2[i,t] ~ wiener(a[i], tau[i], 0.5, delta2);
            PE1=(Qs[3,2]-Qs[1,2]);
            PE2=(reward[i,t]-Qs[3,2]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,2]=Qs[3,2]+alpha[i]*PE2;
            per[2] = 1;
          }
        }
      }
    }
  }
}

generated quantities {
  real mu_alpha;
  real mu_pp;
  real mu_omega;
  real mu_lambda;
  real mu_a;
  real mu_b;
  real mu_tau;
  
    
  real trial[N,T];
  real trn[N,T];
  real lower_boundary;
  real upper_boundary;
  real rwrdprob[N,T,4]; //reward probability
  real rwrd[N,T,4]; //number of subjects, trials and fractals
  real arwrd[N,T];
  real asecstagestate[N,T];
  real rt1[N,T];
  real rt2[N,T];
  real ch1[N,T];
  real ch2[N,T];
  real rand;
  
  mu_alpha = Phi_approx(mu_p[1]);
  mu_pp = mu_p[2];
  mu_omega = Phi_approx(mu_p[3]);
  mu_lambda = Phi_approx(mu_p[4]);
  mu_a = exp(mu_p[5]);
  mu_b = mu_p[6];
  mu_tau = Phi_approx(mu_p[7]) * (mean(minRT)-RTbound) + RTbound;
  
  //generate reward probabilities (according to gaussian random walk)
  for (i in 1:N){
    lower_boundary = 0.2;
    upper_boundary = 0.8;
    for (f in 1:4){
      rwrdprob[i,1,f]=uniform_rng(lower_boundary,upper_boundary); //initial probability
      for (t in 2:T){
        for(h in 1:4){
          rwrdprob[i,t,h]=rwrdprob[i,t-1,h]+0.025*normal_rng(0,1);
          if (rwrdprob[i,t,h]>upper_boundary){
            rwrdprob[i,t,h]=upper_boundary;
          }else if (rwrdprob[i,t,h]<lower_boundary){
            rwrdprob[i,t,h]=lower_boundary;
          }
        }
      }
    }
  }
  //generate rewards
  for (i in 1:N){
    for (t in 1:T){
      for (f in 1:4){
        rand = uniform_rng(0,1);
        if (rand<=rwrdprob[i,t,f]){
          rwrd[i,t,f]=1;
        }else if (rand>rwrdprob[i,t,f]){
          rwrd[i,t,f]=0;
        }
      }
    }
  }
  
  
  //simulate
  for (i in 1:N){
    real Qs[3,2];
    real per[2];
    real Qmb[2];
    real Qnet[2];
    real delta;
    real delta2;
    real PE1;
    real PE2;
    for (e in 1:3){
      for (r in 1:2){
        Qs[e,r] = 0;
      }
    }
    per[1] = 0;
    per[2] = 0;
    for (t in 1:T){
      trial[i,t] = t;
      //print(Qs);
      //first stage choice
      Qmb[1]=.7*max(Qs[2,1:2])+.3*max(Qs[3,1:2]);
      Qmb[2]=.3*max(Qs[2,1:2])+.7*max(Qs[3,1:2]);
      Qnet[1]=(1-omega[i])*Qs[1,1] + omega[i]*Qmb[1]+p[i]*per[1];
      Qnet[2]=(1-omega[i])*Qs[1,2] + omega[i]*Qmb[2]+p[i]*per[2];
      per[1] = 0;
      per[2] = 0;
      ch1[i,t] = wiener_rng(a[i], tau[i], 0.5, b[i]*(Qnet[1]-Qnet[2]))[1]; //1 for left and 2 for right
      trn[i,t]=(uniform_rng(0,1)>=.7)*1; // 1 for rare, 0 for common
      if (ch1[i,t] == 1){
        rt1[i,t] = wiener_rng(a[i], tau[i], 0.5, b[i]*(Qnet[1]-Qnet[2]))[2];
        if (trn[i,t] == 0){
          asecstagestate[i,t] = 1;
          ch2[i,t] = wiener_rng(a[i], tau[i], 0.5, b[i]*(Qs[2,1]-Qs[2,2]))[1];
          if(ch2[i,t] == 1){
            delta2=b[i]*(Qs[2,1]-Qs[2,2]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,1];
            PE1=(Qs[2,1]-Qs[1,1]);
            PE2=(rwrd[i,t,1]-Qs[2,1]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,1]=Qs[2,1]+alpha[i]*PE2;
            per[1] = 1;
          }else{
            delta2=b[i]*(Qs[2,2]-Qs[2,1]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,2];
            PE1=(Qs[2,2]-Qs[1,1]);
            PE2=(rwrd[i,t,2]-Qs[2,2]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,2]=Qs[2,2]+alpha[i]*PE2;
            per[1] = 1;
          }
        }else if (trn[i,t] == 1){
          asecstagestate[i,t] = 2;
          ch2[i,t] = wiener_rng(a[i], tau[i], 0.5, b[i]*(Qs[3,1]-Qs[3,2]))[1];
          if(ch2[i,t] == 1){
            delta2=b[i]*(Qs[3,1]-Qs[3,2]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,3];
            PE1=(Qs[3,1]-Qs[1,1]);
            PE2=(rwrd[i,t,3]-Qs[3,1]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,1]=Qs[3,1]+alpha[i]*PE2;
            per[1] = 1;
          }else{
            delta2=b[i]*(Qs[3,2]-Qs[3,1]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,4];
            PE1=(Qs[3,2]-Qs[1,1]);
            PE2=(rwrd[i,t,4]-Qs[3,2]);
            Qs[1,1]=Qs[1,1]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,2]=Qs[3,2]+alpha[i]*PE2;
            per[1] = 1;
          }
        }
      }else if (ch1[i,t] == 2){
        rt1[i,t] = wiener_rng(a[i], tau[i], 1-0.5, -b[i]*(Qnet[1]-Qnet[2]))[2];
        if (trn[i,t] == 1){
          asecstagestate[i,t] = 1;
          ch2[i,t] = wiener_rng(a[i], tau[i], 0.5, b[i]*(Qs[2,1]-Qs[2,2]))[1];
          if(ch2[i,t] == 1){
            delta2=b[i]*(Qs[2,1]-Qs[2,2]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,1];
            PE1=(Qs[2,1]-Qs[1,2]);
            PE2=(rwrd[i,t,1]-Qs[2,1]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,1]=Qs[2,1]+alpha[i]*PE2;
            per[2] = 1;
          }else{
            delta2=b[i]*(Qs[2,2]-Qs[2,1]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,2];
            PE1=(Qs[2,2]-Qs[1,2]);
            PE2=(rwrd[i,t,2]-Qs[2,2]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[2,2]=Qs[2,2]+alpha[i]*PE2;
            per[2] = 1;
          }
        }else if (trn[i,t] == 0){
          asecstagestate[i,t] = 2;
          ch2[i,t] = wiener_rng(a[i], tau[i], 0.5, b[i]*(Qs[3,1]-Qs[3,2]))[1];
          if(ch2[i,t] == 1){
            delta2=b[i]*(Qs[3,1]-Qs[3,2]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,3];
            PE1=(Qs[3,1]-Qs[1,2]);
            PE2=(rwrd[i,t,3]-Qs[3,1]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,1]=Qs[3,1]+alpha[i]*PE2;
            per[2] = 1;
          }else{
            delta2=b[i]*(Qs[3,2]-Qs[3,1]);
            rt2[i,t] = wiener_rng(a[i], tau[i], 0.5, delta2)[2];
            arwrd[i,t] = rwrd[i,t,4];
            PE1=(Qs[3,2]-Qs[1,2]);
            PE2=(rwrd[i,t,4]-Qs[3,2]);
            Qs[1,2]=Qs[1,2]+alpha[i]*PE1+alpha[i]*lambda[i]*PE2;
            Qs[3,2]=Qs[3,2]+alpha[i]*PE2;
            per[2] = 1;
          }
        }
      }
    }
  }
}
